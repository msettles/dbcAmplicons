#!/usr/bin/env python
"""
Tag reads from dual barcoded experiments
"""

import os
import sys
import gzip
import time
#import re
from collections import Counter
from Bio import SeqIO

from optparse import OptionParser  # http://docs.python.org/library/optparse.html

# ---------------- Set up option parser and input ----------------
usage = "usage: %prog [options] input_prefix"
parser = OptionParser(usage=usage)
parser.add_option('-b', '--barcodediff', help="max hamming dist from barcode [default: %default]",
                  type="int", dest="barcodediff", default=1)
parser.add_option('-p', '--primerdiff', help="max hamming dist from primer [default: %default]",
                  type="int", dest="primerdiff", default=4)
parser.add_option('-e', '--primerend', help="required number of matching bases at end of primer [default: %default]",
                  type="int", dest="primerend", default=1)
parser.add_option('-u', '--uncompressed', help="leave output files uncompressed [default: %default]",
                  action="store_true", dest="uncompressed", default=False)
parser.add_option('-o', '--output_prefix', help="output file basename",
                  action="store", type="str", dest="output_base",default=None)
parser.add_option('-v', '--verbose', help="verbose output [default: %default]",
                  action="store_true", dest="verbose", default=False)

(options,  args) = parser.parse_args()  # uncomment this line for command line support

barcodeFile = "MetaData/barcodeLookupTable.txt"
primerForward = "MetaData/fwd_primer.fasta"
primerReverse = "MetaData/rev_primer.fasta"

if len(args) == 1:
    R1_File = args[0] + '_R1_001.fastq.gz'
    R2_File = args[0] + '_R2_001.fastq.gz'
    R3_File = args[0] + '_R3_001.fastq.gz'
    R4_File = args[0] + '_R4_001.fastq.gz'
    #Check input/output files:
    if not os.path.exists(R1_File):
        print "Error, can't find input file %s" % R1_File
        sys.exit()
    if not os.path.exists(R2_File):
        print "Error, can't find input file %s" % R2_File
        sys.exit()
    if not os.path.exists(R3_File):
        print "Error, can't find input file %s" % R3_File
        sys.exit()
    if not os.path.exists(R4_File):
        print "Error, can't find input file %s" % R4_File
        sys.exit()
else:
    print "Input file not specified on command line"
    sys.exit()

# ----------------------- check output files -----------------------
if options.output_base is None:
    Output_prefix = args[0]
else:
    Output_prefix = options.output_base

# ----------------------- open input files -----------------------
if R1_File.split(".")[-1] == "gz":
    R1 = SeqIO.parse(gzip.open(R1_File, 'rb'), 'fastq')
else:
    R1 = SeqIO.parse(open(R1_File, 'r'), 'fastq')

if R2_File.split(".")[-1] == "gz":
    R2 = SeqIO.parse(gzip.open(R2_File, 'rb'), 'fastq')
else:
    R2 = SeqIO.parse(open(R2_File, 'r'), 'fastq')

if R3_File.split(".")[-1] == "gz":
    R3 = SeqIO.parse(gzip.open(R3_File, 'rb'), 'fastq')
else:
    R3 = SeqIO.parse(open(R3_File, 'r'), 'fastq')

if R4_File.split(".")[-1] == "gz":
    R4 = SeqIO.parse(gzip.open(R4_File, 'rb'), 'fastq')
else:
    R4 = SeqIO.parse(open(R4_File, 'r'), 'fastq')

# ----------------------- open output files ------------
if options.uncompressed is True:
    outf = {'identified':[open(Output_prefix + '_R1.fastq', 'w'), open(Output_prefix + '_R2.fastq', 'w')], 'unidentified':[open(Output_prefix + '_Unidentified_R1.fastq', 'w'), open(Output_prefix + '_Unidentified_R2.fastq', 'w')]}
else:
    outf = {'identified':[gzip.open(Output_prefix + '_R1.fastq.gz', 'wb'), gzip.open(Output_prefix + '_R2.fastq.gz', 'wb')], 'unidentified':[gzip.open(Output_prefix + '_Unidentified_R1.fastq.gz', 'wb'), gzip.open(Output_prefix + '_Unidentified_R2.fastq.gz', 'wb')]}

barcodesFile = open(Output_prefix + '_report.txt', 'w')
uniquesFile = open(Output_prefix + '_uniques.txt', 'w')

# ----------------------- barcode/primer files ------------
barcodeMaxDiff = options.barcodediff
primerMaxDiff = options.primerdiff
primerEndMatch = options.primerend

# ---------------- read in barcodes and make a dictionary for lookup ----------------
barcodes = {}
P5 = []
P7 = []
IDS = []
bcfile = open(barcodeFile,'rb')
f = bcfile.readlines()[1:] # skip the first header line
for row in f:
    row = row.rstrip()
    ID, P5Name, P5BC, P7Name, P7BC = row.split('\t')
    if P5BC not in P5:
        P5.extend([P5BC])
    if P7BC not in P7:
        P7.extend([P7BC])
    IDS.extend([ID])
    barcodes["%s %s" % (P7BC, P5BC)] = ID

bcfile.close()

# ------- build primer dictionary --------
primersP5 = {}
primersP7 = {}

for seq in SeqIO.parse(open(primerForward, 'r'), 'fasta'):
    primersP5[seq.seq.tostring()] = seq.id

for seq in SeqIO.parse(open(primerReverse, 'r'), 'fasta'):
    primersP7[seq.seq.tostring()] = seq.id

# ------- make some counters for reporting ------------
counters = {}
for bc in barcodes:
    counters[barcodes[bc]] = [0, 0, 0]  # first value is perfect matches, second is 1bp mismatch in bc1, third is 1bp mismatch in bc2
goodReadsCounter = 0
otherCounter = 0
uniqueCounter = 0
reads = 0
lasttime = time.time()

# ---------- Unique Dictionary ------------- #
uniqueDict = {}

#------------------- functions ------------------------------
def barcodeDist(b_1, b_2):
    'counts mismatches between two equal-length strings'
    if len(b_1) == len(b_2) and len(b_1) > 0:
        return sum(map(lambda x: x[0] != x[1], zip(b_1, b_2) )) # python is bad-ass
    else:
        print "ERROR lengths of barcodes and index read do not match!"
        print "Target", b_1
        print "Index read:", b_2
        sys.exit()

#------------------- functions ------------------------------
def primerDist(read, primer):
    'counts mismatches between primer and sequence'
    read = read[0:(len(primer)-1)]
    return sum(map(lambda x: x[0] != x[1], zip(read, primer) )) # python is bad-ass

# ------- work horse function ---------------
try:
    while 1:
        read1 = R1.next()
        read2 = R2.next()
        read3 = R3.next()
        read4 = R4.next()
        reads += 1

        ### Barcode One Matching ###
        bc1 = None
        bc1Mismatch = False
        if read2.seq.tostring() in P7:
            bc1 = read2.seq.tostring()
        else:
            for key in P7:
                bcdist = barcodeDist(key, read2.seq.tostring())
                if bcdist <= barcodeMaxDiff:
                    bc1 = key
                    bc1Mismatch = True

        ### Barcode Two Matching ###
        bc2 = None
        bc2Mismatch = False
        if read3.seq.tostring() in P5:
            bc2 = read3.seq.tostring()
        else:
            for key in P5:
                bcdist = barcodeDist(key, read3.seq.tostring())
                if bcdist <= barcodeMaxDiff:
                    bc2 = key
                    bc2Mismatch = True

        ### Barcode Pair Matching ###
        combined_bc = None
        if "%s %s" % (bc1, bc2) in barcodes:
            combined_bc = barcodes["%s %s" % (bc1, bc2)]
            counters[combined_bc][0] += 1
            if bc1Mismatch:
                counters[combined_bc][1] += 1
            if bc2Mismatch:
                counters[combined_bc][2] += 1


        ### Primer Matching ###
        primer1 = None
        primer1Mismatch = 10
        for primer in primersP5.keys():
            pmismatches = primerDist(read1.seq.tostring(), primer)
            if pmismatches < primer1Mismatch:
                primer1 = primer
                primer1Mismatch = pmismatches
        if (primer1 != None):
            read1 = read1[len(primer1):]
            primer1 = primersP5[primer1]


        primer2 = None
        primer2Mismatch = 10
        for primer in primersP7.keys():
            pmismatches = primerDist(read4.seq.tostring(), primer)
            if pmismatches < primer2Mismatch:
                primer2 = primer
                primer2Mismatch = pmismatches
        if (primer2 != None):
            read4 = read4[len(primer2):]
            primer2 = primersP7[primer2]

        primer_id = None
        if primer1 == primer2 and primer1Mismatch <= primerMaxDiff and primer2Mismatch <= primerMaxDiff:
            primer_id = primer1

        ### Output Reads ###
        if combined_bc is not None and primer_id is not None:
            read1.id = read1.name =  "%s 1:N:0:%s-%s|%s|%s|%s|%s" % (read1.id.split()[0], read2.seq.tostring(), read3.seq.tostring(), combined_bc, primer_id, primer1Mismatch, primer2Mismatch)
            read1.description = ""
            SeqIO.write(read1, outf['identified'][0], "fastq")
            read4.id = read4.name =  "%s 2:N:0:%s-%s|%s|%s|%s|%s" % (read4.id.split()[0], read2.seq.tostring(), read3.seq.tostring(), combined_bc, primer_id, primer1Mismatch, primer2Mismatch)
            read4.description = ""
            SeqIO.write(read4, outf['identified'][1], "fastq")
            txt = '\t'.join([read1.id.split()[0], read2.seq.tostring(), read3.seq.tostring(), combined_bc, primer_id, str(primer1Mismatch), str(primer2Mismatch)]) + '\n'
            barcodesFile.write(txt)
            goodReadsCounter += 1
            comb_read = read1 + read4
            comb_read = "".join([comb_read.seq.tostring(), primer_id])
            if (comb_read in uniqueDict):
                uniqueDict[comb_read]['bc_counts'][combined_bc] += 1
            else:
                uniqueDict[comb_read] = {}
                uniqueDict[comb_read]['read1'] = read1.seq.tostring()
                uniqueDict[comb_read]['read2'] = read4.seq.tostring()
                uniqueDict[comb_read]['primer_id'] = primer_id
                uniqueDict[comb_read]['bc_counts'] = Counter()
                uniqueDict[comb_read]['bc_counts'][combined_bc] += 1
                uniqueCounter += 1

        else:
            read1.id = read1.name =  "%s 1:N:0:%s-%s|%s|%s|%s|%s|%s" % (read1.id.split()[0], read2.seq.tostring(), read3.seq.tostring(), combined_bc, primer1, primer1Mismatch, primer2, primer2Mismatch)
            read1.description = ""
            SeqIO.write(read1, outf['unidentified'][0], "fastq")
            read4.id = read4.name =  "%s 2:N:0:%s-%s|%s|%s|%s|%s|%s" % (read4.id.split()[0], read2.seq.tostring(), read3.seq.tostring(), combined_bc, primer1, primer1Mismatch, primer2, primer2Mismatch)
            read4.description = ""
            SeqIO.write(read4, outf['unidentified'][1], "fastq")
            otherCounter += 1


        ### Report every 100K ###
        if reads % 10000 == 0:
            print"-------------------"
            print "%s reads processed.Reads/second %s" % (reads, reads/(time.time() - lasttime))
            countersums = [0, 0, 0]
            for k in counters:
                countersums[0] += counters[k][0]
                countersums[1] += counters[k][1]
                countersums[2] += counters[k][2]
            print "identified:%s,\tmismatch bc1:%s,\tmismatch bc2:%s" % (countersums[0], countersums[1], countersums[2])
            print "good reads: %s" % goodReadsCounter
            print "other: %s" % otherCounter
            print "uniques: %s" % uniqueCounter
            print"-------------------"
except StopIteration:
    pass
finally:
    print "Finished processing "

    print"-------------------"
    print "%s reads processed in %s hours" % (reads,(time.time()-lasttime)/(60*60))
#    for k in IDS:
#        print "%s:\tidentified:%s,\tmismatch bc1:%s,\tmismatch bc2:%s" % (k, counters[k][0], counters[k][1], counters[k][2])
    print "good reads: %s" % goodReadsCounter
    print "other: %s" % otherCounter
    print "uniques: %s" % uniqueCounter
    print"-------------------"
    txt = '\t'.join(["read1" ,"read2" ,"primer_id"])
    for bc in IDS:
        txt = '\t'.join([txt, bc])
    txt = txt  + '\n'
    uniquesFile.write(txt)
    for k in uniqueDict:
        txt = '\t'.join([uniqueDict[k]['read1'], uniqueDict[k]['read2'], uniqueDict[k]['primer_id']])
        for bc in IDS:
            txt = '\t'.join([txt, str(uniqueDict[k]['bc_counts'][bc])])
        txt = txt  + '\n'
        uniquesFile.write(txt)
    for k in outf:
        outf[k][0].close()
        outf[k][1].close()
    barcodesFile.close()
    uniquesFile.close()
