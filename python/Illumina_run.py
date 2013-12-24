#!/usr/bin/env python
import sys
import os
import glob
import gzip
import time

# https://py-editdist.googlecode.com/files/py-editdist-0.3.tar.gz
try: 
	import editdist
	editdist_loaded = True 
except ImportError: 
	editdist_loaded = False 

# ---------------- barcodes class and make a dictionary for lookup ----------------
class barcodeTable:
	""" class to hold barcode information """
	def __init__(self, P5, P7, IDS, barcodes):
		self.P5 = P5
		self.P7 = P7
		self.IDS = IDS
		self.barcodes = barcodes

def parseBarcodes(bfile):
	barcodes = {}
	P5 = []
	P7 = []
	IDS = []
	bcfile = open(bfile,'rb')
	f = bcfile.readlines()[1:] # skip the first header line
	for row in f:
		row = row.rstrip()
		ID, P5Name, P5BC, P7Name, P7BC = row.split('\t')
		if P5BC not in P5:
			P5.extend([P5BC])
		if P7BC not in P7:
			P7.extend([P7BC])
		IDS.extend([ID])
		barcodes["%s%s" % (P7BC, P5BC)] = [ID,0,0]
	bcfile.close()
	return barcodeTable(P5, P7, IDS, barcodes)

#------------------- functions ------------------------------
def barcodeDist(b_1, b_2):
	'counts mismatches between two equal-length strings'
	if len(b_1) == len(b_2) and len(b_1) > 0:
		return sum(map(lambda x: x[0] != x[1], zip(b_1, b_2) ))
	else:
		print "ERROR lengths of barcodes and index read do not match!"
		print "Target", b_1
		print "Index read:", b_2
		sys.exit()

# ---------------- Class for 4 read sequence data with double barcode set ----------------
class SequenceReadSet:
	""" class to hold a read set """
	bc_ID =None
	def __init__(self,name,read_1,qual_1,read_2,qual_2,bc_1,bc_2):
		self.name = name
		self.read_1 = read_1
		self.qual_1 = qual_2
		self.read_2 = read_2
		self.qual_2 = qual_2
		self.bc_1 = bc_1
		self.bc_2 = bc_2
	def getBarcode(self, bcTable, max_diff):
		bc_pair = "%s%s" % (self.bc_1, self.bc_2)
		if bc_pair in bcTable.barcodes:
			self.bc_ID = bcTable.barcodes[bc_pair]
			return bcTable
		else:
		### Barcode One Matching ###
			bc1 = None
			bc1Mismatch = 0
			if self.bc_1 in bcTable.P7:
				bc1 = self.bc_1
			else:
				bc1Mismatch = len(self.bc_1)
				for key in bcTable.P7:
					if editdist_loaded :
						bcdist = editdist.distance(key, self.bc_1)
					else:
						bcdist = barcodeDist(key, self.bc_1)
					if bcdist < bc1Mismatch:
						bc1 = key
						bc1Mismatch = bcdist
			### Barcode Two Matching ###
			bc2 = None
			bc2Mismatch = 0
			if self.bc_2 in bcTable.P5:
				bc2 = self.bc_2
			else:
				bc2Mismatch = len(self.bc_2)
				for key in bcTable.P5:
					if editdist_loaded :
						bcdist = editdist.distance(key, self.bc_2)
					else:
						bcdist = barcodeDist(key, self.bc_2)
					if bcdist <= bc2Mismatch:
						bc2 = key
						bc2Mismatch = bcdist
			### Barcode Pair Matching ###
			combined_bc = [None,len(self.bc_1),len(self.bc_1)]
			if "%s%s" % (bc1, bc2) in bcTable.barcodes and bc1Mismatch <= max_diff and bc2Mismatch <= max_diff:
				combined_bc = [bcTable.barcodes["%s%s" % (bc1, bc2)][0],bc1Mismatch,bc2Mismatch]
			bcTable.barcodes[bc_pair] = combined_bc
			self.bc_ID = combined_bc
			return bcTable

class IlluminaRun:
	isOpen = False
	count = 0
	""" A seuqencing run """
	def __init__(self,file_prefix):
		self.prefix = file_prefix
		self.fread1 = glob.glob("%s*R1_[0-9][0-9][0-9].fastq.gz" % file_prefix)
		if len(self.fread1) == 0 or not all(os.path.exists(f) for f in self.fread1):
			raise Exception("R1 file not found")
		self.fread2 = glob.glob("%s*R2_[0-9][0-9][0-9].fastq.gz" % file_prefix)
		if len(self.fread2) == 0 or not all(os.path.exists(f) for f in self.fread2):
			raise Exception("R2 file not found")
		self.fread3 = glob.glob("%s*R3_[0-9][0-9][0-9].fastq.gz" % file_prefix)
		if len(self.fread3) == 0 or not all(os.path.exists(f) for f in self.fread3):
			raise Exception("R3 file not found")
		self.fread4 = glob.glob("%s*R4_[0-9][0-9][0-9].fastq.gz" % file_prefix)
		if len(self.fread4) == 0 or not all(os.path.exists(f) for f in self.fread4):
			raise Exception("R4 file not found")
		if any(len(f) != len(self.fread1) for f in [self.fread2,self.fread3,self.fread4]):
			raise Exception("Inconsistant number of files for each read")
		# record the number of files per read
		self.numberoffiles = len(self.fread1)
		if self.numberoffiles != 1 :
			raise Exception("Number of files per read exceeds 1")
	def open(self):
		if self.fread1[0].split(".")[-1] == "gz":
			self.R1 = gzip.open(self.fread1[0], 'rb')
		else:
			self.R1 = open(self.fread1[0], 'r')
		if self.fread2[0].split(".")[-1] == "gz":
			self.R2 = gzip.open(self.fread2[0], 'rb')
		else:
			self.R2 = open(self.fread2[0], 'r')
		if self.fread3[0].split(".")[-1] == "gz":
			self.R3 = gzip.open(self.fread3[0], 'rb')
		else:
			self.R3 = open(self.fread3[0], 'r')
		if self.fread4[0].split(".")[-1] == "gz":
			self.R4 = gzip.open(self.fread4[0], 'rb')
		else:
			self.R4 = open(self.fread4[0], 'r')
		self.isOpen = True
		self.count=0
	def close(self):
		self.R1.close()
		self.R2.close()
		self.R3.close()
		self.R4.close()
		self.isOpen = False
	def next(self, count=1):
		if not self.isOpen:
			self.open()
		reads = []
		for i in range(0,count):
			try:
				#process read_1 (Read1)
				name = self.R1.next().split(" ")[0] # name
				read_1 = self.R1.next().rstrip()	# read
				self.R1.next()	# '+'
				qual_1 = self.R1.next().rstrip() # qual
				#process read _2 (BC1)
				self.R2.next() # name
				bc_1 = self.R2.next().rstrip() # read 
				self.R2.next() # '+'
				self.R2.next() # qual
				#process reed_3 (BC2)
				self.R3.next() # name
				bc_2 = self.R3.next().rstrip() # read
				self.R3.next() # '+'
				self.R3.next() # qual
				#process read_4 (Read2)
				self.R4.next() # name
				read_2 = self.R4.next().rstrip() # read
				self.R4.next() # '+'
				qual_2 = self.R4.next().rstrip() # qual
				#add it to the stack
			except StopIteration:
				break
			reads.append(SequenceReadSet(name=name,read_1=read_1,qual_1=qual_1,read_2=read_2,qual_2=qual_2,bc_1=bc_1,bc_2=bc_2))
			self.count += 1
		return(reads)

barcodeFile = "barcodeLookupTable.txt"
bcTable = parseBarcodes(barcodeFile)
print "barcode table length: %s" % len(bcTable.barcodes)
seqsFile = "Amplicon_Raw_fastq/test40k"
run = IlluminaRun(seqsFile)
run.open()

lasttime = time.time()

while 1:
	reads = run.next(10000)
	if len(reads) == 0:
		break
	for read in reads:
		bcTable = read.getBarcode(bcTable,1)
		print "processed read: %s" % read.bc_ID[0]

print "%s reads processed: Reads/second %s" % (run.count, run.count/(time.time() - lasttime))
print "barcode table length: %s" % len(bcTable.barcodes)

run.close()

