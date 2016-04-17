
import os
import sys
import traceback
import time
import signal

from subprocess import Popen
from subprocess import PIPE

from dbcAmplicons import IlluminaTwoReadOutput
from dbcAmplicons import IlluminaOneReadOutput


def sp_bowtie2_index(ref, overwrite=False):
    if os.path.isfile(ref):
        if os.path.isfile(ref + '.rev.2.bt2') and not overwrite:
            print 'Found existing bowtie2 index for %s' % ref
            return 0
        else:
            FNULL = open(os.devnull, 'w')
            call = 'bowtie2-build'
            call = call + ' ' + ref + ' ' + ref
            p = Popen(['bowtie2-build', ref, ref],
                      stdout=FNULL,
                      stderr=FNULL,
                      bufsize=-1,
                      preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
            p.communicate()
            if p.returncode:
                print 'Something in bowtie2-build went wrong'
                raise
            # system call, check for return
            print 'Successfully indexed %s' % ref
            return 0
    else:
        print "%s Reference file not found" % ref
        return 1
    print 'Something in bowtie2-build went wrong'
    raise


# Template
# bowtie2 -x caplanStuff -U <(zcat ../../CaplanShit/00-RawData/Sample_GCCAAT/GCCAAT_R1.fastq.gz| sed 's, ,_,g')
def sp_bowtie2_screen(pe1, pe2, se, ref, overwrite=False, sensitivity=0, procs=1):
    # build the call,
    # each file must first go through awk to replace spaces with a parsable character
    if sp_bowtie2_index(ref, overwrite) != 0:
        sys.exit(1)

    sensitivity_switch = ['--very-fast-local', '--fast-local', '--sensitive-local', '--very-sensitive-local']
    call = 'bowtie2 -I 0 -X 1500 ' + sensitivity_switch[sensitivity] + ' -p ' + str(procs) + ' -x ' + ref
    if ((pe1 is not None) and (pe2 is not None) and (len(pe1) == len(pe2))):
        pe1_gz = "gunzip -c"
        pe2_gz = "gunzip -c"
        pe1_gz_true = False
        pe2_gz_true = False
        pe1_ngz = "cat"
        pe2_ngz = "cat"
        pe1_ngz_true = False
        pe2_ngz_true = False
        for pe_read in pe1:
            if pe_read.split(".")[-1] == "gz":
                pe1_gz = pe1_gz + " " + pe_read
                pe1_gz_true = True
            else:
                pe1_ngz = pe1_ngz + " " + pe_read
                pe1_ngz_true = True
        for pe_read in pe2:
            if pe_read.split(".")[-1] == "gz":
                pe2_gz = pe2_gz + " " + pe_read
                pe2_gz_true = True
            else:
                pe2_ngz = pe2_ngz + " " + pe_read
                pe2_ngz_true = True
        if pe1_gz_true is True:
            call = call + " -1 <(" + pe1_gz + "| sed 's, ,_:_,g')"
        if pe2_gz_true is True:
            call = call + " -2 <(" + pe2_gz + "| sed 's, ,_:_,g')"
        if pe1_ngz_true is True:
            call = call + " -1 <(" + pe1_ngz + "| sed 's, ,_:_,g')"
        if pe2_ngz_true is True:
            call = call + " -2 <(" + pe2_ngz + "| sed 's, ,_:_,g')"
    if (se is not None):
        se_gz = "gunzip -c"
        se_gz_true = False
        se_ngz = "cat"
        se_ngz_true = False
        for se_read in se:
            if se_read.split(".")[-1] == "gz":
                se_gz = se_gz + " " + se_read
                se_gz_true = True
            else:
                se_ngz = se_ngz + " " + se_read
                se_ngz_true = True
        if se_gz_true is True:
            call = call + " -U <(" + se_gz + "| sed 's, ,_:_,g')"
        if se_ngz_true is True:
            call = call + " -U <(" + se_gz + "| sed 's, ,_:_,g')"
    print call
    FNULL = open(os.devnull, 'w')
    p = Popen(call,
              stdout=PIPE,
              stderr=None,
              bufsize=-1,
              shell=True,
              executable='/bin/bash',
              preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    if p.returncode:
        raise
    return p.stdout


def reverseComplement(s):
    """
    given a seqeucne with 'A', 'C', 'T', and 'G' return the reverse complement
    """
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    letters = list(s)
    try:
        letters = [basecomplement[base] for base in letters]
    except:
        raise
    return ''.join(letters[::-1])


def reverse(s):
    """
    given a sequence return the reverse
    """
    letters = list(s)
    return ''.join(letters[::-1])


class screeningApp:

    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, fastq_file3, reference, overwrite, sensitivity, output_prefix, strict, procs, uncompressed=False, verbose=True, debug=False):
        """
            screen reads against a reference fasta file
        """
        self.verbose = verbose
        try:
            mapped_pairs_count = 0
            unmapped_pairs_count = 0
            mapped_singles_count = 0
            unmapped_singles_count = 0
            secondary_alignment = 0

            # Set up output
            self.run_out = {}
            self.run_out["mapped_pairs"] = IlluminaTwoReadOutput(output_prefix + '.mapped', uncompressed)
            self.run_out["unmapped_pairs"] = IlluminaTwoReadOutput(output_prefix + '.unmapped', uncompressed)
            self.run_out["mapped_singles"] = IlluminaOneReadOutput(output_prefix + '.mapped', uncompressed)
            self.run_out["unmapped_singles"] = IlluminaOneReadOutput(output_prefix + '.unmapped', uncompressed)

            # 0x1 template having multiple segments in sequencing
            # 0x2 each segment properly aligned according to the aligner
            # 0x4 segment unmapped
            # 0x8 next segment in the template unmapped
            # 0x10 SEQ being reverse complemented
            # 0x20 SEQ of the next segment in the template being reversed
            # 0x40 the first segment in the template
            # 0x80 the last segment in the template
            # 0x100 secondary alignment
            # 0x200 not passing quality controls
            # 0x400 PCR or optical duplicate
            PE1 = {}
            PE2 = {}

            i = 0
            for line in sp_bowtie2_screen(fastq_file1, fastq_file2, fastq_file3, reference, overwrite, sensitivity, procs):
                if i % 100000 == 0 and i > 0:
                    for key in self.run_out:
                        self.run_out[key].writeReads()
                    if self.verbose:
                        print "Processed: %s, PE in ref: %s, SE in ref: %s" % (i, mapped_pairs_count, mapped_singles_count)
                if line[0] != "@":  # header line
                    i += 1

                    line2 = line.strip().split()
                    flag = int(line2[1])
                    # Secondary alignment
                    if (flag & 0x100):
                        secondary_alignment += 1
                        continue

                    # Handle SE:
                    # mapped SE reads have 0x1 set to 0, and 0x4 (third bit) set to 1
                    if not (flag & 0x1):  # SE READ
                        if not (flag & 0x4):  # MAPPED
                            ID = line2[0]
                            if (flag & 0x10):  # reverse complement
                                line2[9] = reverseComplement(line2[9])
                                line2[10] = reverse(line2[10])
                            self.run_out["mapped_singles"].addRead(['\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])])
                            mapped_singles_count += 1
                        else:  # UNMAPPED
                            ID = line2[0]
                            if (flag & 0x10):  # reverse complement
                                line2[9] = reverseComplement(line2[9])
                                line2[10] = reverse(line2[10])
                            self.run_out["unmapped_singles"].addRead(['\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])])
                            unmapped_singles_count += 1
                        continue
                    # Handle PE:
                    # logic:  0x1 = multiple segments in sequencing,   0x4 = segment unmapped,  0x8 = next segment unmapped
                    if (flag & 0x1):  # PE READ
                        if ((strict and not (flag & 0x4) and not (flag & 0x8)) or  # both pairs mapped
                                (not strict and (not (flag & 0x4) or not (flag & 0x8)))):  # at least one of the pair mapped
                            if (flag & 0x40):  # is this PE1 (first segment in template)
                                # PE1 read, check that PE2 is in dict and write out
                                ID = line2[0].split('_:_')[0]
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r1 = '\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])  # sequence + qual
                                if ID in PE2:
                                    self.run_out["mapped_pairs"].addRead([r1, PE2[ID]])
                                    del PE2[ID]
                                    mapped_pairs_count += 1
                                else:
                                    PE1[ID] = r1
                            elif (flag & 0x80):  # is this PE2 (last segment in template)
                                # PE2 read, check that PE1 is in dict and write out
                                ID = line2[0].split('_:_')[0]
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r2 = '\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])
                                if ID in PE1:
                                    self.run_out["mapped_pairs"].addRead([PE1[ID], r2])
                                    del PE1[ID]
                                    mapped_pairs_count += 1
                                else:
                                    PE2[ID] = r2
                        else:  # an 'unmapped' pair
                            if (flag & 0x40):  # is this PE1 (first segment in template)
                                # PE1 read, check that PE2 is in dict and write out
                                ID = line2[0].split('_:_')[0]
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r1 = '\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])  # sequence + qual
                                if ID in PE2:
                                    self.run_out["unmapped_pairs"].addRead([r1, PE2[ID]])
                                    del PE2[ID]
                                    unmapped_pairs_count += 1
                                else:
                                    PE1[ID] = r1
                            elif (flag & 0x80):  # is this PE2 (last segment in template)
                                # PE2 read, check that PE1 is in dict and write out
                                ID = line2[0].split('_:_')[0]
                                if (flag & 0x10):  # reverse complement
                                    line2[9] = reverseComplement(line2[9])
                                    line2[10] = reverse(line2[10])
                                r2 = '\n'.join(['@' + line2[0].replace('_:_', ' '), line2[9], '+', line2[10]])
                                if ID in PE1:
                                    self.run_out["unmapped_pairs"].addRead([PE1[ID], r1])
                                    del PE1[ID]
                                    unmapped_pairs_count += 1
                                else:
                                    PE2[ID] = r2
            # Write out reads
            for key in self.run_out:
                self.run_out[key].writeReads()

            print "Records processed: %s, PE in ref: %s, SE in ref: %s" % (i, mapped_pairs_count, mapped_singles_count)

            self.clean()
            return 0

        except (KeyboardInterrupt, SystemExit):
            self.clean()
            sys.stderr.write("%s unexpectedly terminated\n" % (__name__))
            return 1
        except:
            self.clean()
            if not debug:
                sys.stderr.write("A fatal error was encountered. trying turning on debug\n")
            if debug:
                sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
            return 1

    def clean(self):
        if self.verbose:
            sys.stderr.write("Cleaning up.\n")
        try:
            pass
        except:
            pass
