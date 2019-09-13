#!/usr/bin/env python

import os
import sys
import argparse
import traceback
from dbcAmplicons import barcodeTable
from dbcAmplicons import sampleTable
from dbcAmplicons import TwoReadIlluminaRun
from dbcAmplicons import IlluminaTwoReadOutput
from dbcAmplicons import OneReadIlluminaRun
from dbcAmplicons import IlluminaOneReadOutput

profile = False

# version 0.0.1
# initial release

version_num = "v0.0.1"


class splitApp:
    """
    Split two read Illumina files (barcodes processed from dbcAmplicons) into separate projects
    """

    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, fastq_fileU, basespace, barcodes_file, barcodediff, I1_orientation, I2_orientation, barcode1, barcode2, samples_file, output_prefix, batchsize=100000, uncompressed=False, verbose=True, debug=False):
        """
        Split double barcoded Illumina sequencing run from two to four reads by sample identifier
        """
        self.verbose = verbose
        evalBarcode = barcodes_file is not None
        evalSample = samples_file is not None
        if evalBarcode:
            bcTable = barcodeTable(barcodes_file, I1_orientation, I2_orientation)
            if self.verbose:
                sys.stdout.write("barcode table length: %s\n" % bcTable.getLength())
        if evalSample:
            sTable = sampleTable(samples_file)
            if verbose:
                sys.stdout.write("sample table length: %s, and %s projects.\n" % (sTable.getSampleNumber(), len(sTable.getProjectList())))

        if fastq_fileU is not None and (fastq_file1 is not None and fastq_file2 is not None):
            sys.stderr.write("ERROR:[SplitBySample] cannot have both paired and single reads\n")
            return 1
        try:
            if fastq_file1 is not None and fastq_file2 is not None:
                self.runPairs = TwoReadIlluminaRun(fastq_file1, fastq_file2)
                self.runPairs.open()
            else:
                self.runPairs = None
            if fastq_fileU is not None:
                self.runSingle = OneReadIlluminaRun(fastq_fileU)
                self.runSingle.open()
            else:
                self.runSingle = None
            if self.runPairs is None and self.runSingle is None:
                sys.stderr.write("ERROR:[SplitBySample] input reads not specified, or incorrect pairs\n")
                raise Exception

            self.barcode_table = open(output_prefix + ".table.txt", 'w')

            self.run_out = {}
            sample_bc = {}

            if (self.runPairs is not None):
                while 1:
                    if self.verbose:
                        sys.stderr.write("Processing sequence files.\n")
                    # get next batch of reads
                    reads = self.runPairs.next(batchsize)
                    if len(reads) == 0:
                        break
                    # process individual reads, check to see if sample was already added to the library of self.run_out
                    for read in reads:
                        if basespace:
                            read_bc = read.getBarcodeBaseSpace(bc1_length=barcode1, bc2_length=barcode2)
                            if evalBarcode:
                                read.assignBarcode(bcTable, barcodediff)  # barcode
                            if evalSample:  # sample
                                read.assignRead(sTable, override_primer=True)  # barcode + primer

                            if read_bc in sample_bc:
                                sample_bc[read_bc] += 1
                            else:
                                sample_bc[read_bc] = 1
                        if read.goodRead is True:
                            sample = read.sample
                            if sample in self.run_out:
                                self.run_out[sample].addRead(read.getFastqSRA())
                            else:
                                self.run_out[sample] = IlluminaTwoReadOutput(os.path.join(output_prefix, sample), uncompressed)
                                self.run_out[sample].addRead(read.getFastqSRA())
                    # Write out reads for each key in dictionary
                    for key in self.run_out:
                        self.run_out[key].writeReads()
                if self.verbose:
                    sys.stderr.write("\nSplit out %s total samples in %s.\n" % (len(self.run_out), output_prefix))
                # write out the barcode
                max_bc = max(sample_bc, key=sample_bc.get)
                self.barcode_table.write("%s\t%s\t%s\n" % (sample_name, max_bc[0], max_bc[1]))
                self.barcode_table.flush()
                total_count += self.run.count()
                return 0
            if (self.runSingle is not None):
                while 1:
                    if self.verbose:
                        sys.stderr.write("Processing sequence files.\n")
                    # get next batch of reads
                    reads = self.runSingle.next(batchsize)
                    if len(reads) == 0:
                        break
                    # process individual reads, check to see if sample was already added to the library of self.run_out
                    for read in reads:
                        if basespace:
                            read_bc = read.getBarcodeBaseSpace(bc1_length=barcode1, bc2_length=barcode2)
                            if evalBarcode:
                                read.assignBarcode(bcTable, barcodediff)  # barcode
                            if evalSample:  # sample
                                read.assignRead(sTable, override_primer=True)  # barcode + primer

                            if read_bc in sample_bc:
                                sample_bc[read_bc] += 1
                            else:
                                sample_bc[read_bc] = 1

                        if read.goodRead is True:
                            sample = read.sample
                            if sample in self.run_out:
                                self.run_out[sample].addRead(read.getFastqSRA())
                            else:
                                self.run_out[sample] = IlluminaOneReadOutput(os.path.join(output_prefix, sample), uncompressed)
                                self.run_out[sample].addRead(read.getFastqSRA())
                    # Write out reads for each key in dictionary
                    for key in self.run_out:
                        self.run_out[key].writeReads()
                if self.verbose:
                    sys.stderr.write("\nSplit out %s total samples in %s.\n" % (len(self.run_out), output_prefix))
                # write out the barcode
                max_bc = max(sample_bc, key=sample_bc.get)
                self.barcode_table.write("%s\t%s\t%s\n" % (sample_name, max_bc[0], max_bc[1]))
                self.barcode_table.flush()
                total_count += self.run.count()
                return 0
        except (KeyboardInterrupt, SystemExit):
            self.clean()
            sys.stderr.write("%s unexpectedly terminated.\n" % (__name__))
            return 1
        except:
            self.clean()
            sys.stderr.write("A fatal error was encountered.\n")
            if debug:
                sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
            return 1

    def clean(self):
        if self.verbose:
            sys.stderr.write("Cleaning up.\n")
        try:
            self.runPairs.close()
            self.runSingle.close()
            self.run_out.close()
            self.barcode_table.close()
        except:
            pass


class splitCMD:
    """
    validate splitApp parser parameters and run the conversion
    """
    def __init__(self):
        pass

    def execute(self, args):
        # ----------------------- options input files -----------------------
        verbose = args.verbose
        if args.barcodes_file is None:
            bcFile = None
            if verbose:
                sys.stderr.write("No barcodes file provided\n")
        else:
            bcFile = args.barcodes_file
        if args.samples_file is None:
            sFile = None
            if verbose:
                sys.stderr.write("No sample file provided\n")
        else:
            sFile = args.samples_file
        # ----------------------- options output prefix -----------------------
        if args.output_base is None:
            output_prefix = "DBCsample"
        elif args.output_base is not None:
            output_prefix = args.output_base
        uncompressed = args.uncompressed
        # ----------------------- other options ------------
        debug = args.debug
        verbose = not args.verbose
        batchsize = args.batchsize

        app = splitApp()

        if profile:
            import cProfile
            cProfile.runctx('app.start(args.fastq_file1, arg.fastq_file2, arg.fastq_fileU,  args.basespace, bcFile, args.barcodediff, args.I1_orientation, args.I2_orientation, args.barcode1, args.barcode2, sFile, output_prefix, batchsize, uncompressed, verbose, debug)', globals(), locals())
            return 255
        else:
            return app.start(args.fastq_file1, args.fastq_file2, args.fastq_fileU, args.basespace, bcFile, args.barcodediff, args.I1_orientation, args.I2_orientation, args.barcode1, args.barcode2, sFile, output_prefix, batchsize, uncompressed, verbose, debug)


def parseArgs():
    """
    generate main parser
    """
    parser = argparse.ArgumentParser(
        description='SplitReadsBySample, a python script for splitting out Illumina sequence reads already processed by dbcAmplicons by sample',
        epilog='For questions or comments, please contact Alida Gerritsen <alida@uidaho.edu>', add_help=True)
    parser.add_argument('--version', action='version', version="%(prog)s Version " + version_num)
    parser.add_argument('-b', '--batchsize', help='batch size to process reads in [default: %(default)s]',
                        type=int, dest='batchsize', default=100000)
    parser.add_argument('-O', '--output_path', help='filename for output files [default: %(default)s]',
                        action='store', type=str, dest='output_base', metavar='PREFIX', default=None)
    parser.add_argument('-u', '--uncompressed', help='leave output files uncompressed [default: %(default)s]',
                        action='store_true', dest='uncompressed', default=False)
    parser.add_argument('-v', '--silent', help='verbose output [default: %(default)s]',
                        action='store_true', dest='verbose', default=False)
    parser.add_argument('-f', '--bcl2fastq', help='Parse bcl2fastq/basespace barcodes [default: %(default)s]',
                        action='store_true', dest='basespace', default=False)
    parser.add_argument('-B', '--barcodes_file', help='file with barcodes',
                                   action='store', type=str, dest='barcodes_file', metavar='FILENAME', required=False)
    parser.add_argument('-d', '--barcodediff', help='max hamming dist from barcode [default: %(default)s]',
                                   type=int, dest='barcodediff', default=1)
    parser.add_argument('--I1', help='processing of Index 1 read, leave as is, or reverse complement. [For four-read datasets, default: %(default)s. If --inline is specified, default: asis]',
                                   action='store', type=str, dest='I1_orientation', metavar='asis/rc', choices=['asis', 'rc'], default='rc')
    parser.add_argument('--I2', help='processing of Index 2 read, leave as is (asis), or reverse complement (rc) [default: %(default)s]',
                                   action='store', type=str, dest='I2_orientation', metavar='asis/rc', choices=['asis', 'rc'], default='asis')
    parser.add_argument('-p', '--barcode1', help='Barcode 1 length [default: %(default)s], will become Read 2',
                        type=int, dest='barcode1', default=8)
    parser.add_argument('-q', '--barcode2', help='Barcode 2 length [default: %(default)s], will become Read 3',
                        type=int, dest='barcode2', default=8)
    parser.add_argument('-S', '--sample_metadata', help='file with sample metadata',
                                   action='store', type=str, dest='samples_file', metavar='FILENAME', default=None)
    parser.add_argument('-1', metavar="read1", dest='fastq_file1', help='read1 of an amplicon fastq two file set',
                        action='store', type=str, required=False, nargs='+')
    parser.add_argument('-2', metavar="read2", dest='fastq_file2', help='read2 of an amplicon fastq two file set',
                        action='store', type=str, required=False, nargs='+')
    parser.add_argument('-U', metavar="singleRead", dest='fastq_fileU', help='fastq from a 1 file Illlumina run',
                        action='store', type=str, required=False, nargs='+')
    parser.add_argument('--debug', help='show traceback on error',
                        action='store_true', dest="debug", default=False)

    args = parser.parse_args()

    return args


def main():
    """
    main function
    """
    lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../')
    if lib_path not in sys.path:
        sys.path.insert(0, lib_path)
    args = parseArgs()
    split = splitCMD()
    split.execute(args)

if __name__ == '__main__':
    main()
