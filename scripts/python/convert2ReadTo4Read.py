#!/usr/bin/env python

import os
import sys
import time
import argparse
import traceback
from dbcAmplicons import TwoReadIlluminaRun
from dbcAmplicons import IlluminaFourReadOutput

profile = False

# version 0.0.1
# initial release
# version 0.0.2
# modify for generic barcode lengths when processing reads/barcodes with Casava
version_num = "v0.0.2-11242015"


class convertApp:
    """
    Convert two read Illumina files (barcodes processed) back to a four read set to processed with dbcAmplicons
    """

    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, barcode1, barcode2, output_prefix, batchsize=100000, uncompressed=False, verbose=True, debug=False):
        """
        Start conversion of double barcoded Illumina sequencing run from two to four reads
        """
        self.verbose = verbose
        try:
            # setup output files
            self.run_out = IlluminaFourReadOutput(output_prefix, uncompressed)

            # establish and open the Illumin run
            self.run = TwoReadIlluminaRun(fastq_file1, fastq_file2)
            self.run.open()
            lasttime = time.time()
            while 1:
                # get next batch of reads
                reads = self.run.next(batchsize)
                if len(reads) == 0:
                    break
                # process individual reads
                for read in reads:
                    self.run_out.addRead(read.getFourReads(bc1_length=barcode1, bc2_length=barcode2))
                # Write out reads
                self.run_out.writeReads()
                if self.verbose:
                    sys.stderr.write("processed %s total reads, %s Reads/second\n" % (self.run.count(), round(self.run.count()/(time.time() - lasttime), 0)))
            if self.verbose:
                sys.stdout.write("%s reads processed in %s minutes\n" % (self.run.count(), round((time.time()-lasttime)/(60), 2)))
            # write out project table
            self.clean()
            return 0
        except (KeyboardInterrupt, SystemExit):
            self.clean()
            sys.stderr.write("%s unexpectedly terminated\n" % (__name__))
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
            self.run.close()
            self.run_out.close()
        except:
            pass


class convertCMD:
    """
    validate convertApp parser parameters and run the converion
    """
    def __init__(self):
        pass

    def execute(self, args):
        # ----------------------- options output prefix -----------------------
        output_prefix = args.output_base
        if args.output_base is None:
            output_prefix = "DBCreads"
        uncompressed = args.uncompressed
        # ----------------------- other options ------------
        debug = args.debug
        verbose = not args.verbose
        batchsize = args.batchsize
        barcode1 = args.barcode1
        barcode2 = args.barcode2

        app = convertApp()

        if profile:
            import cProfile
            cProfile.runctx('app.start(args.fastq_file1, arg.fastq_file2, barcode1, barcode2, output_prefix, batchsize, uncompressed, verbose, debug)', globals(), locals())
            return 255
        else:
            return app.start(args.fastq_file1, args.fastq_file2, barcode1, barcode2, output_prefix, batchsize, uncompressed, verbose, debug)


def parseArgs():
    """
    generate main parser
    """
    parser = argparse.ArgumentParser(
        description='convert2ReadTo4Read, a python script for back converting 2 Illumina sequence reads that have been processed for their barcodes back to a four read set to be reprocessed using dbcAmplicons',
        epilog='For questions or comments, please contact Matt Settles <msettles@uidaho.edu>', add_help=True)
    parser.add_argument('--version', action='version', version="%(prog)s Version " + version_num)
    parser.add_argument('-p', '--barcode1', help='Barcode 1 length [default: %(default)s], will become Read 2',
                        type=int, dest='barcode1', default=8)
    parser.add_argument('-q', '--barcode2', help='Barcode 2 length [default: %(default)s], will become Read 3',
                        type=int, dest='barcode2', default=8)
    parser.add_argument('-b', '--batchsize', help='batch size to process reads in [default: %(default)s]',
                        type=int, dest='batchsize', default=100000)
    parser.add_argument('-O', '--output_path', help='path for output files [default: %(default)s]',
                        action='store', type=str, dest='output_base', metavar='PREFIX', default=None)
    parser.add_argument('-u', '--uncompressed', help='leave output files uncompressed [default: %(default)s]',
                        action='store_true', dest='uncompressed', default=False)
    parser.add_argument('-v', '--silent', help='verbose output [default: %(default)s]',
                        action='store_true', dest='verbose', default=False)
    parser.add_argument('-1', metavar="read1", dest='fastq_file1', help='read1 of an amplicon fastq two file set',
                        action='store', type=str, required=True, nargs='+')
    parser.add_argument('-2', metavar="read2", dest='fastq_file2', help='read2 of an amplicon fastq two file set',
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

    convert = convertCMD()
    convert.execute(args)

if __name__ == '__main__':
    main()
