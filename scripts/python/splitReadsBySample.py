#!/usr/bin/env python

# Copyright 2014, Institute for Bioninformatics and Evolutionary Studies
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an 'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys
import time
import argparse
import traceback
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
    def start(self, fastq_file1, fastq_file2, fastq_fileU, output_prefix, batchsize=100000, uncompressed=False, verbose=True, debug=False):
        """
        Split double barcoded Illumina sequencing run from two to four reads by sample identifier
        """
        self.verbose = verbose
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
            self.run_out = {}
            if (self.runPairs is not None):
                while 1:
                    ## get next batch of reads
                    reads = self.runPairs.next(batchsize)
                    if len(reads) == 0:
                        break
                    ## process individual reads, check to see if sample was already added to the library of self.run_out
                    for read in reads:
                        sample = read.sample
                        if sample in self.run_out:
                            self.run_out[sample].addRead(read.getFastq())
                        else:
                            self.run_out[sample] = IlluminaTwoReadOutput(os.path.join(output_prefix, sample),uncompressed)
                            self.run_out[sample].addRead(read.getFastq())
                    ### Write out reads for each key in dictionary
                    for key in self.run_out:
                        self.run_out[key].writeReads()
                        if self.verbose:
                            sys.stderr.write("processing %s total samples\n" % len(self.run_out))
                return 0
                if self.verbose:
                    sys.stderr.write("processing %s total samples\n" % len(self.run_out))
            if (self.runSingle is not None):
                while 1:
                    ## get next batch of reads
                    reads = self.runSingle.next(batchsize)
                    if len(reads) == 0:
                        break
                    ## process individual reads, check to see if sample was already added to the library of self.run_out
                    for read in reads:
                        sample = read.sample
                        if sample in self.run_out:
                            self.run_out[sample].addRead(read.getFastq())
                        else:
                            self.run_out[sample] = IlluminaOneReadOutput(os.path.join(output_prefix, sample),uncompressed)
                            self.run_out[sample].addRead(read.getFastq())
                    ### Write out reads for each key in dictionary
                    for key in self.run_out:
                        self.run_out[key].writeReads()
                        if self.verbose:
                            sys.stderr.write("processing %s total samples\n" % len(self.run_out))
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
            self.runPairs.close()
            self.runSingle.close()
            self.run_out.close()
        except:
            pass

class splitCMD:
    """
    validate splitApp parser parameters and run the conversion
    """
    def __init__(self):
        pass
    def execute (self,args):
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
            cProfile.runctx('app.start(args.fastq_file1, arg.fastq_file2, arg.fastq_fileU, output_prefix, batchsize, uncompressed, verbose, debug)', globals(), locals())
            return 255
        else:
            return app.start(args.fastq_file1, args.fastq_file2, args.fastq_fileU, output_prefix, batchsize, uncompressed, verbose, debug)

#
#####################################################################################
##  Master parser arguments
def parseArgs():
    """
    generate main parser
    """
    parser = argparse.ArgumentParser( \
        description = 'SplitReadsBySample, a python script for splitting out Illumina sequence reads already processed by dbcAmplicons by sample', \
        epilog ='For questions or comments, please contact Alida Gerritsen <alida@uidaho.edu>', add_help=True)
    parser.add_argument('--version', action='version', version="%(prog)s Version " + version_num)
    parser.add_argument('-b', '--batchsize', help='batch size to process reads in [default: %(default)s]',
                        type=int, dest='batchsize', default=100000)    
    parser.add_argument('-O', '--output_path', help='filename for output files [default: %(default)s]',
                        action='store', type=str, dest='output_base', metavar='PREFIX', default=None)
    parser.add_argument('-u', '--uncompressed', help='leave output files uncompressed [default: %(default)s]',
                        action='store_true', dest='uncompressed', default=False)
    parser.add_argument('-v', '--silent', help='verbose output [default: %(default)s]',
                        action='store_true', dest='verbose', default=False)
    parser.add_argument('-1', metavar="read1", dest='fastq_file1', help='read1 of an amplicon fastq two file set',
                        action='store',type=str, required=False, nargs='+')
    parser.add_argument('-2', metavar="read2", dest='fastq_file2', help='read2 of an amplicon fastq two file set',
                        action='store',type=str, required=False, nargs='+')
    parser.add_argument('-U', metavar="singleRead", dest='fastq_fileU', help='fastq from a 1 file Illlumina run',
                        action='store',type=str, required=False, nargs='+')
    parser.add_argument('--debug', help='show traceback on error',
                        action='store_true', dest="debug", default = False)

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
 