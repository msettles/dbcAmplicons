#!/usr/bin/env python

# Copyright 2013, Institute for Bioninformatics and Evolutionary Studies
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import time
import gzip
from dbcAmplicons import barcodeTable
from dbcAmplicons import IlluminaRun

class App:
    def start(self, input_prefix, output_prefix, barcodesFile="barcodeLookupTable.txt",barcodeMaxDiff=1, uncompressed=False,verbose=True):
        """
            Process double barcoded Illumina Sequencing Run
        """
        try:
            lasttime = time.time()
            bcTable = barcodeTable(barcodesFile)
            if verbose:
                print "barcode table length: %s" % len(bcTable.barcodes)
            run = IlluminaRun(input_prefix)
            run.open()
            if uncompressed is True:
                outf = {'identified':[open(output_prefix + '_R1.fastq', 'w'), open(output_prefix + '_R2.fastq', 'w')], 'unidentified':[open(output_prefix + '_Unidentified_R1.fastq', 'w'), open(output_prefix + '_Unidentified_R2.fastq', 'w')]}
            else:
                outf = {'identified':[gzip.open(output_prefix + '_R1.fastq.gz', 'wb'), gzip.open(output_prefix + '_R2.fastq.gz', 'wb')], 'unidentified':[gzip.open(output_prefix + '_Unidentified_R1.fastq.gz', 'wb'), gzip.open(output_prefix + '_Unidentified_R2.fastq.gz', 'wb')]}
            while 1:
                goodReads = [[],[]]
                badReads = [[],[]]
                reads = run.next(10000)
                if len(reads) == 0:
                    break
                for read in reads:
                    bcTable = read.getBarcode(bcTable,barcodeMaxDiff)
                    readOut = read.writeRead()
                    if read.bc_ID[0] is not None :
                        goodReads[0].append(readOut[0])
                        goodReads[1].append(readOut[1])
                    else:
                        badReads[0].append(readOut[0])
                        badReads[1].append(readOut[1])
                outf['unidentified'][0].write('\n'.join(badReads[0]) + '\n')
                outf['unidentified'][1].write('\n'.join(badReads[1]) + '\n')
                outf['identified'][0].write('\n'.join(goodReads[0]) + '\n')
                outf['identified'][1].write('\n'.join(goodReads[1]) + '\n')
            outf['unidentified'][0].close()
            outf['unidentified'][1].close()
            outf['identified'][0].close()
            outf['identified'][1].close()
            if verbose:
                print "%s reads processed: Reads/second %s" % (run.count, run.count/(time.time() - lasttime))
                print "barcode table length: %s" % len(bcTable.barcodes)
                print("Cleaning up.")
            run.close()
            self.clean()
            return 0	
        except Exception:
            print("A fatal error was encountered.")
            return 1
        except (KeyboardInterrupt, SystemExit):
            self.clean()
            print("%s unexpectedly terminated" % (__name__))
            return 1

    def clean(self):
        pass


