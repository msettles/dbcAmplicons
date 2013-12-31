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
import sys
import time
from dbcAmplicons import barcodeTable
from dbcAmplicons import primerTable
from dbcAmplicons import IlluminaRun
from dbcAmplicons import IlluminaOutput

class App:
    verbose = False
    evalPrimer = False
    def start(self, input_prefix, output_prefix, barcodesFile, primerFile, barcodeMaxDiff=1, primerMaxDiff=4, primerEndMatch=4, batchsize=10000, uncompressed=False,verbose=True):
        """
            Process double barcoded Illumina Sequencing Run
        """
        self.verbose = verbose
        self.evalPrimer = primerFile != None
        try:
            lasttime = time.time()
            ## read in primer sequences
            bcTable = barcodeTable(barcodesFile)
            if self.verbose:
                print "barcode table length: %s" % len(bcTable.barcodes)
            ## read in primer sequences if present
            if self.evalPrimer:
                prTable = primerTable(primerFile)
                if verbose:
                    print "primer table length P5 Primer Sequences:%s, P7 Primer Sequences:%s" % (len(prTable.P5sequences),len(prTable.P7sequences))
            ## setup output files
            self.run_out = IlluminaOutput(output_prefix,uncompressed)
            ## establish and open the Illumin run
            self.run = IlluminaRun(input_prefix)
            self.run.open()
            while 1:
                ## get next batch of reads
                reads = self.run.next(batchsize)
                if len(reads) == 0:
                    break
                ## process individual reads
                goodReads = [[],[]]
                badReads = [[],[]]
                for read in reads:
                    bcTable = read.getBarcode(bcTable,barcodeMaxDiff) ## barcode
                    if self.evalPrimer: ## primer
                        read.getPrimer(prTable,primerMaxDiff,primerEndMatch)
                    readOut = read.writeRead()
                    if read.goodRead == True:
                        self.run_out.identified_count += 1
                        goodReads[0].append(readOut[0])
                        goodReads[1].append(readOut[1])
                    else:
                        self.run_out.unidentified_count += 1
                        badReads[0].append(readOut[0])
                        badReads[1].append(readOut[1])

                self.run_out.writeGoodReads('\n'.join(goodReads[0]) + '\n','\n'.join(goodReads[1]) + '\n')
                self.run_out.writeBadReads('\n'.join(badReads[0]) + '\n','\n'.join(badReads[1]) + '\n')
                if self.verbose:
                    print "processed %s total reads, %s Reads/second, %s identified reads, %s unidentified reads" % (self.run.count, round(self.run.count/(time.time() - lasttime),0), self.run_out.identified_count,self.run_out.unidentified_count)
            if self.verbose:
                print "%s reads processed in %s minutes" % (reads,round((time.time()-lasttime)/(60),2))
                print "final barcode table length: %s" % len(bcTable.barcodes)
            self.clean()
            return 0	
        except Exception:
            self.clean()
            print("A fatal error was encountered.")
            print sys.exc_info()[0]
            return 1
        except (KeyboardInterrupt, SystemExit):
            self.clean()
            print("%s unexpectedly terminated" % (__name__))
            return 1

    def clean(self):
        if self.verbose:
            print("Cleaning up.")
        try:
            self.run.close()
            self.run_out.close()
        except:
            pass
