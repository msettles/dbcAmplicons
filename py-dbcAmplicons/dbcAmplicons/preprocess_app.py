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
import sys, os, traceback
import time
from dbcAmplicons import barcodeTable
from dbcAmplicons import primerTable
from dbcAmplicons import sampleTable
from dbcAmplicons import FourReadIlluminaRun
from dbcAmplicons import IlluminaTwoReadOutput

class preprocessApp:
    verbose = False
    evalPrimer = False
    def start(self, input_prefix, output_prefix, barcodesFile, primerFile, samplesFile, barcodeMaxDiff=1, primerMaxDiff=4, primerEndMatch=4, batchsize=10000, uncompressed=False, output_unidentified=False, verbose=True):
        """
            Process double barcoded Illumina Sequencing Run
        """
        self.verbose = verbose
        self.evalPrimer = primerFile != None
        self.evalSample = samplesFile != None
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
            if self.evalSample:
                sTable = sampleTable(samplesFile)
                if verbose:
                    print "sample table length: %s, and %s projects." % (sTable.sampleCount,len(sTable.projectList))
            ## setup output files
            barcode_counts = {}
            identified_count = 0
            unidentified_count = 0
            self.run_out = {}
            if self.evalSample:
                for project in sTable.projectList:
                    self.run_out[project] = IlluminaTwoReadOutput(os.path.join(output_prefix,project),uncompressed)
            else:
                self.run_out["Identified"] = IlluminaTwoReadOutput(output_prefix,uncompressed)
            if output_unidentified:
                self.run_out["Unidentified"] = IlluminaTwoReadOutput(output_prefix+"_Unidentified",uncompressed)
            ## establish and open the Illumin run
            self.run = FourReadIlluminaRun(input_prefix)
            self.run.open()
            while 1:
                ## get next batch of reads
                reads = self.run.next(batchsize)
                if len(reads) == 0:
                    break
                ## process individual reads
                for read in reads:
                    # bcTable = read.getBarcode(bcTable,barcodeMaxDiff) ## barcode
                    read.assignBarcode(bcTable,barcodeMaxDiff) ## barcode
                    if self.evalPrimer and read.goodRead: ## primer
                        read.assignPrimer(prTable,primerMaxDiff,primerEndMatch)
                    if self.evalSample and read.goodRead: ## sample
                        read.assignRead(sTable,True) ## barcode
                    if read.goodRead == True:
                        identified_count += 1
                        if self.evalSample:
                            self.run_out[read.getProject()].appendRead(read.getRead())
                        else:
                            self.run_out["Identified"].appendRead(read.getRead())
                        if read.getBarcode() in barcode_counts:
                            barcode_counts[read.getBarcode()]["Total"] += 1
                            if self.evalPrimer:
                                barcode_counts[read.getBarcode()][read.getPrimer()] += 1
                        else:
                            # setup blank primer count table
                            barcode_counts[read.getBarcode()] = {}
                            barcode_counts[read.getBarcode()]["Total"] = 1
                            if self.evalPrimer:
                                for pr in prTable.primers:
                                    barcode_counts[read.getBarcode()][pr] = 0
                                barcode_counts[read.getBarcode()][read.getPrimer()] += 1
                    else:
                        unidentified_count += 1
                        if output_unidentified:
                            self.run_out["Unidentified"].appendRead(read.getRead())
                ### Write out reads
                for key in self.run_out:
                    self.run_out[key].writeReads()
                if self.verbose:
                    print "processed %s total reads, %s Reads/second, %s identified reads, %s unidentified reads" % (self.run.count, round(self.run.count/(time.time() - lasttime),0), identified_count,unidentified_count)
            if self.verbose:
                print "%s reads processed in %s minutes" % (self.run.count,round((time.time()-lasttime)/(60),2))
            # Write out barcode and primer table
            if (identified_count > 0):
                bcFile = open(output_prefix + '_Identified_Barcodes.txt', 'w')
                bckeys = barcode_counts.keys()
                txt = '\t'.join(barcode_counts[bckeys[0]].keys())
                txt = "Barcode\t" + txt + '\n'
                bcFile.write(txt)
                for bc in bckeys:
                    txt = str(bc)
                    for pr in barcode_counts[bc]:
                        txt = '\t'.join([txt, str(barcode_counts[bc][pr])])
                    txt = txt  + '\n'
                    bcFile.write(txt)
            self.clean()
            return 0    
        except (KeyboardInterrupt, SystemExit):
            self.clean()
            print("%s unexpectedly terminated" % (__name__))
            return 1
        except:
            self.clean()
            print("A fatal error was encountered.")
            print "".join(traceback.format_exception(*sys.exc_info()))
            #log(str(e) + "".join(traceback.format_exception(*sys.exc_info())), out)
            return 1

    def clean(self):
        if self.verbose:
            print("Cleaning up.")
        try:
            self.run.close()
            for key in self.run_out:
                self.run_out[key].close()
        except:
            pass
