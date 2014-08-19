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
import os, sys, traceback
import time
from dbcAmplicons import sampleTable
from dbcAmplicons import TwoReadIlluminaRun
from dbcAmplicons import IlluminaTwoReadOutput

class splitreadsApp:
    verbose = False
    evalPrimer = False
    def start(self, fastq_file1, fastq_file2, output_prefix, samplesFile, batchsize = 10000, uncompressed = False, output_unidentified = False, verbose = True, debug = False):
        """
            split a double barcoded Illumina Sequencing Run by project
        """
        self.verbose = verbose
        try:
            ## read in primer sequences
            sTable = sampleTable(samplesFile)
            if self.verbose:
                sys.stdout.write("sample table length: %s, and %s projects.\n" % (sTable.getSampleNumber(),len(sTable.getProjectList())))
            ## read in primer sequences if present
            ## setup output files
            identified_count = 0
            unidentified_count = 0
            self.run_out = {}
            for project in sTable.getProjectList():
                self.run_out[project] = IlluminaTwoReadOutput(os.path.join(output_prefix,project),uncompressed)
            if output_unidentified:
                self.run_out["Unidentified"] = IlluminaTwoReadOutput(os.path.join(output_prefix,'UnidentifiedProject'),uncompressed)
            ## establish and open the Illumin run
            self.run = TwoReadIlluminaRun(fastq_file1, fastq_file2)
            self.run.open()
            lasttime = time.time()
            while 1:
                ## get next batch of reads
                reads = self.run.next(batchsize)
                if len(reads) == 0:
                    break
                ## process individual reads
                for read in reads:
                    read.assignRead(sTable) ## barcode
                    if read.goodRead == True:
                        self.run_out[read.project].addRead(read.getFastq())
                        identified_count +=1
                    else:
                        unidentified_count += 1
                        if output_unidentified:
                            self.run_out["Unidentified"].addRead(read.getFastq())
                ### Write out reads
                for key in self.run_out:
                    self.run_out[key].writeReads()
                if self.verbose:
                    sys.stderr.write("processed %s total reads, %s Reads/second, %s identified reads, %s unidentified reads (%s%%)\n" % (self.run.count(), round(self.run.count()/(time.time() - lasttime),0), identified_count,unidentified_count,round((float(identified_count)/float(self.run.count()))*100)))
            if self.verbose:
                sys.stdout.write("%s reads processed in %s minutes, %s (%s%%) identified\n\n" % (self.run.count(),round((time.time()-lasttime)/(60),2),identified_count,round((float(identified_count)/float(self.run.count()))*100,1)))
            for key in self.run_out:
                sys.stdout.write("%s (%s%%)\treads found for project\t%s\n" % (self.run_out[key].count(), round((float(self.run_out[key].count())/float(identified_count))*100,1), key))
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
            self.run.close()
            for key in self.run_out:
                self.run_out[key].close()
        except:
            pass

