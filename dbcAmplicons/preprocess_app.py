# preprocess_app.py
#
import sys
import os
import traceback
import time
from dbcAmplicons import barcodeTable
from dbcAmplicons import primerTable
from dbcAmplicons import sampleTable
from dbcAmplicons import FourReadIlluminaRun
from dbcAmplicons import IlluminaTwoReadOutput
from dbcAmplicons import validateApp
from dbcAmplicons import misc


class preprocessApp:
    """
    Preprocess raw Illumina four read amplicon data
    """

    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, fastq_file3, fastq_file4, output_prefix, barcodesFile, primerFile, samplesFile, barcodeMaxDiff=1, dedup_float=4, primerMaxDiff=4, primerEndMatch=4, batchsize=10000, uncompressed=False, output_unidentified=False, minQ=None, minL=0, verbose=True, debug=False, kprimer=False, test=False):
        """
        Start preprocessing double barcoded Illumina sequencing run, perform
        """
        self.verbose = verbose
        evalPrimer = primerFile is not None
        evalSample = samplesFile is not None
        try:
            v = validateApp()
            # read in primer sequences
            bcTable = barcodeTable(barcodesFile)
            if self.verbose:
                sys.stdout.write("barcode table length: %s\n" % bcTable.getLength())
            # read in primer sequences if present
            if evalPrimer:
                prTable = primerTable(primerFile)
                if verbose:
                    sys.stdout.write("primer table length P5 Primer Sequences:%s, P7 Primer Sequences:%s\n" % (len(prTable.getP5sequences()), len(prTable.getP7sequences())))
                if v.validatePrimer(prTable, debug) != 0:
                    sys.stderr.write("Failed validation\n")
                    self.clean()
                    return 1
            if evalSample:
                sTable = sampleTable(samplesFile)
                if verbose:
                    sys.stdout.write("sample table length: %s, and %s projects.\n" % (sTable.getSampleNumber(), len(sTable.getProjectList())))
                if v.validateSample(bcTable, prTable, sTable, debug) != 0:
                    sys.stderr.write("Failed validation\n")
                    self.clean()
                    return 1

            # output table
            try:
                if evalSample:
                    bctable_name = os.path.join(output_prefix, 'Identified_Barcodes.txt')
                else:
                    bctable_name = output_prefix + '_Identified_Barcodes.txt'
                misc.make_sure_path_exists(os.path.dirname(bctable_name))
                bcFile = open(bctable_name, 'w')
            except:
                sys.stderr.write("ERROR: Can't open file %s for writing\n" % bctable_name)
                raise
            barcode_counts = {}
            bcsuccesscount = 0
            prsuccesscount = 0
            sampsuccesscount = 0
            trimsuccesscount = 0
            identified_count = 0
            # setup output files
            self.run_out = {}
            if evalSample:
                for project in sTable.getProjectList():
                    self.run_out[project] = IlluminaTwoReadOutput(os.path.join(output_prefix, project), uncompressed)
            else:
                self.run_out["Identified"] = IlluminaTwoReadOutput(output_prefix, uncompressed)
            if output_unidentified:
                if evalSample:
                    self.run_out["Unidentified"] = IlluminaTwoReadOutput(os.path.join(output_prefix, 'UnidentifiedProject'), uncompressed)
                else:
                    self.run_out["Unidentified"] = IlluminaTwoReadOutput(output_prefix+"_Unidentified", uncompressed)
            # establish and open the Illumina run
            self.run = FourReadIlluminaRun(fastq_file1, fastq_file2, fastq_file3, fastq_file4)
            self.run.open()
            totaltime = time.time()
            while 1:
                lasttime = time.time()
                # get next batch of reads
                reads = self.run.next(batchsize)
                if len(reads) == 0:
                    break
                # process individual reads
                for read in reads:
                    bcsuccesscount += read.assignBarcode(bcTable, barcodeMaxDiff)  # barcode
                    if evalPrimer:  # primer
                        prsuccesscount += read.assignPrimer(prTable, dedup_float, primerMaxDiff, primerEndMatch)
                    if evalSample:  # sample
                        sampsuccesscount += read.assignRead(sTable)  # barcode + primer
                    if minQ is not None and read.goodRead:
                        trimsuccesscount += read.trimRead(minQ, minL)
                    if read.goodRead is True:
                        identified_count += 1
                        if evalSample:
                            self.run_out[read.getProject()].addRead(read.getFastq(kprimer))
                        else:
                            self.run_out["Identified"].addRead(read.getFastq(kprimer))
                    else:
                        if output_unidentified:
                            self.run_out["Unidentified"].addRead(read.getFastq(True))
                    ###############################################
                    # Record data for final barcode table
                    if read.getBarcode() == None and '-' in barcode_counts:
                        if evalPrimer and read.getPrimer() == None:
                            barcode_counts['-']['-'] += 1
                        elif evalPrimer:
                            barcode_counts['-'][read.getPrimer()] += 1
                        else:
                            barcode_counts['-']["Total"] += 1
                    elif read.getBarcode() in barcode_counts:
                        if evalPrimer and read.getPrimer() == None:
                            barcode_counts[read.getBarcode()]['-'] += 1
                        elif evalPrimer:
                            barcode_counts[read.getBarcode()][read.getPrimer()] += 1
                        else:
                            barcode_counts[read.getBarcode()]["Total"] += 1
                    else:
                        # setup blank primer count table for the new barcode
                        if read.getBarcode() == None:
                            barcode_counts['-'] = {}
                            if evalPrimer:
                                for pr in prTable.getPrimers():
                                    barcode_counts['-'][pr] = 0
                                barcode_counts['-']['-'] = 0
                                if read.getPrimer() == None:
                                    barcode_counts['-']['-'] = 1
                                else:
                                    barcode_counts['-'][read.getPrimer()] = 1
                            else:
                                barcode_counts['-']["Total"] = 1
                        else:
                            barcode_counts[read.getBarcode()] = {}
                            if evalPrimer:
                                for pr in prTable.getPrimers():
                                    barcode_counts[read.getBarcode()][pr] = 0
                                barcode_counts[read.getBarcode()]['-'] = 0
                                if read.getPrimer() == None:
                                    barcode_counts[read.getBarcode()]['-'] = 1
                                else:
                                    barcode_counts[read.getBarcode()][read.getPrimer()] = 1
                            else:
                                barcode_counts[read.getBarcode()]["Total"] = 1

                # Write out reads
                for key in self.run_out:
                    self.run_out[key].writeReads()
                if self.verbose:
                    sys.stderr.write("processed %s total reads, %s Reads/second, %s identified reads(%s%%), %s unidentified reads\n" % (self.run.count(), round(batchsize/(time.time() - lasttime), 0), identified_count, round((float(identified_count)/float(self.run.count()))*100, 1), self.run.count()-identified_count))
                if test:  # exit after the first batch to test the inputs
                    break
            if self.verbose:
                    sys.stdout.write("%s reads processed in %s minutes, %s (%s%%) identified\n\n" % (self.run.count(), round((time.time()-totaltime)/(60), 2), identified_count, round((float(identified_count)/float(self.run.count()))*100, 1)))
            # Write out barcode and primer table
            if (identified_count > 0):
                # write out header line
                if evalPrimer:
                    txt = 'Barcode\t' + '\t'.join(prTable.getPrimers()) + '\tNone' + '\n'
                else:
                    txt = 'Barcode\tTotal\n'
                bcFile.write(txt)
                bckeys = barcode_counts.keys()
                for bc in bcTable.getBarcodes():
                    if bc in bckeys and evalPrimer:
                        txt = str(bc)
                        for pr in prTable.getPrimers():
                            txt = '\t'.join([txt, str(barcode_counts[bc][pr])])
                        txt = "\t".join([txt, str(barcode_counts[bc]['-'])])
                    elif bc in bckeys:
                        txt = "\t".join([str(bc), str(barcode_counts[bc]["Total"])])
                    else:
                        continue
                    bcFile.write(txt + '\n')
                if '-' in bckeys:
                    if evalPrimer:
                        txt = 'None'
                        for pr in prTable.getPrimers():
                            txt = '\t'.join([txt, str(barcode_counts['-'][pr])])
                        txt = "\t".join([txt, str(barcode_counts['-']['-'])])
                    else:
                        txt = "\t".join(['None', str(barcode_counts['-']["Total"])])
                    bcFile.write(txt + '\n')

            # write out project table
            sys.stdout.write("%s reads (%s%% of total run) successfully identified barcode\n" % (bcsuccesscount, round((float(bcsuccesscount)/float(self.run.count()))*100, 1)))
            if evalPrimer:  # primer
                sys.stdout.write("%s reads (%s%% of total run) successfully identified barcode and primer\n" % (prsuccesscount, round((float(prsuccesscount)/float(self.run.count()))*100, 1)))
            if evalSample:  # sample
                sys.stdout.write("%s reads (%s%% of total run) successfully assigned to sample\n" % (sampsuccesscount, round((float(sampsuccesscount)/float(self.run.count()))*100, 1)))
            if minQ is not None:
                sys.stdout.write("%s reads (%s%% of total run) successfully pass trimming criteria\n" % (trimsuccesscount, round((float(trimsuccesscount)/float(self.run.count()))*100, 1)))

            sys.stdout.write("%s reads (%s%% of total run) unidentified\n\n" % (self.run.count()-identified_count, round((float(self.run.count()-identified_count)/float(self.run.count()))*100, 1)))

            if evalSample and self.verbose:
                for key in self.run_out:
                    sys.stdout.write("%s reads (%s%% of total run) found for project\t%s\n" % (self.run_out[key].count(), round((float(self.run_out[key].count())/float(self.run.count()))*100, 1), key))
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
