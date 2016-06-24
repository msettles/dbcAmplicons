# classify_app.py
#
import sys
import os
import traceback
import time
from dbcAmplicons import TwoReadIlluminaRun
from dbcAmplicons import OneReadIlluminaRun
from dbcAmplicons import IlluminaFastaOutput

from subprocess import call
from multiprocessing import Pool


def rdpCall(query, output, gene, train, rdpPath, verbose):
    '''
    rdpCall takes a query fasta and generates the rdp call for the file, gene should be one of 16srrna or fungallsu
    rdpPath should point to the classifier.jar as a part of RDPTools
    '''
    if train is not None:
        rdp_call = ['java', '-Xmx2048M', '-Xms256M', '-XX:+UseParallelGC', '-XX:ParallelGCThreads=2', '-jar', rdpPath, 'classify', '-q', query, '-o', output, '-f', 'fixrank', '-t', train]
    else:
        if gene != "16srrna" and gene != "fungallsu" and gene != "fungalits_warcup" and gene != "fungalits_unite":
            sys.stderr.write('ERROR:[rdpCall] incorrect gene string provided to rdp\n')
            raise
        rdp_call = ['java', '-Xmx2048M', '-Xms256M', '-XX:+UseParallelGC', '-XX:ParallelGCThreads=2', '-jar', rdpPath, 'classify', '-q', query, '-o', output, '-f', 'fixrank', '-g', gene]
    starttime = time.time()
    if verbose:
        sys.stderr.write("Starting rdp for file %s\n" % query)
    res = call(rdp_call)
    if res == 0:
        try:
            os.remove(query)
        except OSError, e:  # if failed, report it back to the user ##
            sys.stderr.write("ERROR:[rdpCall] %s - %s.\n" % (e.filename, e.strerror))
            raise
        if verbose:
            sys.stderr.write("Finished processing %s in %s minutes\n" % (query, round((time.time() - starttime) / 60, 2)))
        return res
    else:
        sys.stderr.write("ERROR:[rdpCall] RDP did not finish properly, returned: %s\n" % res)
        raise


def check_status(results):
    """
    Checks the stats of a dictionary full of jobs and returns the number of jobs which
    have not finished.
    """
    unfinished = 0
    finished = 0
    for r in results:
        if not results[r].ready():
            unfinished += 1
        else:
            finished += 1
    return unfinished


class classifyApp:
    """
    Classify preprocessed Illumina Reads using RDP
    Takes fastq files and outputs batched fasta files that are each processed by rdp in parrallel
    """
    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, fastq_fileU, output_prefix, rdpPath, gene='16srrna', train=None, batchsize=10000, minQ=None, minL=0, procs=1, test=False, verbose=True, debug=False):
        """
        Start classifying double barcoded Illumina sequencing run
        """
        results = {}
        self.verbose = verbose
        try:
            if (train is None and gene != '16srrna' and gene != 'fungallsu' and gene != "fungalits_warcup" and gene != "fungalits_unite"):
                sys.stderr.write("ERROR:[classify] parameter -g (--gene) must be one of 16srrna or fungallsu or fungalits_warcup or fungalits_unite \n")
                raise Exception
            # establish and open the Illumina run
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
                sys.stderr.write("ERROR:[classify] input reads not specified, or incorrect pairs\n")
                raise Exception
            lasttime = time.time()
            batch = 0
            pool = Pool(procs, maxtasksperchild=1)
            if (self.runSingle is not None):
                while 1:
                    # get next batch of reads
                    reads = self.runSingle.next(batchsize)
                    batch = batch + len(reads)
                    if len(reads) == 0:
                        break
                    run_out = IlluminaFastaOutput(output_prefix + "." + str(batch))
                    # process individual reads
                    for read in reads:
                        run_out.addRead(read.getFasta())
                    # Write out reads
                    run_out.writeReads()
                    rdp_out = output_prefix + "." + str(batch) + ".fixrank"
                    results[rdp_out] = pool.apply_async(rdpCall, (run_out.output_prefix, rdp_out, gene, train, rdpPath, self.verbose))
                    if test:
                        break
            if (self.runPairs is not None):
                while 1:
                    # get next batch of reads
                    reads = self.runPairs.next(batchsize)
                    batch = batch + len(reads)
                    if len(reads) == 0:
                        break
                    run_out = IlluminaFastaOutput(output_prefix + "." + str(batch))
                    # process individual reads
                    for read in reads:
                        # if minQ != None: # TODO: add trim read to TwoSequenceReadSet
                            # read.trimRead(minQ, minL)
                        # if read.goodRead == True:
                        run_out.addRead(read.getJoinedFasta())
                    # Write out reads
                    run_out.writeReads()
                    rdp_out = output_prefix + "." + str(batch) + ".fixrank"
                    results[rdp_out] = pool.apply_async(rdpCall, (run_out.output_prefix, rdp_out, gene, train, rdpPath, self.verbose))
                    if test:
                        break

            allfinished = False
            while not allfinished:
                time.sleep(1)
                np = check_status(results)
                if np == 0:
                    allfinished = True
            if self.verbose:
                sys.stderr.write("Combining temporary files\n")
            with open(output_prefix + ".fixrank", "wb") as outfile:
                for f in results.keys():
                    with open(f, "rb") as infile:
                        outfile.write(infile.read())
                    os.remove(f)
            if self.verbose:
                sys.stdout.write("%s reads processed in %s minutes\n" % (batch, round((time.time() - lasttime) / (60), 2)))
            self.clean(results)
            return 0
        except (KeyboardInterrupt, SystemExit):
            self.clean(results)
            sys.stderr.write("%s unexpectedly terminated\n" % (__name__))
            return 1
        except:
            self.clean(results)
            sys.stderr.write("A fatal error was encountered.\n")
            if debug:
                sys.stderr.write("".join(traceback.format_exception(*sys.exc_info())))
            return 1

    def clean(self, results):
        if self.verbose:
            sys.stderr.write("Cleaning up.\n")
        try:
            self.runSingle.close()
            self.runPairs.close()
            for f in results.keys():
                os.remove(f)
        except:
            pass
