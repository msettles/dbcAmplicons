
import os
import sys
import tempfile
from contextlib import contextmanager
from subprocess import Popen
from subprocess import PIPE
import signal
import shlex


index = 'chloroplast_16S.fasta'
fastq1 = 'TestData10K_R1.fastq.gz'
fastq2 = 'TestData10K_R2.fastq.gz'

# Template
# bowtie2 -x caplanStuff -U <(zcat ../../CaplanShit/00-RawData/Sample_GCCAAT/GCCAAT_R1.fastq.gz| sed 's, ,_,g')


def sp_bowtie2_index(ref, overwrite=False):
    if os.path.isfile(ref):
        if os.path.isfile(ref + '.rev.2.bt2') and not overwrite:
            print 'Found existing bowtie2 index for %s' % ref
            return 0
        else:
            FNULL = open(os.devnull, 'w')
            call = 'bowtie2-build'
            call = call + ' ' + ref + ' ' + ref
            p = Popen(shlex.split(call),
                      stdout=FNULL,
                      stderr=FNULL,
                      bufsize=-1,
                      preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
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


def sp_bowtie2_screen(pe1, pe2, se, ref, overwrite=False):
    # build the call,
    # each file must first go through awk to replace spaces with a parsable character
    call = 'bowtie2'
    if sp_bowtie2_index(ref, overwrite) != 0:
        sys.exit(1)
    FNULL = open(os.devnull, 'w')
    call = call + ' -x ' + ref
    if (pe1 is not None and pe2 is not None):
        call = call + " -1 " + pe1 + " -2 " + pe2
    if (se is not None):
        call = call + " -U " + se
    p = Popen(shlex.split(call),
              stdout=PIPE,
              stderr=FNULL,
              bufsize=-1,
              preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    if p.returncode:
        raise
    return p.stdout


class screeningApp:

    def __init__(self):
        self.verbose = False

    def start(self, fastq_file1, fastq_file2, fastq_file3, reference, output_prefix, strict, procs, verbose=True, debug=False):
        """
            screen reads against a reference fasta file
        """
        self.verbose = verbose
        try:
            mapped_pairs_count = 0
            unmapped_pairs_count = 0
            mapped_singles_count = 0
            unmapped_singles_count = 0
            self.run_out = {}

            for project in sTable.getProjectList():
                self.run_out[project] = IlluminaTwoReadOutput(os.path.join(output_prefix, project), uncompressed)
            if output_unidentified:
                self.run_out["Unidentified"] = IlluminaTwoReadOutput(os.path.join(output_prefix, 'UnidentifiedProject'), uncompressed)

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
                    read.assignRead(sTable)  # barcode
                    if read.goodRead is True:
                        self.run_out[read.project].addRead(read.getFastq())
                        identified_count += 1
                    else:
                        unidentified_count += 1
                        if output_unidentified:
                            self.run_out["Unidentified"].addRead(read.getFastq())
                # Write out reads
                for key in self.run_out:
                    self.run_out[key].writeReads()
                if self.verbose:
                    sys.stderr.write("processed %s total reads, %s Reads/second, %s identified reads, %s unidentified reads (%s%%)\n" % (self.run.count(), round(self.run.count()/(time.time() - lasttime), 0), identified_count, unidentified_count, round((float(identified_count)/float(self.run.count()))*100)))
            if self.verbose:
                sys.stdout.write("%s reads processed in %s minutes, %s (%s%%) identified\n\n" % (self.run.count(), round((time.time()-lasttime)/(60), 2), identified_count, round((float(identified_count)/float(self.run.count()))*100, 1)))
            for key in self.run_out:
                sys.stdout.write("%s (%s%%)\treads found for project\t%s\n" % (self.run_out[key].count(), round((float(self.run_out[key].count())/float(self.run.count()))*100, 1), key))
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
