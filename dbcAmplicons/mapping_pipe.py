
import os
import shutil
import tempfile
from contextlib import contextmanager
from subprocess import Popen
from subprocess import PIPE
import signal
import shlex

from multiprocessing import Process

index = 'chloroplast_16S.fasta'
fastq1 = 'TestData10K_R1.fastq.gz'
fastq2 = 'TestData10K_R2.fastq.gz'


@contextmanager
def create_named_pipe():
    path = tempfile.mktemp()
    try:
        os.mkfifo(path)
        yield path
    finally:
        os.remove(path)


def sp_gzip_read(file):
    p = Popen(shlex.split('gzip --decompress --to-stdout') + [file],
        stdout=PIPE,
        stderr=None,
        bufsize=-1,
        preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    if p.returncode:
        raise
    return p.stdout


class fastqIter:
    " A simple file iterator that returns 4 lines for fast fastq iteration. "
    def __init__(self, handle):
        self.inf = handle
    def __iter__(self):
        return self
    def next(self):
        lines = {'id': self.inf.readline().strip(),
                 'seq': self.inf.readline().strip(),
                 '+': self.inf.readline().strip(),
                 'qual': self.inf.readline().strip()}
        assert(len(lines['seq']) == len(lines['qual']))
        if lines['id'] == '' or lines['seq'] == '' or lines['+'] == '' or lines['qual'] == '':
            raise StopIteration
        else:
            return lines
    @staticmethod
    def parse(handle):
        return fastqIter(handle)
    def close(self):
        self.inf.close()



def sp_bowtie2_read(indx, fq1, fq2):
    p = Popen(['bowtie2', '-x', indx, '-1', fq1, '-2', fq2],
              stdout=PIPE,
              stderr=PIPE,
              bufsize=-1,
              preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    if p.returncode:
        raise
    return (p.stdout, p.stderr)


def fastq_rename(fq1, fq2, out1, out2):
    iterator1 = fastqIter(sp_gzip_read(fq1))
    iterator2 = fastqIter(sp_gzip_read(fq2))
    try:
        while 1:
            seq1 = iterator1.next()
            seq2 = iterator2.next()
            split_name1 = seq1['id'].split(" ")
            split_name2 = seq2['id'].split(" ")
            if split_name1[0] == split_name2[0]:
                name = split_name1[0]
            else:
                raise Exception("read1 and read2 IDS do not match")
            sample = split_name1[1].split(":")[3]
            if (len(split_name1) == 4):
                primer = split_name1[1].split(":")[4]
            seq1['id'] = name + ":" + sample + ":" + primer + " 1:Y:0:"
            seq2['id'] = name + ":" + sample + ":" + primer + " 2:Y:0:"
            out1.write('\n'.join([seq1['id'],seq1['seq'],seq1['+'],seq1['qual']]))
            out2.write('\n'.join([seq2['id'],seq2['seq'],seq2['+'],seq2['qual']]))
    except StopIteration:
        pass
    finally:
        print "Finished processing one pair of files."


outR1 = create_named_pipe()
outR2 = create_named_pipe()

wpipe1 = open(outR1)
with create_named_pipe() as outR1:
    with create_named_pipe() as outR2:
        print outR1
        print outR2
        with open(outR1, 'wb') as wpipe1:
            with open(outR2, 'wb') as wpipe2:
                p = Process(target=fastq_rename, args=(fastq1, fastq2, wpipe1, wpipe2, ))
                p.start()
                sam_out = sp_bowtie2_read(index, outR1, outR2)
                output = sam_out[0]






outR2 = create_named_pipe()


def openFile(f):
    return(open(f, 'r'))


# Res[0] is the sam file output, Res[1] is the stderr bowtie2 'report'
Res = sp_bowtie2_read(index, fastq1, fastq2)



with named_pipe() as path:
    p = Popen(["pram_axdnull", str(kmer), input_filename, path],
              stdout=PIPE) # read from path
    with open(path, 'wb') as wpipe:
        kmer_proc = Popen(["generate_kmers", str(kmer)],
                          stdout=wpipe) # write to path
    output = p.communicate()[0]
    kmer_proc.wait()