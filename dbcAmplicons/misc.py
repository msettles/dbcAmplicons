#### Misc functions

import os, errno
from subprocess import Popen, PIPE
import glob

'''
Gzip utilities, run gzip in a subprocess
'''
def sp_gzip_read(file):
    p = Popen(['gzip', '--decompress', '--to-stdout', file], stdout = PIPE, stderr = PIPE)
    return p.stdout

def sp_gzip_write(file):
    p = Popen('gzip > ' + file,stdin=PIPE,shell=True)
    return p.stdin


### Test Gzip file threaded gzip
#ifile = sp_gzip_read('Hkates_R1.fastq.gz')

#ofile = sp_gzip_write("Hkates_R1.2.fastq.gz")

#while 1:
#    line = ifile.readline()
#    ofile.write(line)


def infer_read_file_name(baseread, seakread):
    ''' Find other read filenames (ex. R1, R2, R3, R4) in the directory based on Read 1 filename '''
    basename = os.path.basename(baseread)
    path = os.path.dirname(os.path.realpath(baseread))
    testname = glob.glob(path + '/*' + os.path.splitext(baseread)[1])
    count = 0
    pos = -1
    read = ''
    for name in testname:
        count = 0
        if os.path.basename(name) == basename:
            continue
        elif len(os.path.basename(name)) != len(basename):
            continue
        else:
            for i, (ch1, ch2) in enumerate(zip(os.path.basename(name), basename)):
                if ch1 != ch2:
                    count += 1
                    if ch2 == '1' and ch1 == seakread:
                        pos = i
                    print str(i) + '\t' + ch1 + '\t' +  ch2
            if count == 1 and pos != -1:
                read = path + '/' + basename[0:pos] + seakread + basename[pos+1:]
                break
    return read


def make_sure_path_exists(path):
    """
    Try and create a path, if not error
    """
    if path != '':
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def reverseComplement(s):
    """
    given a seqeucne with 'A', 'C', 'T', and 'G' return the reverse complement
    """
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    letters = list(s)
    try:
        letters = [basecomplement[base] for base in letters]
    except:
        raise
    return ''.join(letters[::-1])

iupacdict = {
    'A':['A'],
    'C':['C'],
    'G':['G'],
    'T':['T'],
    'M':['A','C'],
    'R':['A','G'],
    'W':['A','T'],
    'S':['C','G'],
    'Y':['C','T'],
    'K':['G','T'],
    'V':['A','C','G'],
    'H':['A','C','T'],
    'D':['A','G','T'],
    'B':['C','G','T'],
    'X':['A','C','G','T'],
    'N':['A','C','G','T']}


def expand_iupac(seq):
    res = ['']
    for m in seq:
        newres = []
        for i in res:
            for j in iupacdict.get(m):
                newres.append(i + j)
        res = newres
    return res

#test = 'NNGATARNG'
#expand_iupac(test)
#test2 = 'AGGATAATG'
#expand_iupac(test2)
