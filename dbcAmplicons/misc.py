# Misc functions

import sys
import os
import errno
from subprocess import Popen
from subprocess import PIPE
import glob
import shlex

import re


'''
Gzip utilities, run gzip in a subprocess
'''


def sp_gzip_read(file):
    p = Popen(shlex.split('gzip --decompress --to-stdout') + [file], stdout=PIPE, stderr=None, bufsize=-1)
    return p.stdout


def sp_gzip_write(file):
    filep = open(file, 'wb')
    p = Popen('gzip', stdin=PIPE, stdout=filep, shell=True, bufsize=-1)
    return p.stdin


def parse_flash(fileinput_stream, verbose=True):
    skip = 4
    for i, line in enumerate(fileinput_stream):
        if skip == 4:
            # parse version
            sys.stdout.write('Using Flash_version:' + re.split(r' +', line.rstrip())[3] + '\n')
            skip = 3
            continue
        if skip == 3 and "Parameters" not in line:
            continue
        elif skip == 3 and "Parameters" in line:
            skip = 2
            continue
        elif skip == 2 and "Starting FASTQ readers and writer threads" not in line:
            # parse Parameters
            data = re.split(': +', re.sub(r'\[FLASH\] +', '', line.rstrip()))
            if len(data) == 2:
                name = re.sub(r' ', '_', data[0])
                if verbose:
                    sys.stdout.write(name + ':' + data[1] + '\n')
            continue
        elif skip == 2 and "Starting FASTQ readers and writer threads" in line:
            skip = 1
            continue
        elif skip == 1 and "Read combination statistics" not in line:
            sys.stderr.write(re.sub(r'\[FLASH\] +', '', line))
            continue
        elif skip == 1 and "Read combination statistics" in line:
            skip = 0
            continue
        elif skip == 0 and "Writing histogram files" not in line:
            # parse read combination statistics
            data = re.split(': +', re.sub(r'\[FLASH\] +', '', line.rstrip()))
            if len(data) == 2:
                name = re.sub(r' ', '_', data[0])
                sys.stdout.write(name + ':' + data[1] + '\n')
            continue
        elif skip == 0 and "Writing histogram files" in line:
            return(0)


def infer_read_file_name(baseread, seakread):
    ''' Find other read filenames (ex. R1, R2, R3, R4) in the directory based on Read 1 filename '''
    basename = os.path.basename(baseread)
    path = os.path.dirname(os.path.realpath(baseread))
    testname = glob.glob(path + '/*' + os.path.splitext(baseread)[1])
    count = 0
    pos = -1
    read = []
    for name in testname:
        count = 0
        if os.path.basename(name) == basename:  # ignore the same file
            continue
        elif len(os.path.basename(name)) != len(basename):  # must be the same length
            continue
        else:
            for i, (ch1, ch2) in enumerate(zip(os.path.basename(name), basename)):  # calculate the hamming distance
                if ch1 != ch2 and ch2 == '1' and ch1 == seakread:
                    count += 1
                    pos = i
            if count == 1:
                read.append(path + '/' + basename[0:pos] + seakread + basename[pos+1:])
                continue
    if len(read) == 1:
        return read[0]
    else:
        raise Exception("Error inferring read " + seakread + " from read 1, found " + str(len(read)) + " suitable matches.")


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
    return path


def expand_path(list):
    newlist = []
    for file in list:
        newlist.append(os.path.realpath(file))
    return newlist


def reverseComplement(s):
    """
    given a seqeucne with 'A', 'C', 'T', and 'G' return the reverse complement
    """
    basecomplement = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c', 'T': 'A', 't': 'a', 'N': 'N', 'n': 'n'}
    letters = list(s)
    try:
        letters = [basecomplement[base] for base in letters]
    except:
        raise
    return ''.join(letters[::-1])


def expand_iupac(seq):
    iupacdict = {
        'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'T': ['T'],
        'M': ['A', 'C'],
        'R': ['A', 'G'],
        'W': ['A', 'T'],
        'S': ['C', 'G'],
        'Y': ['C', 'T'],
        'K': ['G', 'T'],
        'V': ['A', 'C', 'G'],
        'H': ['A', 'C', 'T'],
        'D': ['A', 'G', 'T'],
        'B': ['C', 'G', 'T'],
        'X': ['A', 'C', 'G', 'T'],
        'N': ['A', 'C', 'G', 'T'],
        'a': ['a'],
        'c': ['c'],
        'g': ['g'],
        't': ['t'],
        'm': ['a', 'c'],
        'r': ['a', 'g'],
        'w': ['a', 't'],
        's': ['c', 'g'],
        'y': ['c', 't'],
        'k': ['g', 't'],
        'v': ['a', 'c', 'g'],
        'h': ['a', 'c', 't'],
        'd': ['a', 'g', 't'],
        'b': ['c', 'g', 't'],
        'x': ['a', 'c', 'g', 't'],
        'n': ['a', 'c', 'g', 't']}
    res = ['']
    for m in seq:
        newres = []
        for i in res:
            for j in iupacdict.get(m):
                newres.append(i + j)
        res = newres
    return res
