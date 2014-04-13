#### Misc functions

import os, errno

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
