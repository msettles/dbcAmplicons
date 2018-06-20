# primers.py
#
# primer lookup file should look like, where Read is P5 or R1 or READ1 and P7 or R2 or READ2,
# the '#' character represents a comments and will be ignored
# #Read    Pair_ID Primer_ID   Sequence
# P5  16S 27F_YM1 GTAGAGTTTGATCCTGGCTCAG
# P5  16S 27F_YM2 CGTAGAGTTTGATCATGGCTCAG
# P5  16S 27F_YM3 ACGTAGAGTTTGATTCTGGCTCAG
# P5  16S 27F_YM4 TACGTAGAGTTTGATTATGGCTCAG
# P5  16S 27F_Bif GTACGTAGGGTTCGATTCTGGCTCAG
# P5  16S 27F_Bor CGTACGTAGAGTTTGATCCTGGCTTAG

"""
primer.py parses and stores primer information associated with a double barcoded illumina amplicon project
"""
import sys
from dbcAmplicons import misc


class primerTable:
    # ---------------- primer class ----------------
    """
    Class to read in and hold amplicon pcr primer table information associated with an Illumina double
    barcoded amplicon project
    """
    def __init__(self, primerfile):
        """
        Initialize a new primerTable object with the file primer table, parses and stores the primer information
        """
        self.P5sequences = []
        self.P5id = {}
        self.P5pair = {}
        self.P7sequences = []
        self.P7id = {}
        self.P7pair = {}
        self.primers = []
        # TODO: add in check for presense of both P5 and P7 in pair
        try:
            prfile = open(primerfile, 'r')
        except IOError:
            sys.stderr.write('ERROR:[Primers] Error cannot open %s\n' % primerfile)
            raise
        f = prfile.readlines()
        line = 0
        for row in f:
            line += 1
            if row[0] == "#" or row[0] == "\n":  # comment or blank line
                continue
            row = row.rstrip()
            try:
                READ, PAIR, ID, SEQ = row.split('\t')[0:4]
            except ValueError as e:
                sys.stderr.write("ERROR:[Primers] Error reading line %s of primer file: %s\n" % (str(line), str(e)))
                raise
            except Exception:
                sys.stderr.write("ERROR:[Primers] Unexpected error on line %s of the primers file: %s\n" % (line, sys.exc_info()[0]))
                raise
            self.primers.extend([PAIR])
            pseqs = misc.expand_iupac(SEQ.upper())
            if READ in ["P5", "R1", "READ1", "F", "FORWARD"]:
                for pseq in pseqs:
                    if pseq in self.P5sequences:
                        self.P5id[pseq].extend([ID])
                        self.P5pair[pseq].extend([PAIR])
                    else:
                        self.P5sequences.extend([pseq])
                        self.P5id[pseq] = [ID]
                        self.P5pair[pseq] = [PAIR]
            if READ in ["P7", "R2", "READ2", "R", "REVERSE"]:
                for pseq in pseqs:
                    if pseq in self.P7sequences:
                        self.P7id[pseq].extend([ID])
                        self.P7pair[pseq].extend([PAIR])
                    else:
                        self.P7sequences.extend([pseq])
                        self.P7id[pseq] = [ID]
                        self.P7pair[pseq] = [PAIR]
        self.primers = sorted(list(set(self.primers)))

        prfile.close()

    def getPrimers(self):
        """
        Return the list of primer names
        """
        return self.primers

    def getP5sequences(self):
        """
        Return the list of P5 sequence
        """
        return self.P5sequences

    def getP7sequences(self):
        """
        Return the list of P7 sequences
        """
        return self.P7sequences

    def getMatch(self, seq1, seq2):
        """
        Determine if two primers are matching primers and return the
        primer pair name and two id names
        """
        if seq1 is not None and seq1 in self.P5pair.keys():
            pair1 = self.P5pair[seq1]
            id1 = self.P5id[seq1]
        else:
            pair1 = None
            id1 = [None]
        if seq2 is not None and seq2 in self.P7pair.keys():
            pair2 = self.P7pair[seq2]
            id2 = self.P7id[seq2]
        else:
            pair2 = None
            id2 = [None]
        # at least one primer not id
        if pair1 is None or pair2 is None:
            return [None, id1[0], id2[0]]
        # simple case of single length and matching
        elif len(pair1) == 1 and len(pair2) == 1 and pair1 == pair2:
            return [pair1[0], id1[0], id2[0]]
        # simple case of single length and not matching
        elif len(pair1) == 1 and len(pair2) == 1 and pair1 != pair2:
            return [None, id1[0], id2[0]]
        # at least one of two are > 1
        else:
            for i in range(len(pair1)):
                for j in range(len(pair2)):
                    if pair1[i] == pair2[j]:
                        return[pair1[i], id1[i], id2[j]]
        # Catch all
        return[None, id1, id2]
