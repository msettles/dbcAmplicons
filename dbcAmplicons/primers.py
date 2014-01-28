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

## primer lookup file should look like, where Read is P5 or R1 or READ1 and P7 or R2 or READ2, 
## the '#' character represents a comments and will be ignored
## #Read    Pair_ID Primer_ID   Sequence
## P5  16S 27F_YM1 GTAGAGTTTGATCCTGGCTCAG
## P5  16S 27F_YM2 CGTAGAGTTTGATCATGGCTCAG
## P5  16S 27F_YM3 ACGTAGAGTTTGATTCTGGCTCAG
## P5  16S 27F_YM4 TACGTAGAGTTTGATTATGGCTCAG
## P5  16S 27F_Bif GTACGTAGGGTTCGATTCTGGCTCAG
## P5  16S 27F_Bor CGTACGTAGAGTTTGATCCTGGCTTAG

"""
primer.py parses and stores primer information associated with a double barcoded illumina amplicon project
"""
import sys
# ---------------- primer class ----------------
class primerTable:
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
        try:
            prfile = open(primerfile, 'r')
        except IOError:
            print 'ERROR:[Primers] Error cannot open', primerfile
            raise
        f = prfile.readlines()
        line = 0
        for row in f:
            line += 1
            if row[0] == "#" or row[0] == "\n": # comment or blank line
                continue
            row = row.rstrip()
            try:
                READ, PAIR, ID, SEQ = row.split('\t')
            except ValueError as e:
                print "ERROR:[Primers] Error reading line %s of primer file: %s" % (str(line), str(e))
                raise
            except:
                print "ERROR:[Primers] Unexpected error on line %s of the primers file: %s" % (line,sys.exc_info()[0])
                raise
            self.primers.extend([PAIR])
            if READ in ["P5","R1","READ1"]:
                if SEQ in self.P5sequences:
                    self.P5id[SEQ].extend([ID])
                    self.P5pair[SEQ].extend([PAIR])
                else:
                    self.P5sequences.extend([SEQ])
                    self.P5id[SEQ] = [ID]
                    self.P5pair[SEQ] = [PAIR]
            if READ in ["P7","R2","READ2"]:
                if SEQ in self.P7sequences:
                    self.P7id[SEQ].extend([ID])
                    self.P7pair[SEQ].extend([PAIR])
                else:
                    self.P7sequences.extend([SEQ])
                    self.P7id[SEQ] = [ID]
                    self.P7pair[SEQ] = [PAIR]
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
    def getMatch(self,seq1,seq2):
        """
        Determine if two primers are matching primers and return the 
        primer pair name and two id names
        """
        if seq1 in self.P5pair.keys():
            pair1 = self.P5pair[seq1]
            id1 = self.P5id[seq1]
        else:
            pair1 = None
            id1 = [None]
        if seq2 in self.P7pair.keys():
            pair2 = self.P7pair[seq2]
            id2 = self.P7id[seq2]
        else:
            pair2 = None
            id2 = [None]
        # at least one primer not id
        if pair1 == None or pair2 == None:
            return [None,id1[0],id2[0]]
        # simple case of single length and matching
        elif len(pair1) == 1 and len(pair2) == 1 and pair1 == pair2:
            return [pair1[0],id1[0],id2[0]]
        # simple case of single length and not matching
        elif len(pair1) == 1 and len(pair2) == 1 and pair1 != pair2:
            return [None,id1[0],id2[0]]
        # at least one of two are > 1
        else:
            for i in range(len(pair1)):
                for j in range(len(pair2)):
                    if pair1[i] == pair2[j]:
                        return[pair1[i],id1[i],id2[j]]
        ## Catch all
        return[None,id1,id2]
