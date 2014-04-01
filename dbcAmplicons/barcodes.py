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

# ---------------- reverse complement a sequence s ----------------

## barcode lookup file should look like
## #Barcode_name   P5_barcode  P7_barcode
## Alpha1  TAGATCGC    TAAGGCGA
## Alpha2  CTCTCTAT    TAAGGCGA
## Alpha3  TATCCTCT    TAAGGCGA
## Alpha4  AGAGTAGA    TAAGGCGA
##

"""
barcodes.py parses and stores barcode information associated with a double barcoded illumina amplicon project
"""
import sys
from dbcAmplicons import misc
# ---------------- barcodes class ----------------
class barcodeTable:
    """
    Class to read in and hold barcode table information associated with an Illumina double
    barcoded amplicon project
    """
    def __init__(self, barcodefile):
        """
        Initialize a new barcodeTable object with the file barcode table, parses and stores the barcode information
        """
        self.barcodes = {}
        self.P5 = []
        self.P7 = []
        self.IDS = []
        try:
            bcfile = open(barcodefile, 'r')
        except IOError:
            print 'ERROR:[Barcodes] Error cannot open', barcodefile
            raise
        f = bcfile.readlines()
        line = 0
        for row in f:
            line += 1
            if row[0] == "#" or row[0] == "\n": # comment or blank line
                continue
            row = row.rstrip()
            try:
                ID, P5BC, P7BC = row.split('\t')
                # P7 barcode shows up as the reverse complement in the sequencing run
                P7BC = misc.reverseComplement(P7BC)
            except ValueError as e:
                print 'ERROR:[Barcodes] Error reading line %s of barcode file: %s' % (str(line), str(e))
                raise
            except KeyError:
                print 'ERROR:[Barcodes] Error reverse complementing P7 barcode %s, unknown character' % P7BC
                raise
            except:
                print 'ERROR:[Barcodes] Unexpected error on line %s of the barcodes file: %s' % (line,sys.exc_info()[0])
                raise
            if P5BC not in self.P5:
                self.P5.extend([P5BC])
            if P7BC not in self.P7:
                self.P7.extend([P7BC])
            self.IDS.extend([ID])
            self.barcodes["%s%s" % (P7BC, P5BC)] = [ID,0,0]
        bcfile.close()
    def getLength(self):
        """
        get the length (number of barcodes) in the barcode table
        """
        return len(self.barcodes)
    def getP5(self):
        """
        get the P5 barcode sequences
        """
        return self.P5
    def getP7(self):
        """
        get the P7 barcode sequences
        """
        return self.P7
    def getBarcodes(self):
        """
        get the barcode pair ID available
        """
        return self.IDS
    def getMatch(self,bc1,bc2):
        """ 
        Determine if two barcodes have a matching barcode pair id, else return None
        """
        try:
            return(self.barcodes["%s%s" % (bc1, bc2)][0])
        except KeyError:
           return (None)
