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
import sys

def reverseComplement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = list(s)
    try:
        letters = [basecomplement[base] for base in letters]
    except:
        raise
    return ''.join(letters[::-1])

# ---------------- barcodes class ----------------
class barcodeTable:
    """ class to hold barcode information """
    def __init__(self, barcodefile):
        self.barcodes = {}
        self.P5 = []
        self.P7 = []
        self.IDS = []
        try:
            bcfile = open(barcodefile, 'r')
        except IOError:
            print '[Barcodes] Error cannot open', barcodefile
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
                P7BC = reverseComplement(P7BC)
            except ValueError as e:
                print "[Barcodes] Error reading line %s of barcode file: %s" % (str(line), str(e))
                raise
            except KeyError:
                print "[Barcodes] Error reverse complementing P7 barcode %s, unknown character" % P7BC
                raise
            except:
                print "[Barcodes] Unexpected error on line %s of the barcodes file: %s" % (line,sys.exc_info()[0])
                raise
            if P5BC not in self.P5:
                self.P5.extend([P5BC])
            if P7BC not in self.P7:
                self.P7.extend([P7BC])
            self.IDS.extend([ID])
            self.barcodes["%s%s" % (P7BC, P5BC)] = [ID,0,0]
        bcfile.close()
