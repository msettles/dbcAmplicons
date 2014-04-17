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

## sample sheet files should have at miniumum the 4 columns [SampleID,BarcodeID,PrimerPairID,ProjectID] order doesn't matter and should look something like
## SampleID TubeID  BarcodeID   PrimerPairID    Vol Conc    Quantity    ProjectID   Investigator
## 1    1   Hotel353    ITS-4_5 NA  NA  NA  Anahi-Pollen    Anahi
## 2    2   Hotel354    ITS-4_5 NA  NA  NA  Anahi-Pollen    Anahi
## 3    3   Hotel355    ITS-4_5 NA  NA  NA  Anahi-Pollen    Anahi

import re
import sys
class KeyFoundError(Exception):
    """
    Exception class KeyFoundError, to store barcode and primer that generate the key error
    """
    def __init__(self, barcode, primer):
        self.barcode = barcode
        self.primer = primer
    def getBarcode(self):
        return repr(self.barcode)
    def getPrimer(self):
        return repr(self.primer)

# ---------------- samples class ----------------
class sampleTable:
    """
    Class to read in and hold sample table information associated with an Illumina double
    barcoded amplicon project
    """
    def __init__(self, samplefile):
        """
        Initialize a new sampleTable object with the file sample table, parses and stores the sample information and 
        associated project. Class assumes the input file samplefile contains the following 4 columns
        'SampleID','BarcodeID','PrimerPairID','ProjectID' (defined in the header) others columns in the file are allowed and ignored
        """
        self.sampleCount = 0
        self.sampleTable = {}
        projects = []
        samples = []
        try:
            sfile = open(samplefile, 'r')
        except IOError:
            print 'ERROR:[Samples] cannot open', samplefile
            raise
        f = sfile.next() ## read the file
        header = f.rstrip()
        vheader = header.split('\t')
        try:
            sampleID_index = vheader.index("SampleID")
        except ValueError:
            print 'ERROR:[Samples] Column "SampleID" was not found in the samples table'
            raise
        try:
            barcodeID_index = vheader.index("BarcodeID")
        except ValueError:
            print 'ERROR:[Samples] Column "BarcodeID" was not found in the samples table'
            raise
        try:
            primerID_index = vheader.index("PrimerPairID")
        except ValueError:
            print 'ERROR:[Samples] Column "PrimerPairID" was not found in the samples table'
            raise
        try:
            projectID_index = vheader.index("ProjectID")
        except ValueError:
            print 'ERROR:[Samples] Column "ProjectID" was not found in the samples table'
            raise
        try:
            line=1
            for row in sfile:
                line += 1
                if row[0] == "#" or row[0] == "\n": # comment or blank line
                    continue
                row = row.rstrip() # strip off newline
                row = row.split('\t') # split by tab
                barcode = row[barcodeID_index]
                sid = re.sub(r'[\\/:"\' ]+', ".",row[sampleID_index]) # replace unsafe characters from sampleID with '.'
                pid = re.sub(r'[:"\'*?<>| ]+', ".",row[projectID_index]) # replace unsafe characters from projectID with '.'
                projects.append(pid)
                samples.append(sid)
                for primer in row[primerID_index].split(','):
                    primer = primer.strip()
                    if (primer == '*' or primer == '-') and barcode in self.sampleTable.keys():
                        raise KeyFoundError(barcode,primer)
                    elif barcode in self.sampleTable.keys():
                        if primer in self.sampleTable[barcode].keys():
                            raise KeyFoundError(barcode,primer)
                        else:
                            self.sampleTable[barcode][primer] = [sid,pid]
                    else:
                        self.sampleTable[barcode] = {}
                        self.sampleTable[barcode][primer] = [sid,pid]
                self.projectList = sorted(set(projects))
                self.samplesList = sorted(set(samples))
                self.sampleCount += 1
        except KeyFoundError as e:
            print 'ERROR:[Samples] Error in sample sheet: Barcode [%s] and Primer [%s] pair previously specified' % (e.getBarcode(), e.getPrimer())
            sfile.close()
            raise
        except:
            print 'ERROR:[Samples] Unexpected error on line %s of the samples file: %s' % (line,sys.exc_info()[0])
            sfile.close()
            raise
        sfile.close()
    def getSampleNumber(self):
        """
        Get the number of samples read in from the sampleFile
        """
        return repr(self.sampleCount)
    def getProjectList(self):
        """
        Get the list of projects the samples represent
        """
        return self.projectList
    def getSampleList(self):
        """
        Get the list of projects the samples represent
        """
        return self.samplesList
    def getSampleID(self,barcode,primer):
        """
        Given a barcode and primer, return the associated sampleID, "*" is allowed in the primer for 'any' primer match
        """
        try:
            sid = self.sampleTable[barcode]
            if sid.keys() == ['*'] and primer != None:
                return sid['*'][0]
            elif sid.keys() == ['-'] and primer == None:
                return sid['-'][0]
            else:
                return(sid[primer][0])
        except KeyError:
            return(None)
        except:
            print 'ERROR:[Samples] Unexpected error in getSampleID:', sys.exc_info()[0]
            raise
    def getProjectID(self,barcode,primer):
        """
        Given a barcode and primer, return the associated project, "*" is allowd in the primer for 'any' primer match
        """
        try:
            sid = self.sampleTable[barcode]
            if sid.keys() == ['*'] and primer != None:
                return sid['*'][1]
            elif sid.keys() == ['-'] and primer == None:
                print "Made it into dash", sid['-'][1]
                return sid['-'][1]
            else:
                return(sid[primer][1])
        except KeyError:
            return(None)
        except:
            print 'ERROR:[Samples] Unexpected error getProjectID:', sys.exc_info()[0]
            raise
