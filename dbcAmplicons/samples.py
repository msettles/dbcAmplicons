# samples.py
#
# sample sheet files should have at miniumum the 4 columns [SampleID,BarcodeID,PrimerPairID,ProjectID] order doesn't matter and should look something like
# SampleID TubeID  BarcodeID   PrimerPairID    Vol Conc    Quantity    ProjectID   Investigator
# 1    1   Hotel353    ITS-4_5 NA  NA  NA  Anahi-Pollen    Anahi
# 2    2   Hotel354    ITS-4_5 NA  NA  NA  Anahi-Pollen    Anahi
# 3    3   Hotel355    ITS-4_5 NA  NA  NA  Anahi-Pollen    Anahi

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


class sampleTable:
    # ---------------- samples class ----------------
    """
    Class to read in and hold sample table information associated with an Illumina double
    barcoded amplicon project
    """
    def __init__(self, samplefile):
        """
        Initialize a new sampleTable object with the file sample table, parses and stores the sample
        information and associated project. Class assumes the input file samplefile contains the
        following 4 columns 'SampleID','BarcodeID','PrimerPairID','ProjectID' (defined in the header)
        others columns in the file are allowed and placed in Metadata
        """
        self.sampleCount = 0
        self.sampleTable = {}
        self.sampleMetadata = {}
        self.vhasMetadata = False
        projects = []
        samples = []
        try:
            sfile = open(samplefile, 'r')
        except IOError:
            sys.stderr.write('ERROR:[Samples] cannot open %s\n' % samplefile)
            raise
        f = sfile.next()  # read the file
        header = f.rstrip()
        vheader = header.split('\t')
        metadata_index = range(0, len(vheader))
        try:
            sampleID_index = vheader.index("SampleID")
            del metadata_index[sampleID_index]
        except ValueError:
            try:
                sampleID_index = vheader.index("#SampleID")
                del metadata_index[metadata_index.index(sampleID_index)]
            except ValueError:
                sys.stderr.write('ERROR:[Samples] Column "SampleID" was not found in the samples table\n')
                raise
        try:
            barcodeID_index = vheader.index("BarcodeID")
            del metadata_index[metadata_index.index(barcodeID_index)]
        except ValueError:
            sys.stderr.write('ERROR:[Samples] Column "BarcodeID" was not found in the samples table\n')
            raise
        try:
            primerID_index = vheader.index("PrimerPairID")
            del metadata_index[metadata_index.index(primerID_index)]
        except ValueError:
            sys.stderr.write('ERROR:[Samples] Column "PrimerPairID" was not found in the samples table\n')
            raise
        try:
            projectID_index = vheader.index("ProjectID")
            del metadata_index[metadata_index.index(projectID_index)]
        except ValueError:
            sys.stderr.write('ERROR:[Samples] Column "ProjectID" was not found in the samples table\n')
            raise
        metadata_headers = []
        if len(metadata_index) > 0:
            for index in metadata_index:
                metadata_headers.append(vheader[index])
        try:
            line = 1
            for row in sfile:
                line += 1
                if row[0] == "#" or row[0] == "\n":  # comment or blank line
                    continue
                row = row.rstrip()  # strip off newline
                row = row.split('\t')  # split by tab
                barcode = row[barcodeID_index]
                # sid = re.sub(r'[\\/:"\' ]+', ".", row[sampleID_index])  # replace unsafe characters from sampleID with '.'
                sid = re.sub(r'[^a-zA-Z0-9_-]', ".", row[sampleID_index])  # replace unsafe characters from sampleID with '.'
                pid = re.sub(r'[:"\'*?<>| ]+', ".", row[projectID_index])  # replace unsafe characters from projectID with '.'
                # pid = re.sub(r'[^a-zA-Z0-9_-]', ".", row[projectID_index])  # replace unsafe characters from projectID with '.'
                projects.append(pid)
                samples.append(sid)
                for primer in row[primerID_index].split(','):
                    primer = primer.strip()
                    if barcode in self.sampleTable.keys():
                        if primer in self.sampleTable[barcode].keys():
                            raise KeyFoundError(barcode, primer)
                        elif '*' in self.sampleTable[barcode].keys() and primer != '-':
                            raise KeyFoundError(barcode, primer)
                        else:
                            self.sampleTable[barcode][primer] = [sid, pid]
                            if sid not in self.sampleMetadata:
                                self.sampleMetadata[sid] = {}
                                self.sampleMetadata[sid][(sid, primer)] = {"PrimerPairID": primer, "BarcodeID": barcode, "SampleID": sid, "ProjectID": pid}
                            else:
                                self.sampleMetadata[sid][(sid, primer)] = {"PrimerPairID": primer, "BarcodeID": barcode, "SampleID": sid, "ProjectID": pid}
                    else:
                        self.sampleTable[barcode] = {}
                        self.sampleTable[barcode][primer] = [sid, pid]
                        if sid not in self.sampleMetadata:
                            self.sampleMetadata[sid] = {}
                            self.sampleMetadata[sid][(sid, primer)] = {"PrimerPairID": primer, "BarcodeID": barcode, "SampleID": sid, "ProjectID": pid}
                        else:
                            self.sampleMetadata[sid][(sid, primer)] = {"PrimerPairID": primer, "BarcodeID": barcode, "SampleID": sid, "ProjectID": pid}
                # parse metadata
                if len(metadata_index) > 0:
                    self.vhasMetadata = True
                    self.sampleMetadata[sid]["Metadata"] = {}
                    for i, index in enumerate(metadata_index):
                        self.sampleMetadata[sid]["Metadata"][metadata_headers[i]] = row[index]
                self.sampleCount += 1
            self.projectList = sorted(set(projects))
            self.samplesList = sorted(set(samples))
        except KeyFoundError as e:
            sys.stderr.write('ERROR:[Samples] Error in sample sheet: Barcode [%s] and Primer [%s] pair previously specified\n' % (e.getBarcode(), e.getPrimer()))
            sfile.close()
            raise
        except Exception:
            sys.stderr.write('ERROR:[Samples] Unexpected error on line %s of the samples file: %s\n' % (line, sys.exc_info()[0]))
            sfile.close()
            raise
        sfile.close()

    def hasMetadata(self):
        """
        whether or not the samplesheet also contains additional metadata
        """
        return self.vhasMetadata

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
        Get the list of samples
        """
        return self.samplesList

    def getSampleID(self, barcode, primer, override_primer=False):
        """
        Given a barcode and primer, return the associated sampleID, "*" is allowed in the primer for 'any' primer match
        """
        if barcode in self.sampleMetadata:
            # barcode, primer has already been assigned
            return self.sampleMetadata[barcode][(barcode, primer)]["SampleID"]
        try:
            sid = self.sampleTable[barcode]
            if primer is not None and '*' in sid.keys():
                return sid['*'][0]
            elif primer is None and '-' in sid.keys():
                return sid['-'][0]
            elif override_primer:
                if len(sid) > 1:
                    raise Exception('ERROR:[Samples]. override_primer was specified, and barcode %s is not unique in the sample table' % barcode)
                return(sid[sid.keys()[0]][0])
            else:
                return(sid[primer][0])
        except KeyError:
            return(None)
        except Exception:
            sys.stderr.write('ERROR:[Samples] Unexpected error in getSampleID: %s\n' % sys.exc_info()[0])
            raise

    def getProjectID(self, barcode, primer, override_primer=False):
        """
        Given a barcode and primer, return the associated project, "*" is allowed in the primer for 'any' primer match
        """
        # TODO Error occurs here when barcode id == sample id
        if barcode in self.sampleMetadata:
            # barcode, primer has already been assigned and barcode is sample_id
            return self.sampleMetadata[barcode][(barcode, primer)]["ProjectID"]
        try:
            sid = self.sampleTable[barcode]
            if primer is not None and '*' in sid.keys():
                return sid['*'][1]
            elif primer is None and '-' in sid.keys():
                return sid['-'][1]
            elif override_primer:
                if len(sid) > 1:
                    raise Exception('ERROR:[Samples]. override_primer was specified, and barcode %s is not unique in the sample table' % barcode)
                return(sid[sid.keys()[0]][1])
            else:
                return(sid[primer][1])
        except KeyError:
            return(None)
        except Exception:
            sys.stderr.write('ERROR:[Samples] Unexpected error getProjectID:\n', sys.exc_info()[0])
            raise
