# barcodes.py
#
# barcode lookup file should look like
# #Barcode_name I1_barcode  I2_barcode
# Alpha1  TAAGGCGA TAGATCGC
# Alpha2  TAAGGCGA CTCTCTAT
# Alpha3  TAAGGCGA TATCCTCT
# Alpha4  TAAGGCGA AGAGTAGA
#

"""
barcodes.py parses and stores barcode information associated with a double barcoded illumina amplicon project
"""
import sys
from dbcAmplicons import misc


class barcodeTable:
    # ---------------- barcodes class ----------------
    """
    Class to read in and hold barcode table information associated with an Illumina double
    barcoded amplicon project
    """
    def __init__(self, barcodefile, i1_rc=True, i2_rc=False):
        """
        Initialize a new barcodeTable object with the file barcode table, parses and stores the barcode information
        """
        self.barcodes = {}
        self.bcPairI1 = {}
        self.bcPairI2 = {}
        self.I2 = []
        self.I1 = []
        self.IDS = []
        try:
            bcfile = open(barcodefile, 'r')
        except IOError:
            sys.stderr.write('ERROR:[Barcodes] Error cannot open %s\n' % barcodefile)
            raise
        f = bcfile.readlines()
        line = 0
        for row in f:
            line += 1
            if row[0] == "#" or row[0] == "\n":  # comment or blank line
                continue

            row = row.rstrip()
            try:
                #Test if file has 2 barcodes
                if len(row.split('\t')) == 3:
                    ID, I1BC, I2BC = row.split('\t')[0:3]
                    # I1 barcode shows up as the reverse complement in the sequencing run
                    if i1_rc:
                        I1BC = misc.reverseComplement(I1BC)
                    if i2_rc:
                        I2BC = misc.reverseComplement(I2BC)
                #Test if file has 1 barcode
                elif len(row.split('\t')) == 2:
                    ID, I1BC, = row.split('\t')[0:2]
                    I2BC = I1BC
                    # I1 barcode shows up as the reverse complement in the sequencing run
                    if i1_rc:
                        I1BC = misc.reverseComplement(I1BC)
                else:
                    print('ERROR: [Barcodes] File has an incorect number of columns')
            except ValueError as e:
                sys.stderr.write('ERROR:[Barcodes] Error reading line %s of barcode file: %s\n' % (str(line), str(e)))
                raise
            except KeyError:
                sys.stderr.write('ERROR:[Barcodes] Error reverse complementing I1 barcode %s, unknown character\n' % I1BC)
                raise
            except:
                sys.stderr.write('ERROR:[Barcodes] Unexpected error on line %s of the barcodes file: %s\n' % (line, sys.exc_info()[0]))
                raise
            #Test if file has 2 barcodes
            if len(row.split('\t')) == 3:
                if I2BC not in self.I2:
                    self.I2.extend([I2BC])
                if I1BC not in self.I1:
                    self.I1.extend([I1BC])
                self.IDS.extend([ID])
                self.barcodes["%s%s" % (I1BC, I2BC)] = [ID, 0, 0]
                self.bcPairI1[I1BC] = I2BC
                self.bcPairI2[I2BC] = I1BC
            #Test if file has 1 barcode
            elif len(row.split('\t')) == 2:
                self.I2.extend([I2BC])
                if I1BC not in self.I1:
                    self.I1.extend([I1BC])
                self.IDS.extend([ID])
                self.barcodes["%s%s" % (I1BC, I2BC)] = [ID, 0, 0]
                self.bcPairI1[I1BC] = I2BC
                self.bcPairI2[I2BC] = I1BC


        bcfile.close()

    def getLength(self):
        """
        get the length (number of barcodes) in the barcode table
        """
        return len(self.barcodes)

    def getI2(self):
        """
        get the I2 barcode sequences
        """
        return self.I2

    def getI1(self):
        """
        get the I1 barcode sequences
        """
        return self.I1

    def getBarcodes(self):
        """
        get the barcode pair ID available
        """
        return self.IDS

    def getMatch(self, bc1, bc2):
        """
        Determine if two barcodes have a matching barcode pair id, else return None
        """
        try:
            return(self.barcodes["%s%s" % (bc1, bc2)][0])
        except KeyError:
            return (None)

    def getMatchI2(self, bc2):
        """
        Determine if two barcodes have a matching barcode pair id, else return None
        """
        try:
            return(self.bcPairI2[bc2])
        except KeyError:
            return (None)

    def getMatchI1(self, bc1):
        """
        Determine if two barcodes have a matching barcode pair id, else return None
        """
        try:
            return(self.bcPairI1[bc1])
        except KeyError:
            return (None)
