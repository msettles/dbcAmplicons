# sequenceReads.py
# sequenceReads.py stores and processes individual DNA sequence reads.

import sys
from dbcAmplicons import misc

try:
    from dbcAmplicons import editdist
    editdist_loaded = True
except ImportError:
    sys.stderr.write("Warning: editdist library not loaded, Insertion/Deletion detection in barcodes and primers will not be performed\n")
    editdist_loaded = False

try:
    from dbcAmplicons import trim
    trim_loaded = True
except ImportError:
    sys.stderr.write("Warning: trim library not loaded, trimming using python\n")
    trim_loaded = False


def barcodeDist(b_l, b_2, max_diff):
    # ------------------- calculate distance between two barcode sequences ------------------------------
    """
    gets edit distance between two equal-length barcode strings
    """
    bc = None
    bc_mismatch = max_diff + 1

    if editdist_loaded:
        bc_i, bc_mismatch = editdist.hamming_distance_list(b_l, b_2, max_diff + 1)

    else:
        for i, key in enumerate(b_l):
            if len(key) == len(b_2) and len(b_2) > 0:
                dist = sum(map(lambda x: x[0] != x[1], zip(key, b_2)))
            else:
                sys.stderr.write("ERROR:[barcodeDist] lengths of barcodes and index read do not match!\n")
                sys.stderr.write("Target: %s" % key)
                sys.stderr.write("Index read: %s" % b_2)
                raise Exception

            if dist < bc_mismatch:
                bc = i
                bc_mismatch = dist
    if bc_mismatch > max_diff:
        bc = None
    else:
        bc = b_l[bc_i]
    return (bc, bc_mismatch)


def primerDist(primer_l, read, dedup_float, max_diff, end_match):
    # ------------------- calculate distance between primer sequence and first part of read ------------------------------
    """
    gets edit distance between primer and read sequence
    """
    pr = None
    prMismatch = max_diff + 1
    prStartPosition = 0
    prEndPosition = 0

    if editdist_loaded:
        pr_i, prMismatch, prStartPosition, prEndPosition = editdist.bounded_distance_list(primer_l, read, dedup_float, max_diff, end_match)

    else:
        for i, key in enumerate(primer_l):
            read = read[0:len(key)]
            if end_match == 0:
                dist = sum(map(lambda x: x[0] != x[1], zip(read, key)))
            else:
                dist = sum(map(lambda x: x[0] != x[1], zip(read[:-end_match], key[:-end_match])))
                for i in range(len(key) - end_match, len(key)):
                    if read[i] != key[i]:
                        dist += 100
            if dist < prMismatch:
                pr_i = i
                prMismatch = dist
                prEndPosition = len(key)

    if prMismatch > max_diff:
        pr = None
        prStartPosition = 0
        prEndPosition = 0
    else:
        pr = primer_l[pr_i]

    return (pr, prMismatch, prStartPosition, prEndPosition)


class FourSequenceReadSet:
    # ---------------- Class for 4 read sequence data with double barcode set ----------------
    """
    Class to hold one Illumina four read (double barocodes) set, only one read names is used for all sequence reads
    Class processes a read by defining the barcode pair and primer pair [optional] and project [optional]. Finally
    class returns a paired read set for output
    """

    def __init__(self, name, read_1, qual_1, read_2, qual_2, bc_1, bc_2):
        """
        Initialize a FourSequenceReadSet with a name, two read sequences
        and cooresponding quality sequence and two barcode sequences. A read
        is initially defined as 'not' a good read and requires processing before
        being labeled as a good read.
        """
        self.barcode = [None, 0, 0]  # when filled, a vector of length 3 [barcodePairID, editdist1, editdist2]
        self.primer = [None, None, 0, 0, 0, None, 0, 0, 0]  # when filled, a vector of length 7 [primerPairID,primerID1,editdist1,readendpos1, primerID2,editdist2,readendpos2]
        self.sample = None
        self.project = None
        self.name = name
        self.read_1 = read_1
        self.qual_1 = qual_1
        self.read_2 = read_2
        self.qual_2 = qual_2
        self.bc_1 = bc_1
        self.bc_2 = bc_2
        self.trim_left = len(read_1)
        self.trim_right = len(read_2)
        self.goodRead = False

    def assignBarcode(self, bcTable, max_diff):
        """
        Given a barcodeTable object and the maximum number of allowed difference (mismatch, insertion, deletions)
        assign a barcode pair ID from the reads barcodes.
        """
        # Barcode One Matching
        bc1, bc1Mismatch = barcodeDist(bcTable.getI1(), self.bc_1, max_diff)

        # Barcode Two Matching
        bc2, bc2Mismatch = barcodeDist(bcTable.getI2(), self.bc_2, max_diff)

        # Barcode Pair Matching
        self.barcode = [bcTable.getMatch(bc1, bc2), bc1Mismatch, bc2Mismatch]
        self.sample = self.barcode[0]
        self.goodRead = self.barcode[0] is not None
        if self.goodRead:
            return 1
        else:
            return 0

    def assignPrimer(self, prTable, dedup_float, max_diff, endmatch):
        """
        Given a primerTable object, the maximum number of allowed difference (mismatch, insertion, deletions) and the
        required number of end match bases (final endmatch bases must match) assign a primer pair ID from the read
        sequences.
        """
        # Primer One Matching
        pr1, pr1Mismatch, pr1StartPosition, pr1EndPosition = primerDist(prTable.getP5sequences(), self.read_1, dedup_float, max_diff, endmatch)

        # Primer One Matching
        pr2, pr2Mismatch, pr2StartPosition, pr2EndPosition = primerDist(prTable.getP7sequences(), self.read_2, dedup_float, max_diff, endmatch)

        # Barcode Pair Matching
        combined_pr = prTable.getMatch(pr1, pr2)
        self.primer = [combined_pr[0], combined_pr[1], pr1Mismatch, pr1StartPosition, pr1EndPosition, combined_pr[2], pr2Mismatch, pr2StartPosition, pr2EndPosition]
        self.goodRead = self.goodRead and combined_pr[0] is not None
        if self.goodRead:
            return 1
        else:
            return 0

    def assignRead(self, sTable):
        """
        Given a samplesTable object, assign a sample ID and project ID using the reads barcode and primer designation
        """
        self.project = sTable.getProjectID(self.barcode[0], self.primer[0])
        self.sample = sTable.getSampleID(self.barcode[0], self.primer[0])
        self.goodRead = self.project is not None
        if self.goodRead:
            return 1
        else:
            return 0

    def trimRead(self, minQ, minL):
        """
        Trim the read by a minQ score
        """
        if (trim_loaded):
            trim_points = trim.trim(self.qual_1, self.qual_2, minQ)
            self.trim_left = trim_points["left_trim"]
            self.trim_right = trim_points["right_trim"]
            if ((self.trim_left - self.primer[4]) < minL or (self.trim_right - self.primer[8]) < minL):
                self.goodRead = False
        if self.goodRead:
            return 1
        else:
            return 0

    def getBarcode(self):
        """
        Return the reads barcode pair ID
        """
        return self.barcode[0]

    def getPrimer(self):
        """
        Return the reads primer pair ID
        """
        return self.primer[0]

    def getSampleID(self):
        """
        Return the reads sample ID
        """
        return self.sample

    def getProject(self):
        """
        Return the reads project ID
        """
        return self.project

    def getFastq(self, kprimer):
        """
        Create four line string ('\n' separator included) for the read pair, returning a length 2 vector (one for each read)
        """
        if self.primer[0] is not None and kprimer is False:
            read1_name = "%s 1:N:0:%s:%s %s|%s|%s|%s %s|%s|%s|%s" % (self.name, self.sample, self.primer[0], self.bc_1, self.barcode[1], self.bc_2, self.barcode[2], self.primer[1], self.primer[2], self.primer[4], self.read_1[0:self.primer[3]])
            read2_name = "%s 2:N:0:%s:%s %s|%s|%s|%s %s|%s|%s|%s" % (self.name, self.sample, self.primer[0], self.bc_1, self.barcode[1], self.bc_2, self.barcode[2], self.primer[5], self.primer[6], self.primer[8], self.read_2[0:self.primer[7]])
            r1 = '\n'.join([read1_name, self.read_1[self.primer[4]:self.trim_left], '+', self.qual_1[self.primer[4]:self.trim_left]])
            r2 = '\n'.join([read2_name, self.read_2[self.primer[8]:self.trim_right], '+', self.qual_2[self.primer[8]:self.trim_right]])
        else:
            read1_name = "%s 1:N:0:%s %s|%s|%s|%s" % (self.name, self.sample, self.bc_1, self.barcode[1], self.bc_2, self.barcode[2])
            read2_name = "%s 2:N:0:%s %s|%s|%s|%s" % (self.name, self.sample, self.bc_1, self.barcode[1], self.bc_2, self.barcode[2])
            r1 = '\n'.join([read1_name, self.read_1[0:self.trim_left], '+', self.qual_1[0:self.trim_left]])
            r2 = '\n'.join([read2_name, self.read_2[0:self.trim_right], '+', self.qual_2[0:self.trim_right]])
        return [r1, r2]


# ---------------- Class for 2 read sequence data processed with dbcAmplicons preprocess ----------------
class TwoSequenceReadSet:
    """
    Class to hold one Illumina two read set, assumed to have already been preprocessed with Barcodes and Primers already
    identified. Class processes a read by defining sample and project ids. Finally class returns a paired read set for output.
    """
    def __init__(self, name_1, read_1, qual_1, name_2, read_2, qual_2):
        """
        Initialize a TwoSequenceReadSet with names, two read sequences and cooresponding quality sequence.
        Barcode and primer sequence is inferred by their placement in the read names.
        A read is initially defined as 'not' a good read and requires processing before being labeled as a good read.
        """
        try:
            split_name = name_1.split(" ")
            self.name = split_name[0]
            self.barcode = split_name[1].split(":")[3]
            self.sample = self.barcode
            # dbcAmplicons processed, Primer and Barcode Identified
            if (len(split_name) == 4):
                self.primer_string1 = split_name[3]
                self.primer_string2 = name_2.split(" ")[3]
                self.primer = split_name[1].split(":")[4]
                self.barcode_string = split_name[2]
            # dbcAmplicons processed, Barcode only Identified
            elif (len(split_name) == 3):
                self.primer_string1 = None
                self.primer_string2 = None
                self.primer = None
                self.barcode_string = split_name[2]
            # Casava processed, check for barcode in the second string
            else:
                self.primer_string1 = None
                self.primer_string2 = None
                self.primer = None
                self.barcode_string = None
            self.read_1 = read_1
            self.qual_1 = qual_1
            self.read_2 = read_2
            self.qual_2 = qual_2
            self.trim_left = len(read_1)
            self.trim_right = len(read_2)
            self.goodRead = False
        except IndexError:
            sys.stderr.write('ERROR:[TwoSequenceReadSet] Read names are not formatted in the expected manner\n')
            raise
        except:
            sys.stderr.write('ERROR:[TwoSequenceReadSet] Unknown error occured initiating read\n')
            raise

    def assignRead(self, sTable):
        """
        Given a samplesTable object, assign a sample ID and project ID using the reads barcode and primer designation
        """
        self.project = sTable.getProjectID(self.barcode, self.primer)
        self.sample = sTable.getSampleID(self.barcode, self.primer)
        self.goodRead = self.project is not None
        return 0

    def trimRead(self, minQ, minL):
        """
        Trim the read by a minQ score
        """
        if (trim_loaded):
            trim_points = trim.trim(self.qual_1, self.qual_2, minQ)
            self.trim_left = trim_points["left_trim"]
            self.trim_right = trim_points["right_trim"]
            if (self.trim_left < minL or self.trim_right < minL):
                self.goodRead = False
            else:
                self.goodRead = True

    def getFastqSRA(self):
        """
        Create four line string ('\n' separator included) for the read pair, returning a length 2 vector (one for each read)
        """
        read1_name = "%s 1:N:0:" % (self.name)
        read2_name = "%s 2:N:0:" % (self.name)
        r1 = '\n'.join([read1_name, self.read_1[0:self.trim_left], '+', self.qual_1[0:self.trim_left]])
        r2 = '\n'.join([read2_name, self.read_2[0:self.trim_right], '+', self.qual_2[0:self.trim_right]])
        return [r1, r2]

    def getFastq(self):
        """
        Create four line string ('\n' separator included) for the read pair, returning a length 2 vector (one for each read)
        """
        if self.primer is not None:
            read1_name = "%s 1:N:0:%s:%s %s %s" % (self.name, self.sample, self.primer, self.barcode_string, self.primer_string1)
            read2_name = "%s 2:N:0:%s:%s %s %s" % (self.name, self.sample, self.primer, self.barcode_string, self.primer_string2)
        else:
            read1_name = "%s 1:N:0:%s %s" % (self.name, self.sample, self.barcode_string)
            read2_name = "%s 2:N:0:%s %s" % (self.name, self.sample, self.barcode_string)
        r1 = '\n'.join([read1_name, self.read_1[0:self.trim_left], '+', self.qual_1[0:self.trim_left]])
        r2 = '\n'.join([read2_name, self.read_2[0:self.trim_right], '+', self.qual_2[0:self.trim_right]])
        return [r1, r2]

    def getFourReads(self, bc1_length=8, bc2_length=8):
        """
        Create four line string ('\n' separator included) for the read set, returning a length 4 vector (one for each read)
        """
        try:
            if self.barcode_string is not None:
                bc1 = self.barcode_string.split('|')[0]
                bc2 = self.barcode_string.split('|')[2]
            elif len(self.barcode) is (bc1_length + bc2_length):
                bc1 = self.barcode[:bc1_length]
                bc2 = self.barcode[bc1_length:]
            else:
                raise Exception("string in the barcode is not %s characters" % str(bc1_length + bc2_length))
            r1 = '\n'.join([self.name + ' 1:Y:0:', self.read_1[0:self.trim_left], '+', self.qual_1[0:self.trim_left]])
            r2 = '\n'.join([self.name + ' 2:Y:0:', bc1, '+', 'C' * len(bc1)])  # Give barcodes and arbitary quality of Cs
            r3 = '\n'.join([self.name + ' 3:Y:0:', bc2, '+', 'C' * len(bc2)])
            r4 = '\n'.join([self.name + ' 4:Y:0:', self.read_2[0:self.trim_right], '+', self.qual_2[0:self.trim_right]])
            return [r1, r2, r3, r4]
        except IndexError:
            sys.stderr.write('ERROR:[TwoSequenceReadSet] unable to exract barocode sequence from the read names\n')
            raise
        except:
            sys.stderr.write('ERROR:[TwoSequenceReadSet] Unknown error occured generating four read set\n')
            raise

    def getFourReadsInline(self, bcTable, bc1_length=8, bc2_length=0, max_diff=1, flip_float=flip_float):
        """
        Create four line string ('\n' separator included) for the read set, returning a length 4 vector (one for each read)
        """
        """
        Given a barcodeTable object and the maximum number of allowed difference (mismatch, insertion, deletions)
        assign a barcode pair ID from the reads barcodes.
        """
        flip_float = False
        # Barcode One Matching
        if bc1_length > 0:
            bc1, bc1Mismatch = barcodeDist(bcTable.getI1(), self.read_1[:bc1_length], max_diff)
            if flip_float:
                bc1f, bc1fMismatch = barcodeDist(bcTable.getI2(), self.read_2[:bc1_length], max_diff)
                if bc1fMismatch < bc1Mismatch:  # flip_float produces a better match reverse the reads
                    flip_float = True
                    bc1 = bc1f
                    bc1Mismatch = bc1fMismatch
                    # Barcode Two Matching (with barcoce 1 present)
                    tmp = self.read_1
                    self.read_1 = self.read_2
                    self.read_2 = tmp
                    tmp = self.qual_1
                    self.qual_1 = self.qual_2
                    self.qual_2 = tmp
            if bc2_length > 0:
                bc2, bc2Mismatch = barcodeDist(bcTable.getI2(), self.read_2[:bc2_length], max_diff)
            else:
                bc2, bc2Mismatch = (bcTable.getMatchI1(bc1), 0)
        # Barcode Two Matching (when no barcode 1 present)
        elif bc2_length > 0:
            bc2, bc2Mismatch = barcodeDist(bcTable.getI2(), self.read_2[:bc2_length], max_diff)
            if flip_float:
                bc2f, bc2fMismatch = barcodeDist(bcTable.getI2(), self.read_1[:bc2_length], max_diff)
                if bc2fMismatch < bc2Mismatch:  # flip_float produces a better match rever the reads
                    flip_float = True
                    bc2 = bc2f
                    bc2Mismatch = bc2fMismatch
                    tmp = self.read_1
                    self.read_1 = self.read_2
                    self.read_2 = tmp
                    tmp = self.qual_1
                    self.qual_1 = self.qual_2
                    self.qual_2 = tmp
            bc1, bc1Mismatch = (bcTable.getMatchI2(bc2), 0)
        else:
            sys.stderr.write('ERROR:[TwoSequenceReadSet] barcodes length both 0 in four read inline bc set\n')
            raise

        # check if bc found
        if bc1 is None or bc2 is None:
            return []

        r1 = '\n'.join([self.name + ' 1:Y:0:', self.read_1[bc1_length:], '+', self.qual_1[bc1_length:]])
        r2 = '\n'.join([self.name + ' 2:Y:0:', bc1, '+', 'C' * len(bc1)])  # Give barcodes and arbitary quality of Cs
        r3 = '\n'.join([self.name + ' 3:Y:0:', bc2, '+', 'C' * len(bc2)])
        r4 = '\n'.join([self.name + ' 4:Y:0:', self.read_2[bc2_length:], '+', self.qual_2[bc2_length:]])
        return [r1, r2, r3, r4, flip_float]

    def getFasta(self):
        """
        Create two line string ('\n' separator included) for the read pair, returning a length 1 vector (one read)
        """
        name = '>' + self.name[1:]
        if self.primer is not None:
            read1_name = "%s 1:N:0:%s:%s:%i" % (name, self.sample, self.primer, len(self.read_1))
            read2_name = "%s 2:N:0:%s:%s:%i" % (name, self.sample, self.primer, len(self.read_2))
        else:
            read1_name = "%s 1:N:0:%s:%i" % (name, self.sample, len(self.read_1))
            read2_name = "%s 2:N:0:%s:%i" % (name, self.sample, len(self.read_2))
        r1 = '\n'.join([read1_name, self.read_1[0:self.trim_left]])
        r2 = '\n'.join([read2_name, self.read_2[0:self.trim_right]])
        return [r1, r2]

    def getJoinedFasta(self):
        """
        Create two line string ('\n' separator included) for the read pair, concatenating the two reads into a single returning length 1 vector (one read)
        """
        name = '>' + self.name[1:]
        joined_seq = self.read_1[0:self.trim_left] + misc.reverseComplement(self.read_2[0:self.trim_right])
        if self.primer is not None:
            read1_name = "%s|%s:%s:PAIR" % (name, self.sample, self.primer)
        else:
            read1_name = "%s|%s:PAIR" % (name, self.sample)
        r1 = '\n'.join([read1_name, joined_seq])
        return [r1]


# ---------------- Class for 2 read sequence data processed with dbcAmplicons preprocess ----------------
class OneSequenceReadSet:
    """
    Class to hold a one Illumina read set, assumes the paired reads produced by dbcAmplicons preprocess have been merged
    """
    def __init__(self, name_1, read_1, qual_1):
        """
        Initialize a OneSequenceReadSet with name, one read sequences and cooresponding quality sequence.
        A read is initially defined as 'not' a good read and requires processing before being labeled as a good read.
        """
        self.goodRead = False
        self.primer = None
        # parse the read name for the sample and primer ids
        try:
            split_name = name_1.split(" ")
            self.name = split_name[0]
            self.sample = split_name[1].split(":")[3]
            if (len(split_name) == 4):
                self.primer = split_name[1].split(":")[4]
        except IndexError:
            sys.stderr.write('ERROR:[OneSequenceReadSet] Read names are not formatted in the expected manner\n')
            raise
        except:
            sys.stderr.write('ERROR:[OneSequenceReadSet] Unknown error occured initiating read\n')
            raise
        self.read_1 = read_1
        self.qual_1 = qual_1

    def getFastqSRA(self):
        """
        Create four line string ('\n' separator included) for the read pair, returning a length 2 vector (one for each read)
        """
        read1_name = "%s 1:N:0:" % (self.name)
        r1 = '\n'.join([read1_name, self.read_1, '+', self.qual_1])
        return [r1]

    def getFastq(self):
        """
        Create four line string ('\n' separator included) for the read, returning a length 1 vector (one read)
        """
        if self.primer is not None:
            read1_name = "%s 1:N:0:%s:%s" % (self.name, self.sample, self.primer)
        else:
            read1_name = "%s 1:N:0:%s" % (self.name, self.sample)
        r1 = '\n'.join([read1_name, self.read_1, '+', self.qual_1])
        return [r1]

    def getFasta(self):
        """
        Create two line string ('\n' separator included) for the read, returning a length 1 vector (one read)
        """
        name = '>' + self.name[1:]
        if self.primer is not None:
            read1_name = "%s|%s:%s:%i" % (name, self.sample, self.primer, len(self.read_1))
        else:
            read1_name = "%s|%s:%i" % (name, self.sample, len(self.read_1))
        r1 = '\n'.join([read1_name, self.read_1])
        return [r1]

    def trimRead(self, minQ, minL):
        """
        Trim the read by a minQ score
        """
        if (trim_loaded):
            trim_points = trim.trim(self.qual_1, self.qual_1, minQ)
            self.trim_left = trim_points["left_trim"]
            if (self.trim_left < minL):
                self.goodRead = False
            else:
                self.goodRead = True
