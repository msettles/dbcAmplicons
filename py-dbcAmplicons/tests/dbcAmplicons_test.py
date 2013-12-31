# tests/dbcAmplicons_test.py

import unittest
from dbcAmplicons import barcodeTable
from dbcAmplicons import primerTable
from dbcAmplicons import editdist

class BarcodeTest(unittest.TestCase):
    def test_barcodes(self):
        d = barcodeTable("barcodeLookupTable.txt")
        self.assert_(len(d.barcodes) == 864)

class PrimerTest(unittest.TestCase):
    def test_primers(self):
        p = primerTable("primerLookupTable.txt")
        self.assert_(len(p.P5sequences) == 7 and len(p.P7sequences) ==  7)

class distanceTest(unittest.TestCase):
    def test_distance(self):
        d = editdist.distance("ABCDE","FBDG")
        self.assert_(d == 3)

class boundedDistanceTest(unittest.TestCase):
    def test_distance(self):
        d = editdist.bounded_distance("ABCDE","FBCDG",0,0)
        self.assert_(d[0] == 2 and d[1] == -1)

if __name__ == "__main__":
    unittest.main()

