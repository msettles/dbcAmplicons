# tests/dbcAmplicons_test.py

import unittest
from dbcAmplicons import barcodeTable

class BarcodeTest(unittest.TestCase):
    def test_barcodes(self):
        d = barcodeTable("barcodeLookupTable.txt")
        self.assert_(len(d.barcodes) == 864)

if __name__ == "__main__":
    unittest.main()

