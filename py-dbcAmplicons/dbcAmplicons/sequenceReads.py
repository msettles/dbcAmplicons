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

import sys
from dbcAmplicons import editdist
#try: 
#	import editdist
editdist_loaded = True
#except ImportError: 
#	editdist_loaded = False 


#------------------- calculate distance between two barcode sequences ------------------------------
def barcodeDist(b_1, b_2):
	'counts mismatches between two equal-length strings'
	if editdist_loaded :
		#print "Using editdist!"
		return editdist.distance(b_1,b_2)
	elif len(b_1) == len(b_2) and len(b_1) > 0:
		return sum(map(lambda x: x[0] != x[1], zip(b_1, b_2) ))
	else:
		print "ERROR lengths of barcodes and index read do not match!"
		print "Target", b_1
		print "Index read:", b_2
		sys.exit()



# ---------------- Class for 4 read sequence data with double barcode set ----------------
class SequenceReadSet:
	""" class to hold a read set """
	bc_ID = [None,8,8]
	def __init__(self,name,read_1,qual_1,read_2,qual_2,bc_1,bc_2):
		self.name = name
		self.read_1 = read_1
		self.qual_1 = qual_2
		self.read_2 = read_2
		self.qual_2 = qual_2
		self.bc_1 = bc_1
		self.bc_2 = bc_2
	def getBarcode(self, bcTable, max_diff):
		bc_pair = "%s%s" % (self.bc_1, self.bc_2)
		if bc_pair in bcTable.barcodes:
			self.bc_ID = bcTable.barcodes[bc_pair]
			return bcTable
		else:
		### Barcode One Matching ###
			bc1 = None
			bc1Mismatch = 0
			if self.bc_1 in bcTable.P7:
				bc1 = self.bc_1
			else:
				bc1Mismatch = len(self.bc_1)
				for key in bcTable.P7:
					bcdist = barcodeDist(key, self.bc_1)
					if bcdist < bc1Mismatch:
						bc1 = key
						bc1Mismatch = bcdist
			### Barcode Two Matching ###
			bc2 = None
			bc2Mismatch = 0
			if self.bc_2 in bcTable.P5:
				bc2 = self.bc_2
			else:
				bc2Mismatch = len(self.bc_2)
				for key in bcTable.P5:
					bcdist = barcodeDist(key, self.bc_2)
					if bcdist <= bc2Mismatch:
						bc2 = key
						bc2Mismatch = bcdist
			### Barcode Pair Matching ###
			combined_bc = [None,len(self.bc_1),len(self.bc_1)]
			if "%s%s" % (bc1, bc2) in bcTable.barcodes and bc1Mismatch <= max_diff and bc2Mismatch <= max_diff:
				combined_bc = [bcTable.barcodes["%s%s" % (bc1, bc2)][0],bc1Mismatch,bc2Mismatch]
			bcTable.barcodes[bc_pair] = combined_bc
			self.bc_ID = combined_bc
			return bcTable
	def writeRead(self):
		""" creat four line string for the read """
		read1_name = "%s 1:N:0:%s%s %s|%s|%s" % (self.name, self.bc_1, self.bc_2, self.bc_ID[0], self.bc_ID[1], self.bc_ID[2])
		read2_name = "%s 2:N:0:%s%s %s|%s|%s" % (self.name, self.bc_1, self.bc_2, self.bc_ID[0], self.bc_ID[1], self.bc_ID[2])
		r1 = '\n'.join([read1_name, self.read_1,'+',self.qual_1])
		r2 = '\n'.join([read2_name, self.read_2,'+',self.qual_2])
		return [r1,r2]

