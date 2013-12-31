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
try: 
	from dbcAmplicons import editdist
	editdist_loaded = True
except ImportError: 
	editdist_loaded = False 


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

#------------------- calculate distance between primer sequence and first part of read ------------------------------
def primerDist(primer,read, max_diff, end_match):
	'counts mismatches between primer and sequence'
	if editdist_loaded:
		return editdist.bounded_distance(primer,read,max_diff, end_match)
	else:
		read = read[0:len(primer)]
		if end_match == 0:
			dist = sum(map(lambda x: x[0] != x[1], zip(read, primer) ))
		else:
			dist = sum(map(lambda x: x[0] != x[1], zip(read[:-end_match], primer[:-end_match]) ))
			for i in range(len(primer)-end_match,len(primer)):
				if read[i] != primer[i]:
					dist += 100
					return [dist,len(primer)]
		return [dist,len(primer)]


# ---------------- Class for 4 read sequence data with double barcode set ----------------
class SequenceReadSet:
	""" class to hold a read set """
	goodRead = False
	bc_ID = None
	pr_ID = None
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
			if self.bc_ID[0] != None:
				self.goodRead = True
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
			combined_bc = [None,bc1Mismatch,bc2Mismatch]
			if "%s%s" % (bc1, bc2) in bcTable.barcodes.keys() and bc1Mismatch <= max_diff and bc2Mismatch <= max_diff:
				combined_bc[0] = bcTable.barcodes["%s%s" % (bc1, bc2)][0]
				self.goodRead = True
			bcTable.barcodes[bc_pair] = combined_bc
			self.bc_ID = combined_bc
			return bcTable
	def getPrimer(self, prTable, max_diff, endmatch):
		### Barcode One Matching ###
		pr1 = None
		pr1Mismatch = max_diff+1
		pr1Position = 0
		for key in prTable.P5sequences:
			prdist = primerDist(key, self.read_1,max_diff,endmatch)
			if prdist[0] < pr1Mismatch:
				pr1 = key
				pr1Mismatch = prdist[0]
				pr1Position = prdist[1]
		### Barcode Two Matching ###
		pr2 = None
		pr2Mismatch = max_diff+1
		pr2Position = 0
		for key in prTable.P7sequences:
			prdist = primerDist(key, self.read_2,max_diff,endmatch)
			if prdist[0] <= pr2Mismatch:
				pr2 = key
				pr2Mismatch = prdist[0]
				pr2Position = prdist[1]
		### Barcode Pair Matching ###
		combined_pr = prTable.getMatch(pr1,pr2)
		self.goodRead = self.goodRead and combined_pr[0] != None and pr1Mismatch <= max_diff and pr2Mismatch <= max_diff
		self.pr_ID = [combined_pr[0],combined_pr[1],pr2Mismatch,pr1Position,combined_pr[2],pr2Mismatch,pr2Position]
		return 1
	def writeRead(self):
		if (self.bc_ID == None):
			raise 
		""" create four line string for the read """
		if self.pr_ID != None:
			read1_name = "%s 1:N:0:%s:%s %s|%s|%s|%s %s|%s|%s" % (self.name, self.bc_ID[0], self.pr_ID[0], self.bc_1, self.bc_ID[1], self.bc_2 , self.bc_ID[2], self.pr_ID[1], self.pr_ID[2], self.pr_ID[3])
			read2_name = "%s 2:N:0:%s:%s %s|%s|%s|%s %s|%s|%s" % (self.name, self.bc_ID[0], self.pr_ID[0], self.bc_1, self.bc_ID[1], self.bc_2 , self.bc_ID[2], self.pr_ID[4], self.pr_ID[5], self.pr_ID[6])
			r1 = '\n'.join([read1_name, self.read_1[self.pr_ID[3]:],'+',self.qual_1[self.pr_ID[3]:]])
			r2 = '\n'.join([read2_name, self.read_2[self.pr_ID[6]:],'+',self.qual_2[self.pr_ID[6]:]])
		else:
			read1_name = "%s 1:N:0:%s %s|%s|%s|%s" % (self.name, self.bc_ID[0], self.bc_1, self.bc_ID[1], self.bc_2 , self.bc_ID[2])
			read2_name = "%s 2:N:0:%s %s|%s|%s|%s" % (self.name, self.bc_ID[0], self.bc_1, self.bc_ID[1], self.bc_2 , self.bc_ID[2])
			r1 = '\n'.join([read1_name, self.read_1,'+',self.qual_1])
			r2 = '\n'.join([read2_name, self.read_2,'+',self.qual_2])
		return [r1,r2]

