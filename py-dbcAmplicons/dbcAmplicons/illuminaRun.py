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

import os
import glob
import gzip
from dbcAmplicons import SequenceReadSet

class IlluminaRun:
	isOpen = False
	count = 0
	""" A sequencing run """
	def __init__(self,file_prefix):
		self.prefix = file_prefix
		self.fread1 = glob.glob("%s*R1_[0-9][0-9][0-9].fastq.gz" % file_prefix)
		if len(self.fread1) == 0 or not all(os.path.exists(f) for f in self.fread1):
			raise Exception("R1 file not found")
		self.fread2 = glob.glob("%s*R2_[0-9][0-9][0-9].fastq.gz" % file_prefix)
		if len(self.fread2) == 0 or not all(os.path.exists(f) for f in self.fread2):
			raise Exception("R2 file not found")
		self.fread3 = glob.glob("%s*R3_[0-9][0-9][0-9].fastq.gz" % file_prefix)
		if len(self.fread3) == 0 or not all(os.path.exists(f) for f in self.fread3):
			raise Exception("R3 file not found")
		self.fread4 = glob.glob("%s*R4_[0-9][0-9][0-9].fastq.gz" % file_prefix)
		if len(self.fread4) == 0 or not all(os.path.exists(f) for f in self.fread4):
			raise Exception("R4 file not found")
		if any(len(f) != len(self.fread1) for f in [self.fread2,self.fread3,self.fread4]):
			raise Exception("Inconsistant number of files for each read")
		# record the number of files per read
		self.numberoffiles = len(self.fread1)
		if self.numberoffiles != 1 :
			raise Exception("Number of files per read exceeds 1")
	def open(self):
		if self.fread1[0].split(".")[-1] == "gz":
			self.R1 = gzip.open(self.fread1[0], 'rb')
		else:
			self.R1 = open(self.fread1[0], 'r')
		if self.fread2[0].split(".")[-1] == "gz":
			self.R2 = gzip.open(self.fread2[0], 'rb')
		else:
			self.R2 = open(self.fread2[0], 'r')
		if self.fread3[0].split(".")[-1] == "gz":
			self.R3 = gzip.open(self.fread3[0], 'rb')
		else:
			self.R3 = open(self.fread3[0], 'r')
		if self.fread4[0].split(".")[-1] == "gz":
			self.R4 = gzip.open(self.fread4[0], 'rb')
		else:
			self.R4 = open(self.fread4[0], 'r')
		self.isOpen = True
		self.count=0
	def close(self):
		self.R1.close()
		self.R2.close()
		self.R3.close()
		self.R4.close()
		self.isOpen = False
	def next(self, count=1):
		if not self.isOpen:
			self.open()
		reads = []
		for i in range(0,count):
			try:
				#process read_1 (Read1)
				name = self.R1.next().split(" ")[0] # name
				read_1 = self.R1.next().rstrip()	# read
				self.R1.next()	# '+'
				qual_1 = self.R1.next().rstrip() # qual
				#process read _2 (BC1)
				self.R2.next() # name
				bc_1 = self.R2.next().rstrip() # read 
				self.R2.next() # '+'
				self.R2.next() # qual
				#process reed_3 (BC2)
				self.R3.next() # name
				bc_2 = self.R3.next().rstrip() # read
				self.R3.next() # '+'
				self.R3.next() # qual
				#process read_4 (Read2)
				self.R4.next() # name
				read_2 = self.R4.next().rstrip() # read
				self.R4.next() # '+'
				qual_2 = self.R4.next().rstrip() # qual
				#add it to the stack
			except StopIteration:
				break
			reads.append(SequenceReadSet(name=name,read_1=read_1,qual_1=qual_1,read_2=read_2,qual_2=qual_2,bc_1=bc_1,bc_2=bc_2))
			self.count += 1
		return reads


