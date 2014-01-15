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


import re
import sys
class KeyFoundError(Exception):
	def __init__(self, barcode, primer):
		self.barcode = barcode
		self.primer = primer
	def getBarcode(self):
		return repr(self.barcode)
	def getPrimer(self):
		return repr(self.primer)

# ---------------- samples class ----------------
class sampleTable:
	sampleCount = 0
	""" class to hold sample information """
	def __init__(self, samplefile):
		self.sampleTable = {}
		try:
			sfile = open(samplefile, 'r')
		except IOError:
			print 'cannot open', samplefile
		f = sfile.readlines() ## read the file
		header = f[0].rstrip()
		vheader = header.split('\t')
		try:
			sampleID_index = vheader.index("SampleID")
		except ValueError:
			print "Samples: Column 'SampleID' was not found in the samples table"
			raise
		try:
			barcodeID_index = vheader.index("BarcodeID")
		except ValueError:
			print "Samples: Column 'BarcodeID' was not found in the samples table"
			raise
		try:
			primerID_index = vheader.index("PrimerPairID")
		except ValueError:
			print "Samples: Column 'PrimerPairID' was not found in the samples table"
			raise
		try:
			projectID_index = vheader.index("ProjectID")
		except ValueError:
			print "Samples: Column 'ProjectID' was not found in the samples table"
			raise
		try:
			line=1
			for row in f:
				line += 1
				if row[0] == "#": # comment line
					continue
				row = row.rstrip() # strip off newline
				row = row.split('\t') # split by tab
				barcode = row[barcodeID_index]
				sid = re.sub(r'[\\/:"\' ]+', ".",row[sampleID_index]) # replace unsafe characters from sampleID with '.'
				pid = re.sub(r'[\\/:"\'*?<>| ]+', ".",row[projectID_index]) # replace unsafe characters from projectID with '.'
				for primer in row[primerID_index].split(','):
					primer = primer.strip()
					if primer == '*' and barcode in self.sampleTable.keys():
						raise KeyFoundError(barcode,primer)
					elif barcode in self.sampleTable.keys():
						if primer in self.sampleTable[barcode].keys():
							raise KeyFoundError(barcode,primer)
						else:
							self.sampleTable[barcode][primer] = [sid,pid]
					else:
						self.sampleTable[barcode] = {}
						self.sampleTable[barcode][primer] = [sid,pid]
				self.sampleCount += 1
		except KeyFoundError as e:
			print "Samples: Error in sample sheet: Barcode [%s] and Primer [%s] pair previously specified" % (e.getBarcode(), e.getPrimer())
			sfile.close()
			raise
		except:
			print "Samples: Unexpected error on line %s: %s" % (line,sys.exc_info()[0])
			sfile.close()
			raise
		sfile.close()
	def getSampleNumber(self):
		return repr(self.sampleCount)
	def getSampleID(self,barcode,primer):
		try:
			sid = self.sampleTable[barcode]
			if sid.keys() == "*":
				return sid["*"][0]
			else:
				return(sid["*"][0])
		except KeyError:
			return(None)
		except:
			print "Samples: Unexpected error:", sys.exc_info()[0]
			raise
	def getProjectID(self,barcode,primer):
		try:
			sid = self.sampleTable[barcode]
			if sid.keys() == "*":
				return sid["*"][1]
			else:
				return(sid["*"][1])
		except KeyError:
			return(None)
		except:
			print "Samples: Unexpected error:", sys.exc_info()[0]
			raise
