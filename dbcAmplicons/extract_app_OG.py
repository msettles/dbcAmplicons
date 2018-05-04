#extract_app

import sys, os

class extractApp():

	def __init__(self):
		self.verbose = False

	def start(self, taxon, fixrank, fastq, threshold, output):
		#Format inputs:
		print("taxon target is: " + taxon)
		fixrank_file = open(fixrank, "r")
		vsequences = open(fastq, "r")

		#generate list of fasta headers to be extracted:
		headers = []
		sampleIDs = []
		with fixrank_file as f:
			for line in f.readlines():
		 		if taxon in line:
		 			if line.split("\t",26)[25] >= float(threshold):
		 				headers.append(line.split("|")[0])
		 				sampleIDs.append(line.split('|')[1].split(':')[0])

		#Test for file presence
		if os.path.isfile(output+"."+taxon+".fasta"):
			open(output+"."+taxon+".fasta", "w").close()
		vOutput = open(output+"."+taxon+".fasta", "a")

		#Extract selected reads from fasta
		with vsequences as file:
			vReads = file.read().replace('\n','')
		vReadList = vReads.split('@M')

		n = 0
		with vOutput as myfile:
			for header in headers:
				for read in vReadList:
					if ('M'+read[:len(headers[n])-1]) == headers[n]:
						myfile.write('>' + sampleIDs[n] + '~' + read.split()[0] + '\n' + read.split()[3].split('|')[3].split('+')[0] + '\n')
				n +=1