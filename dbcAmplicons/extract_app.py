#extract_app

import sys, os, time
from dbcAmplicons import OneReadIlluminaRun

#starttime = time.time()
#sys.stderr.write("Indexed fastq reads in "+ str(round((time.time()-starttime)/60, 4))+ " minutes\n")

class extractApp():

	def __init__(self):
		self.verbose = False

	def start(self, taxon, fixrank, fastq, threshold, output, batchsize):

		#Format inputs:
		print("---\nTaxon target is: " + taxon)
		vFixrankFile 	= open(fixrank, "r")
		vFastqFile 		= open(fastq, "r")

		#Extract/filter fixrank classifications
		starttime = time.time()
		vClassList = []
		with vFixrankFile as file:
			for line in file.read().split('\n'):
				try:
					if line.split()[22] == taxon:
						if float(line.split("\t",26)[25]) >= float(threshold):
							vClassList.append(line.split('|')[0])
				except:
					continue
		sys.stderr.write("---\nFiltered and extracted fixrank classifications in  "+ str(round((time.time()-starttime)/60, 4))+ " minutes\n---")

		#Write output
		starttime = time.time()
		if os.path.isfile(output+"."+taxon+".fasta"):
			open(output+"."+taxon+".fasta", "w").close()
		vOutFile = (output+"."+taxon+".fasta", "a", 1)
		vOut = open(output+"."+taxon+".fasta", "a", 1)
		
		#Match and extract fastq reads
		starttime = time.time()
		vOutput = []
		vBatchNumber = 0
		with vFastqFile as file:
			for read in file.read().split('@M'):
				if vBatchNumber == batchsize:
					print(str(len(vClassList))+' reads left to extract.')
					for i in vOutput:
  						vOut.write(i)
					vOutput = []
					vBatchNumber = 0
				elif len(vClassList) == 0:
					for i in vOutput:
  						vOut.write(i)
					break
				for classification in vClassList:
					if ('M'+read.split('\n')[0].split(' ')[0]) == classification:
						vOutput.append('>M'+read.split('\n')[0]+read.split('\n')[1] + '\n')
						vClassList.remove(classification)
						break
				vBatchNumber +=1
		sys.stderr.write("---\nExtracteds fastq reads in  "+ str(round((time.time()-starttime)/60, 4))+ " minutes\n")