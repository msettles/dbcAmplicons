#extract_app

import sys, os, time

class extractApp():

	def __init__(self):
		self.verbose = False

	def start(self, taxon, fixrank, fastq, threshold, output):
		#Format inputs:
		print("taxon target is: " + taxon)
		fixrank_file = open(fixrank, "r")
		vsequences = open(fastq, "r")

		#ID all fastq reads
		starttime = time.time()
		vReadList = []
		with vsequences as file:
			for read in file.read().split('@M'):
				if read!='':
					vReadList.append('M'+read)
		sys.stderr.write("Finished fastq ID in "+ str(round((time.time()-starttime)/60, 10))+ " minutes\n")

		#Test for output file presence
		if os.path.isfile(output+"."+taxon+".fasta"):
			open(output+"."+taxon+".fasta", "w").close()
		vOutput = open(output+"."+taxon+".fasta", "a")

		#Filter classification results based on taxa and quality score
		starttime = time.time()
		vClassList = []
		with fixrank_file as file:
			for line in file.read().split('\n'):
				if line.split()[22] == taxon:
					if float(line.split("\t",26)[25]) >= float(threshold):
						vClassList.append(line)
		sys.stderr.write("Finished filtering class file in "+ str(round((time.time()-starttime)/60, 10))+ " minutes\n")
						
		#Extract sequence data for classified reads
		starttime = time.time()
		with vOutput as myfile:
			for read in vReadList:
				for classification in vClassList:
					if classification.split('|')[0] == read.split()[0]:
						myfile.write('>' + classification.split('|')[1].split(':')[0] + '~' + read.split('\n')[0].replace(classification.split('|')[1].split(':')[0],'') + '\n' + read.split('\n')[1] + '\n')
						vReadList.remove(read)
						break
		sys.stderr.write("Finished extracting reads in "+ str(round((time.time()-starttime)/60, 10))+ " minutes\n")