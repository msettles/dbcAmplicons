#!/usr/bin/env sh

#dbcAmplicons preprocess -b 15001 -B barcodeLookupTable.txt -P primerLookupTable.txt -1 Amplicon_Raw_fastq/Test100K_16S_R1_001.fastq.gz Amplicon_Raw_fastq/test40k_R1_001.fastq.gz -2 Amplicon_Raw_fastq/Test100K_16S_R2_001.fastq.gz Amplicon_Raw_fastq/test40k_R2_001.fastq.gz -3 Amplicon_Raw_fastq/Test100K_16S_R3_001.fastq.gz Amplicon_Raw_fastq/test40k_R3_001.fastq.gz -4 Amplicon_Raw_fastq/Test100K_16S_R4_001.fastq.gz Amplicon_Raw_fastq/test40k_R4_001.fastq.gz -O preprocess/trimL --debug
#dbcAmplicons preprocess -b 15001 -B barcodeLookupTable.txt -P primerLookupTable.txt -1 Amplicon_Raw_fastq/Test100K_16S_R1_001.fastq.gz Amplicon_Raw_fastq/test40k_R1_001.fastq.gz  -O preprocess/trimL --debug
dbcAmplicons preprocess -b 15001 -S sampleLookupTable.txt -B barcodeLookupTable.txt -P primerLookupTable.txt -1 Amplicon_Raw_fastq/Test100K_16S_R1_001.fastq.gz Amplicon_Raw_fastq/test40k_R1_001.fastq.gz  -U -O preprocess --debug

#cat $cmd
#system $cmd

#barcode table length: 864
#primer table length P5 Primer Sequences:7, P7 Primer Sequences:7
#sample table length: 208, and 5 projects.
#processed 15001 total reads, 2704.0 Reads/second, 568 identified reads(3.8%), 14433 unidentified reads
#processed 30002 total reads, 2635.0 Reads/second, 1138 identified reads(3.8%), 28864 unidentified reads
#processed 45003 total reads, 2704.0 Reads/second, 5218 identified reads(11.6%), 39785 unidentified reads
#processed 60004 total reads, 2908.0 Reads/second, 16506 identified reads(27.5%), 43498 unidentified reads
#processed 75005 total reads, 3044.0 Reads/second, 27832 identified reads(37.1%), 47173 unidentified reads
#processed 90006 total reads, 3142.0 Reads/second, 39140 identified reads(43.5%), 50866 unidentified reads
#processed 105007 total reads, 3217.0 Reads/second, 50455 identified reads(48.0%), 54552 unidentified reads
#processed 120008 total reads, 3274.0 Reads/second, 61793 identified reads(51.5%), 58215 unidentified reads
#processed 135009 total reads, 3228.0 Reads/second, 73156 identified reads(54.2%), 61853 unidentified reads
#processed 140000 total reads, 3203.0 Reads/second, 76970 identified reads(55.0%), 63030 unidentified reads
#140000 reads processed in 0.73 minutes, 76970 (55.0%) identified

#63030 reads (45.0% of total run) found for project	Unidentified
#18244 reads (13.0% of total run) found for project	match_wildcard
#25093 reads (17.9% of total run) found for project	match_twoprimersecond
#31438 reads (22.5% of total run) found for project	match_twoprimer
#2195 reads (1.6% of total run) found for project	subfolder/match_16S
#0 reads (0.0% of total run) found for project	nomatch
#Cleaning up.



#dbcAmplicons splitreads -b 15001 -S sampleLookupTable.txt -1 preprocess/trimL_R1.fastq.gz -2 preprocess/trimL_R2.fastq.gz -O splitreads
#dbcAmplicons splitreads -b 15001 -S sampleLookupTable.txt -1 preprocess/trimL_R1.fastq.gz -2 preprocess/trimL_R2.fastq.gz -O splitreads --debug

# Current output should be

#sample table length: 208, and 5 projects.
#processed 15001 total reads, 3766.0 Reads/second, 10188 identified reads, 4813 unidentified reads
#processed 30002 total reads, 3289.0 Reads/second, 24073 identified reads, 5929 unidentified reads
#processed 45003 total reads, 3167.0 Reads/second, 37964 identified reads, 7039 unidentified reads
#processed 60004 total reads, 3107.0 Reads/second, 51809 identified reads, 8195 unidentified reads
#processed 75005 total reads, 3076.0 Reads/second, 65595 identified reads, 9410 unidentified reads
#processed 87346 total reads, 3030.0 Reads/second, 76970 identified reads, 10376 unidentified reads
#87346 reads processed in 0.48 minutes
#25093	reads found for project	match_twoprimersecond
#31438	reads found for project	match_twoprimer
#2195	reads found for project	subfolder/match_16S
#18244	reads found for project	match_wildcard
#0	reads found for project	nomatch
#Cleaning up.

dbcAmplicons join -t 4 -x 0.25 -1 preprocess/match_twoprimer_R1.fastq.gz -2 preprocess/match_twoprimer_R2.fastq.gz -O join/match_twoprimer -v
#dbcAmplicons join -t 4 -x 0.25 -1 splitreads/match_twoprimer_R1.fastq.gz -O join/match_twoprimer
#dbcAmplicons join -t 4 -x 0.25 -1 preprocess/trimL_R1.fastq.gz -O join/all_primers

# Results should be

#Using Flash_version:v1.2.11
#Min_overlap:10
#Max_overlap:600
#Max_mismatch_density:0.250000
#Allow_"outie"_pairs:true
#Cap_mismatch_quals:false
#Combiner_threads:4
#Input_format:FASTQ, phred_offset=33
#Output_format:FASTQ, phred_offset=33, gzip
#Total_pairs:31438
#Combined_pairs:28132
#Innie_pairs:27910 (99.21% of combined)
#Outie_pairs:222 (0.79% of combined)
#Uncombined_pairs:3306
#Percent_combined:89.48%



dbcAmplicons classify -b 7500 -O join/classify -U join/match_twoprimer.extendedFrags.fastq.gz --debug  --rdpPath /home/alida/Documents/Idaho/bioinformatics/RDPTools/classifier/dist/classifier.jar -p 4 -1 join/match_twoprimer.notCombined_1.fastq.gz -2 join/match_twoprimer.notCombined_2.fastq.gz

# Result should be

#Starting rdp for file join/classify.7500.fasta
#Starting rdp for file join/classify.15000.fasta
#Starting rdp for file join/classify.22500.fasta
#Starting rdp for file join/classify.28132.fasta
#Finished processing join/classify.28132.fasta in 1.21 minutes
#Starting rdp for file join/classify.31438.fasta
#Finished processing join/classify.15000.fasta in 1.35 minutes
#Finished processing join/classify.22500.fasta in 1.52 minutes
#Finished processing join/classify.7500.fasta in 1.54 minutes
#Finished processing join/classify.31438.fasta in 0.44 minutes
#Combining temporary files
#31438 reads processed in 1.66 minutes
#Cleaning up.


dbcAmplicons abundance -O join/abundance -F join/classify.fixrank --debug

# Output should be 
#31438 lines processed in 0.01 minutes
#Writing output
#finished in 0.01 minutes
#Cleaning up.

#convert2ReadTo4Read.py -O test -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz -2 Amplicon_Raw_fastq/Hkates_R2.40k.fastq.gz

convert2ReadTo4Read.py -O backtest/test -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz 

#processed 40000 total reads, 3844.0 Reads/second
#40000 reads processed in 0.17 minutes
#Cleaning up.

SplitReadsBySample.py -O allSamples -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz -2 Amplicon_Raw_fastq/Hkates_R2.40k.fastq.gz --debug
#SplitReadsBySample.py -O allSamples -U Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz --debug

#dbcAmplicons preprocess -B barcodeLookupTable.txt -O test2 -1 test_R1_001.fastq.gz -2 test_R2_001.fastq.gz -3 test_R3_001.fastq.gz -4 test_R4_001.fastq.gz

dbcAmplicons preprocess -B barcodeLookupTable.txt -O backtest/test2 -1 backtest/test_R1_001.fastq.gz
#No sample file identified
#No primer file identified
#barcode table length: 864
#processed 40000 total reads, 3719.0 Reads/second, 39186 identified reads, 814 unidentified reads (98.0%)
#40000 reads processed in 0.18 minutes, 39186 (98.0%) identified
#Cleaning up.
