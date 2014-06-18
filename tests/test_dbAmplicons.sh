#!/usr/bin/env sh

#dbcAmplicons preprocess -b 15001 -B barcodeLookupTable.txt -P primerLookupTable.txt -1 Amplicon_Raw_fastq/Test100K_16S_R1_001.fastq.gz Amplicon_Raw_fastq/test40k_R1_001.fastq.gz -2 Amplicon_Raw_fastq/Test100K_16S_R2_001.fastq.gz Amplicon_Raw_fastq/test40k_R2_001.fastq.gz -3 Amplicon_Raw_fastq/Test100K_16S_R3_001.fastq.gz Amplicon_Raw_fastq/test40k_R3_001.fastq.gz -4 Amplicon_Raw_fastq/Test100K_16S_R4_001.fastq.gz Amplicon_Raw_fastq/test40k_R4_001.fastq.gz -O preprocess/trimL --debug
dbcAmplicons preprocess -b 15001 -B barcodeLookupTable.txt -P primerLookupTable.txt -1 Amplicon_Raw_fastq/Test100K_16S_R1_001.fastq.gz Amplicon_Raw_fastq/test40k_R1_001.fastq.gz  -O preprocess/trimL --debug
#cat $cmd
#system $cmd

# Current output should be

#No sample file identified
#barcode table length: 864
#primer table length P5 Primer Sequences:7, P7 Primer Sequences:7
#processed 15001 total reads, 4345.0 Reads/second, 2103 identified reads, 12898 unidentified reads
#processed 30002 total reads, 4404.0 Reads/second, 4299 identified reads, 25703 unidentified reads
#processed 45003 total reads, 4051.0 Reads/second, 9651 identified reads, 35352 unidentified reads
#processed 60004 total reads, 3367.0 Reads/second, 21796 identified reads, 38208 unidentified reads
#processed 75005 total reads, 3032.0 Reads/second, 34066 identified reads, 40939 unidentified reads
#processed 90006 total reads, 2816.0 Reads/second, 46271 identified reads, 43735 unidentified reads
#processed 105007 total reads, 2643.0 Reads/second, 58547 identified reads, 46460 unidentified reads
#processed 120008 total reads, 2513.0 Reads/second, 70845 identified reads, 49163 unidentified reads
#processed 135009 total reads, 2438.0 Reads/second, 83193 identified reads, 51816 unidentified reads
#processed 140000 total reads, 2421.0 Reads/second, 87346 identified reads, 52654 unidentified reads
#140000 reads processed in 0.96 minutes
#Cleaning up.

#dbcAmplicons splitreads -b 15001 -S sampleLookupTable.txt -1 preprocess/trimL_R1.fastq.gz -2 preprocess/trimL_R2.fastq.gz -O splitreads
dbcAmplicons splitreads -b 15001 -S sampleLookupTable.txt -1 preprocess/trimL_R1.fastq.gz -O splitreads

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

#dbcAmplicons join -t 4 -x 0.25 -1 splitreads/match_twoprimer_R1.fastq.gz -2 splitreads/match_twoprimer_R2.fastq.gz -O join/match_twoprimer
dbcAmplicons join -t 4 -x 0.25 -1 splitreads/match_twoprimer_R1.fastq.gz -O join/match_twoprimer

# Results should be

#[FLASH] Starting FLASH v1.2.8
#[FLASH] Fast Length Adjustment of SHort reads
#[FLASH]
#[FLASH] Input files:
#[FLASH]     splitreads/match_twoprimer_R1.fastq.gz
#[FLASH]     splitreads/match_twoprimer_R2.fastq.gz
#[FLASH]
#[FLASH] Output files:
#[FLASH]     ./join/match_twoprimer.extendedFrags.fastq.gz
#[FLASH]     ./join/match_twoprimer.notCombined_1.fastq.gz
#[FLASH]     ./join/match_twoprimer.notCombined_2.fastq.gz
#[FLASH]     ./join/match_twoprimer.hist
#[FLASH]     ./join/match_twoprimer.histogram
#[FLASH]
#[FLASH] Parameters:
#[FLASH]     Min overlap:          10
#[FLASH]     Max overlap:          65
#[FLASH]     Phred offset:         33
#[FLASH]     Combiner threads:     4
#[FLASH]     Max mismatch density: 0.250000
#[FLASH]     Cap mismatch quals:   false
#[FLASH]     Output format:        gzip
#[FLASH]     Interleaved input:    false
#[FLASH]     Interleaved output:   false
#[FLASH]
#[FLASH] Starting FASTQ readers and writer threads
#[FLASH] Starting 4 combiner threads
#[FLASH] Processed 25000 read pairs
#[FLASH] Processed 31438 read pairs
#[FLASH]
#[FLASH] Read combination statistics:
#[FLASH]     Total reads:      31438
#[FLASH]     Combined reads:   25357
#[FLASH]     Uncombined reads: 6081
#[FLASH]     Percent combined: 80.66%
#[FLASH]
#[FLASH] Writing histogram files.
#[FLASH]
#[FLASH] FLASH v1.2.8 complete!
#[FLASH] 3.178 seconds elapsed

dbcAmplicons classify -b 7500 -O join/classify -U join/match_twoprimer.extendedFrags.fastq.gz --debug  --rdpPath=/Users/msettles/opt/RDPTools/classifier.jar -p 4 -1 join/match_twoprimer.notCombined_1.fastq.gz -2 join/match_twoprimer.notCombined_2.fastq.gz

# Result should be

#Starting rdp for file joint/classify.7500.fasta
#Starting rdp for file joint/classify.15000.fasta
#Starting rdp for file joint/classify.22500.fasta
#Starting rdp for file joint/classify.25357.fasta
#Finished processing joint/classify.25357.fasta in 1.36 minutes
#Starting rdp for file joint/classify.31438.fasta
#Finished processing joint/classify.22500.fasta in 2.97 minutes
#Finished processing joint/classify.15000.fasta in 3.0 minutes
#Finished processing joint/classify.7500.fasta in 3.01 minutes
#Finished processing joint/classify.31438.fasta in 2.04 minutes
#Combining temporary files
#31438 reads processed in 3.43 minutes
#Cleaning up.


#convert2ReadTo4Read.py -O test -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz -2 Amplicon_Raw_fastq/Hkates_R2.40k.fastq.gz
convert2ReadTo4Read.py -O backtest/test -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz 


#dbcAmplicons preprocess -B barcodeLookupTable.txt -O test2 -1 test_R1_001.fastq.gz -2 test_R2_001.fastq.gz -3 test_R3_001.fastq.gz -4 test_R4_001.fastq.gz
dbcAmplicons preprocess -B barcodeLookupTable.txt -O backtest/test2 -1 backtest/test_R1_001.fastq.gz

## 100% should be identified
