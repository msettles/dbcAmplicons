#!/usr/bin/env sh

dbcAmplicons preprocess -B 15001 -b barcodeLookupTable.txt -p primerLookupTable.txt -1 Amplicon_Raw_fastq/Test100K_16S_R1_001.fastq.gz Amplicon_Raw_fastq/test40k_R1_001.fastq.gz -2 Amplicon_Raw_fastq/Test100K_16S_R2_001.fastq.gz Amplicon_Raw_fastq/test40k_R2_001.fastq.gz -3 Amplicon_Raw_fastq/Test100K_16S_R3_001.fastq.gz Amplicon_Raw_fastq/test40k_R3_001.fastq.gz -4 Amplicon_Raw_fastq/Test100K_16S_R4_001.fastq.gz Amplicon_Raw_fastq/test40k_R4_001.fastq.gz -o preprocess/trimL --debug


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
