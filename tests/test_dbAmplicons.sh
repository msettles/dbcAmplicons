#!/usr/bin/env sh

RDP_PATH="~/opt/RDPTools/classifier/dist/classifier.jar"

## testing preprocess
echo "Testing dbcAmplicons preprocess"
dbcAmplicons preprocess -b 15001 -B barcodeLookupTable.txt -P primerLookupTable.txt -1 Amplicon_Raw_fastq/Test100K_16S_R1_001.fastq.gz Amplicon_Raw_fastq/test40k_R1_001.fastq.gz  -O preprocess/trimL --debug &> out.txt

echo "Testing dbcAmplicons splitreads"
dbcAmplicons splitreads -b 15001 -S sampleLookupTable.txt -1 preprocess/trimL_R1.fastq.gz -O splitreads --debug >> out.txt 2>&1

echo "Testing dbcAmplicons join"
dbcAmplicons join -t 4 -x 0.25 -1 splitreads/match_twoprimer_R1.fastq.gz -O join/match_twoprimer >> out.txt 2>&1

echo "Testing dbcAmplicons classify"
dbcAmplicons classify -b 7500 -O join/classify -U join/match_twoprimer.extendedFrags.fastq.gz --debug  --rdpPath $RDP_PATH -p 4 -1 join/match_twoprimer.notCombined_1.fastq.gz -2 join/match_twoprimer.notCombined_2.fastq.gz >> out.txt 2>&1

echo "Testing dbcAmplicons abundance"
dbcAmplicons abundance -O join/abundance -F join/classify.fixrank --debug >> out.txt 2>&1

echo "Testing convert2ReadTo4Read"
convert2ReadTo4Read.py -O backtest/test -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz >> out.txt 2>&1

echo "Testing dbcAmplicons preprocess post convert2ReadTo4Read back test"
dbcAmplicons preprocess -B barcodeLookupTable.txt -O backtest/test2 -1 backtest/test_R1_001.fastq.gz >> out.txt 2>&1

echo "Testing splitReadsBySample paired reads"
splitReadsBySample.py -O allSamples -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz -2 Amplicon_Raw_fastq/Hkates_R2.40k.fastq.gz --debug >> out.txt 2>&1

echo "Testing splitReadsBySample single reads"
splitReadsBySample.py -O allSamples -U Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz --debug >> out.txt 2>&1

echo "DIFFERENCES BETWEEN TESTS AND EXPECTED"
echo "######################################"
diff out.txt expected.txt
echo "######################################"

echo "Removing folders and files"
rm -rf allSamples backtest join preprocess splitreads out.txt

