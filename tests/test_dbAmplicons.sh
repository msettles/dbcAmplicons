#!/usr/bin/env sh

echo '' > test_output.txt
## test validation
echo "Testing dbcAmplicons validate, this should fail"
dbcAmplicons validate -B barcodeLookupTable.txt -P primerLookupTable-bd.txt -S sampleLookupTable-bd.txt >> test_output.txt

## testing preprocess
echo "Testing dbcAmplicons preprocess"
dbcAmplicons preprocess -b 15001 -q 10 -l 200 -S sampleLookupTable.txt -B barcodeLookupTable.txt -P primerLookupTable-dedup.txt -f 6 -1 Amplicon_Raw_fastq/Test100K_16S_R1_001.fastq.gz Amplicon_Raw_fastq/test40k_R1_001.fastq.gz -O preprocess/trimL --debug  >> test_output.txt

#echo "Testing dbcAmplicons splitreads"
#dbcAmplicons splitreads -b 15001 -S sampleLookupTable.txt -1 preprocess/trimL_R1.fastq.gz -O splitreads --debug

echo "Testing dbcAmplicons join"
dbcAmplicons join -t 4 -x 0.25 -1 preprocess/trimL/match_twoprimer_R1.fastq.gz -O join/match_twoprimer  >> test_output.txt

echo "Testing dbcAmplicons classify"
dbcAmplicons classify -b 7500 -q 10 -l 200 -O join/classify -U join/match_twoprimer.extendedFrags.fastq.gz --debug -p 4 -1 join/match_twoprimer.notCombined_1.fastq.gz -2 join/match_twoprimer.notCombined_2.fastq.gz  >> test_output.txt

#echo "Testing dbcAmplicons abundance"
#dbcAmplicons abundance -O join/abundance -F join/classify.fixrank --debug >> test_output.txt 

echo "Testing dbcAmplicons abundance (biom format)"
dbcAmplicons abundance -O join/abundance -F join/classify.fixrank -S sampleLookupTable.txt -b --debug >> test_output.txt

echo "Testing dbcAmplicons abundance (biom format, hdf5)"
dbcAmplicons abundance -O join/abundance -F join/classify.fixrank -S sampleLookupTable.txt --hdf5 --debug >> test_output.txt

echo "Testing convert2ReadTo4Read"
convert2ReadTo4Read.py -O backtest/test -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz >> test_output.txt

echo "Testing dbcAmplicons preprocess post convert2ReadTo4Read back test"
dbcAmplicons preprocess -B barcodeLookupTable.txt -O backtest/test2 -1 backtest/test_R1.fastq.gz >> test_output.txt

echo "Testing splitReadsBySample paired reads"
splitReadsBySample.py -O allSamples -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz -2 Amplicon_Raw_fastq/Hkates_R2.40k.fastq.gz --debug >> test_output.txt

echo "Testing splitReadsBySample single reads"
splitReadsBySample.py -O allSamples -U Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz --debug >> test_output.txt

echo "Testing screening"
dbcAmplicons screen -R Amplicon_Raw_fastq/test_map.fa -1 Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz -2 Amplicon_Raw_fastq/Hkates_R2.40k.fastq.gz -U Amplicon_Raw_fastq/Hkates_R1.40k.fastq.gz --debug -s 2  -o >> test_output.txt

#echo "Removing folders and files"
#rm -rf allSamples backtest join preprocess splitreads out.txt

diff test_validate.txt test_output.txt
