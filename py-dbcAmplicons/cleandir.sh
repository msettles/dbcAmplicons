#!/bin/sh

rm -rf tests/Amplicon_Raw_fastq/*Unidentified*
rm -rf tests/Amplicon_Raw_fastq/*Identified*
rm -rf tests/Amplicon_Raw_fastq/Test100K_16S_Barcodes.txt
sudo rm -rf dist
sudo rm -rf build
sudo rm -rf dbcAmplicons/*.pyc
sudo rm -rf dbcAmplicons.egg-info
