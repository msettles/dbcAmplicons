#!/bin/sh

rm -rf tests/Amplicon_Raw_fastq/*Unidentified*
rm -rf tests/Amplicon_Raw_fastq/*Identified*
rm -rf tests/Amplicon_Raw_fastq/*Barcodes*
sudo rm -rf dist
sudo rm -rf build
sudo rm -rf dbcAmplicons/*.pyc
sudo rm -rf dbcAmplicons/*.so
sudo rm -rf dbcAmplicons.egg-info
