#!/bin/sh
rm -rf tests/backtest
rm -rf tests/join
rm -rf tests/preprocess
rm -rf tests/splitreads

rm -rf tests/Amplicon_Raw_fastq/*Unidentified*
rm -rf tests/Amplicon_Raw_fastq/*Identified*
rm -rf tests/Amplicon_Raw_fastq/*Barcodes*

rm -rf dist
rm -rf build
rm -rf dbcAmplicons/*.pyc
rm -rf dbcAmplicons/*.so
rm -rf dbcAmplicons.egg-info
