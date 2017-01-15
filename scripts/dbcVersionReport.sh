#!/usr/bin/env sh

echo "Analysis version report created."

touch version_report.txt
printf "Report generated on $(date)\n" > version_report.txt
dbcAmplicons --version >> version_report.txt 2>&1
flash2 --version >> version_report.txt 2>&1
printf "RDP database version 11.4 (May 26, 2015)\n" >> version_report.txt
printf "RDP classifier version 2.10.2\n" >> version_report.txt