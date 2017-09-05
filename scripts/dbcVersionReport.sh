#!/usr/bin/env sh

echo "Analysis version report created."

touch version_report.txt
printf "Report generated on $(date)\n" > version_report.txt
dbcAmplicons --version >> version_report.txt 2>&1
flash2 --version >> version_report.txt 2>&1
java -jar $RDP_PATH version >> version_report.txt 2>&1
