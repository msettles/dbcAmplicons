#!/usr/bin/env sh

printf "Report generated on $(date)\n"
dbcAmplicons --version
flash2 --version
java -jar $RDP_PATH/classifier.jar version
