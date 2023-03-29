#!/usr/bin/env bash

#Xander Nuttle
#run_pear.sh
#Call: /data/talkowski/xander/MIPs/analysis_programs/run_pear.sh for_fastq_file rev_fastq_file (long)fileset_number merging_output_directory

REFERENCE_DIR=/var/tmp/xnuttle
CURRENT_DIR=$(pwd)
OUTNAME=$(echo $2|sed 's/FS1_R/FS1_M/g')

mkdir -p $REFERENCE_DIR
chgrp -R miket $REFERENCE_DIR
chmod -R g+wx $REFERENCE_DIR
rsync -a --bwlimit=500 $CURRENT_DIR/$1 $REFERENCE_DIR
rsync -a --bwlimit=500 $CURRENT_DIR/$2 $REFERENCE_DIR
cd $REFERENCE_DIR
/PHShome/adn5/programs/PEAR/pear-0.9.11-linux-x86_64/bin/pear -f $REFERENCE_DIR/$1 -r $REFERENCE_DIR/$2 -o $REFERENCE_DIR/merged_${3}
gzip $REFERENCE_DIR/merged_${3}.assembled.fastq
mv $REFERENCE_DIR/merged_${3}.assembled.fastq.gz $4/$OUTNAME
rm $REFERENCE_DIR/merged_${3}.discarded.fastq $REFERENCE_DIR/merged_${3}.unassembled.forward.fastq $REFERENCE_DIR/merged_${3}.unassembled.reverse.fastq
rm $REFERENCE_DIR/$1 $REFERENCE_DIR/$2

