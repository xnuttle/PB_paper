#!/usr/bin/env bash

#Xander Nuttle
#prep_fastqs.sh
#Call: /data/talkowski/xander/MIPs/analysis_programs/prep_fastqs.sh sampleset_name merging_input_directory molecular_tag_length

REFERENCE_DIR=/var/tmp/xnuttle
CURRENT_DIR=$(pwd)
SAMPLE=$(cut -f1 ${1}.barcodekey)

mkdir -p $REFERENCE_DIR
chgrp -R miket $REFERENCE_DIR
chmod -R g+wx $REFERENCE_DIR
rsync -a --bwlimit=500 $CURRENT_DIR/$1* $REFERENCE_DIR
cd $REFERENCE_DIR
/data/talkowski/xander/MIPs/analysis_programs/dm_fastq_to_fastq_for_pear ${1}.r1.fastq.gz ${1}.bc1.fastq.gz ${1}.bc2.fastq.gz ${1}.r2.fastq.gz 142 $3 250000 ${1}.barcodekey
mv $REFERENCE_DIR/${SAMPLE}_FS*fastq.gz $2
rm $REFERENCE_DIR/$1*

