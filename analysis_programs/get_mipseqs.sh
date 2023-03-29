#!/usr/bin/env bash

#Xander Nuttle
#get_mipseqs.sh
#Call: /data/talkowski/xander/MIPs/analysis_programs/get_mipseqs.sh gzipped_sam_file miptargets_file

REFERENCE_DIR=/var/tmp/xnuttle
CURRENT_DIR=$(pwd)
MFILE_NAME=$(basename ${2})
SAMP_NAME=$(basename ${1} .fastq.gz.sam.gz)

mkdir -p $REFERENCE_DIR
chgrp -R miket $REFERENCE_DIR
rsync -a --bwlimit=500 $CURRENT_DIR/$1 $REFERENCE_DIR
rsync -a --bwlimit=500 $2 $REFERENCE_DIR
cd $REFERENCE_DIR
/data/talkowski/xander/MIPs/analysis_programs/mip_seq_analysis $1 $MFILE_NAME
mv $REFERENCE_DIR/${SAMP_NAME}.mipseqs.gz $CURRENT_DIR
rm $REFERENCE_DIR/${SAMP_NAME}*

