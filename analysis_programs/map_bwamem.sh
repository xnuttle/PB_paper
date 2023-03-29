#!/usr/bin/env bash

#Xander Nuttle
#map_bwamem.sh
#Call: /data/talkowski/xander/MIPs/analysis_programs/map_bwamem.sh fastq_file genome_fasta_file genome_size_file mapping_output_directory

REFERENCE_DIR=/var/tmp/xnuttle
CURRENT_DIR=$(pwd)
FASTA_NAME=$(basename ${2})
INDEX_NAME=$(basename ${3})

mkdir -p $REFERENCE_DIR
chgrp -R miket $REFERENCE_DIR
rsync -a --bwlimit=10000 ${2}* $REFERENCE_DIR 
rsync -a --bwlimit=500 $CURRENT_DIR/$1 $REFERENCE_DIR
/PHShome/adn5/programs/bin/bwa mem -C $REFERENCE_DIR/$FASTA_NAME $REFERENCE_DIR/$1|/apps/lib-osver/samtools/1.10/bin/samtools view -t $REFERENCE_DIR/$INDEX_NAME -F 0x800 -|gzip > $REFERENCE_DIR/${1}.sam.gz
mv $REFERENCE_DIR/${1}.sam.gz $4
rm $REFERENCE_DIR/${1}*

