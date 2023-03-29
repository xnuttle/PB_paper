#!/usr/bin/env bash

#Xander Nuttle
#plot_pb.sh
#Call: /data/talkowski/xander/MIPs/analysis_programs/plot_pb.sh GENE_NAME 

REFERENCE_DIR=/var/tmp/xnuttle
CURRENT_DIR=$(pwd)
SCRIPT_DIR=/data/talkowski/xander/MIPs/analysis_programs
GLOCS_DIR=/data/talkowski/xander/MIPs/genomes/PBINP2C3/guidelocs
guides=$(ls $GLOCS_DIR|grep $1|wc -l)
exptname=$(basename `dirname $CURRENT_DIR`)
mipcounts=${exptname}.mipcounts
barcodekey=${exptname}.barcodekey

mkdir -p $REFERENCE_DIR/
chgrp -R miket $REFERENCE_DIR/
chmod -R g+wx $REFERENCE_DIR/
rsync -a --bwlimit=500 $CURRENT_DIR/$mipcounts $REFERENCE_DIR/
rsync -a --bwlimit=500 $CURRENT_DIR/$barcodekey $REFERENCE_DIR/
if [ $guides -eq 1 ]; then
	rsync -a --bwlimit=500 $GLOCS_DIR/${1}.guidelocs $REFERENCE_DIR/
fi
cd $REFERENCE_DIR
head -1 $mipcounts > ${exptname}_${1}.mipcounts
grep -P "\t$1\t" $mipcounts >> ${exptname}_${1}.mipcounts
module load R/3.3.0
unset R_HOME
Rscript $SCRIPT_DIR/pdf_pb_mips.r $exptname $1
mv ${exptname}_${1}.pdf $CURRENT_DIR
rm $REFERENCE_DIR/${exptname}_${1}.mipcounts
if [ $guides -eq 1 ]; then
	rm $REFERENCE_DIR/${1}.guidelocs
fi

