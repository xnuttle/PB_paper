#!/usr/bin/env bash

#Xander Nuttle
#process_mipseqs.sh
#Call: /data/talkowski/xander/MIPs/analysis_programs/process_mipseqs.sh text_file_listing_mipseqs_files miptargets_file

REFERENCE_DIR=/var/tmp/xnuttle
PROGRAM_DIR=/data/talkowski/xander/MIPs/analysis_programs
CURRENT_DIR=$(pwd)
SAMP_NAME=$(basename ${1} .seqsfiles)
MTARGS_NAME=$(basename ${2})
HYDIN_TARGS=/data/talkowski/xander/MIPs/genomes/PBINP2C3/miptargets/chrHYDIN.miptargets
HYDIN2_TARGS=/data/talkowski/xander/MIPs/genomes/PBINP2C3/miptargets/chrHYDIN2.miptargets
HTARGS=$(basename ${HYDIN_TARGS})
H2TARGS=$(basename ${HYDIN2_TARGS})

mkdir -p $REFERENCE_DIR
chgrp -R miket $REFERENCE_DIR
rsync -a --bwlimit=500 $CURRENT_DIR/$1 $REFERENCE_DIR
rsync -a --bwlimit=500 $2 $REFERENCE_DIR
rsync -a --bwlimit=500 ${SAMP_NAME}*mipseqs.gz $REFERENCE_DIR
rsync -a --bwlimit=500 $HYDIN_TARGS $REFERENCE_DIR
rsync -a --bwlimit=500 $HYDIN2_TARGS $REFERENCE_DIR
cd $REFERENCE_DIR
$PROGRAM_DIR/count_mipseqs $1 $MTARGS_NAME
$PROGRAM_DIR/finalize_mipseqs ${SAMP_NAME}.seqcounts.gz 10 0.1
$PROGRAM_DIR/finalseqs_to_mipcounts ${SAMP_NAME}.dp10.af0.1.finalseqs.gz

grep -v HYDIN ${SAMP_NAME}.mipcounts > ${SAMP_NAME}.temp.mipcounts
for mip in $(cut -f1 ${HTARGS}|grep -v Name); do
	mip2=`echo $mip|sed 's/HYDIN/HYDIN2/g'`
	coord1=`grep $mip $HTARGS|awk '{print $4,$9}'|sed 's/ /+/g'|bc`
	coord2=`grep $mip2 $H2TARGS|awk '{print $4,$9}'|sed 's/ /+/g'|bc`
	echo -e "$SAMP_NAME\tchrHYDIN\t$coord1\tB\t" >> ${SAMP_NAME}.aaa
	count1=`grep chrHYDIN ${SAMP_NAME}.mipcounts|grep -v HYDIN2|grep -P "chrHYDIN\t$coord1\t"|awk '{print $5,$6}'|sed 's/ /+/g'|bc`
	count2=`grep chrHYDIN2 ${SAMP_NAME}.mipcounts|grep -P "chrHYDIN2\t$coord2\t"|awk '{print $5,$6}'|sed 's/ /+/g'|bc`
	good1=`echo $count1|wc -L`
	good2=`echo $count2|wc -L`
	if [ $good1 -gt 0 ]; then echo $count1 >> ${SAMP_NAME}.bbb; else echo "0" >> ${SAMP_NAME}.bbb; fi
	if [ $good2 -gt 0 ]; then echo $count2 >> ${SAMP_NAME}.ccc; else echo "0" >> ${SAMP_NAME}.ccc; fi
done

cp ${SAMP_NAME}.temp.mipcounts ${SAMP_NAME}.mipcounts
paste ${SAMP_NAME}.aaa ${SAMP_NAME}.bbb ${SAMP_NAME}.ccc >> ${SAMP_NAME}.mipcounts

mv $REFERENCE_DIR/${SAMP_NAME}.seqcounts.gz $CURRENT_DIR
mv $REFERENCE_DIR/${SAMP_NAME}.dp10.af0.1.finalseqs.gz $CURRENT_DIR
mv $REFERENCE_DIR/${SAMP_NAME}.mipcounts $CURRENT_DIR
rm $REFERENCE_DIR/${SAMP_NAME}*

