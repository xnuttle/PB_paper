#!/usr/bin/env bash

#Xander Nuttle
#callcn.sh
#Call: /data/talkowski/xander/MIPs/analysis_programs/callcn.sh contig_name

REFERENCE_DIR=/var/tmp/xnuttle
CURRENT_DIR=$(pwd)
PROGRAM_DIR=/data/talkowski/xander/MIPs/analysis_programs
exptname=$(basename `dirname $CURRENT_DIR`)
mcounts=${exptname}.mipcounts
hydin=$(echo ${1}|grep HYDIN|wc -l)

mkdir -p $REFERENCE_DIR/
chgrp -R miket $REFERENCE_DIR/
chmod -R g+wx $REFERENCE_DIR/
rsync -a --bwlimit=500 $CURRENT_DIR/$mcounts $REFERENCE_DIR/

cd $REFERENCE_DIR
for samp in $(grep ${1} ${mcounts}|cut -f1|uniq); do
	head -1 $mcounts > ${samp}_${1}.mipcounts
	echo > ${samp}_${1}.miptargets.v3
	for targ in $(grep ${1} ${mcounts}|grep ${samp}|sed 's/\t/:/g'); do
		echo $targ|sed 's/:/\t/g' >> ${samp}_${1}.mipcounts
		echo $targ|awk -F : '{print $2,$3","$3,$1,"S 0 0 AB + 20"}' >> ${samp}_${1}.miptargets.v3
	done
	if [ $hydin -eq 1 ]; then
		$PROGRAM_DIR/call_mip_pscn ${samp}_${1}.miptargets.v3 ${samp}_${1}.mipcounts 4
	else
		$PROGRAM_DIR/call_mip_hapcn ${samp}_${1}.miptargets.v3 ${samp}_${1}.mipcounts 2
	fi
done

echo -e "Individual\tHap1_CN\tHap2_CN\tLOD_Score\tPossible_Complex_CN_Genotype" > $CURRENT_DIR/${exptname}_${1}.cncalls
cat *cncalls|grep -v Individual >> $CURRENT_DIR/${exptname}_${1}.cncalls
cat *compevents > $CURRENT_DIR/${exptname}_${1}.compevents
cat *simplecalls > $CURRENT_DIR/${exptname}_${1}.simplecalls
rm *${1}.mipcounts
rm *${1}.miptargets.v3
rm *${1}.cncalls
rm *${1}.compevents
rm *${1}.simplecalls

