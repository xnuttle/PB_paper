#Xander Nuttle
#mrmip_pb9_dm_fastq.sh
#Call: bsub -sla miket_sc -q long -R 'hname!=ib024 && hname!=ib028 && hname!=ib031 && hname!=cmu103 && hname!=cmu105 && hname!=cmu234' -J mrmip_pb9_dm_fastq < /data/talkowski/xander/MIPs/analysis_programs/mrmip_pb9_dm_fastq.sh
#
#The script automates the MIP analysis pipeline for genotyping for use in analyzing data from PB megapool 9, with
#the starting sequencing data in demultiplexed gzipped fastq file format. It should be run from the master directory for an experiment.

#instruct shell to exit immediately if any step in pipeline fails
set -e

#set up variable specifying nodes to exclude from executing jobs
avoid_badnodes="-R 'hname!=ib022 && hname!=ib024 && hname!=ib028 && hname!=ib031 && hname!=cmu103 && hname!=cmu105 && hname!=cmu234'"

#set up directory variables
EXP_DIR=$(pwd)
FASTQ_DIR=$EXP_DIR/raw_fastq_files
LOG_DIR=$EXP_DIR/log
PEAR_IN_DIR=$EXP_DIR/pear_input
PEAR_OUT_DIR=$EXP_DIR/pear_output
MAP_IN_DIR=$EXP_DIR/pear_output
MAP_OUT_DIR=$EXP_DIR/bwa_mapping_output
OUT_DIR=$EXP_DIR/final_results
PROGRAM_DIR=/data/talkowski/xander/MIPs/analysis_programs

#set up other variables
barcodefile=`ls $EXP_DIR|grep barcodekey`
datadirfile=`ls $EXP_DIR|grep datadir`
taglength=8
genome=/data/talkowski/xander/MIPs/genomes/PB_indels_and_RGDs/PB_indels_and_RGDs.fasta
gsize=/data/talkowski/xander/MIPs/genomes/PB_indels_and_RGDs/PB_indels_and_RGDs.fasta.fai
miptargets=/data/talkowski/xander/MIPs/genomes/PB_indels_and_RGDs/PB_indels_and_RGDs.miptargets
experiment=$(basename $EXP_DIR)
crisprinfo=/data/talkowski/xander/MIPs/genomes/PB_indels_and_RGDs/PB_indels_and_RGDs.crispr
ptargs=/data/talkowski/xander/MIPs/genomes/PB_indels_and_RGDs/miptargets/chrPLASMIDS.miptargets
pbcode=/data/talkowski/xander/MIPs/genomes/PB_indels_and_RGDs/PB_indels_and_RGDs.pbcode
gtargs=/data/talkowski/xander/MIPs/genomes/PB_indels_and_RGDs/guides.miptargets

#make fastq directory, log directory (if not present), and pear input directory
mkdir -p $FASTQ_DIR $LOG_DIR $PEAR_IN_DIR
echo "Made fastq, log, and pear input directories." > $LOG_DIR/mrmip_pb_dm_fastq.log

#PREPARE READS FOR PEAR
cd $FASTQ_DIR
bash $PROGRAM_DIR/set_up_demultiplexed_fastqs_pb9.sh $EXP_DIR/$barcodefile $EXP_DIR/$datadirfile
bash $PROGRAM_DIR/makejob_fastq_prep_si.sh $PEAR_IN_DIR $taglength
echo "bsub -sla miket_sc $avoid_badnodes < prepjobs.sh"|bash
echo "Cluster fastq prep job submitted." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#monitor status of cluster fastq prep job periodically and do not proceed until it finishes
prepjob=$(bjobs|grep bs|head -1|awk '{print $1}')
prepping=$(bjobs|grep bs|wc -l)
while [ $prepping -gt 0 ]; do
	echo "Cluster fastq prep job not complete, waiting..." >> $LOG_DIR/mrmip_pb_dm_fastq.log
	sleep 5m
	prepping=$(bjobs|grep bs|wc -l)
done
echo "Cluster fastq prep job complete." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#wait for cluster fastq prep job report to be generated
sleep 5m

#ensure cluster fastq prep job exited successfully
nsets=`wc -l samplesets.txt|awk '{print $1}'`
gsets=`bhist -l $prepjob|grep Done|wc -l`
if [ $gsets -eq $nsets ]; then
	echo "CLUSTER FASTQ PREP SUCCESSFUL!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
	echo "CLUSTER FASTQ PREP FAILED" >> $LOG_DIR/mrmip_pb_dm_fastq.log
	exit 1
fi
echo "Cluster fastq prep job exited successfully." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#ensure gzipped fastq files have good integrity
cd $PEAR_IN_DIR
fqcorrupt=0
for file in $(ls|grep 'fastq.gz'); do
	gzip -t $file
	fqcorrupt=`echo "$?+$fqcorrupt"|bc`
done
if [ $fqcorrupt -eq 0 ]; then
	echo "ALL GZIPPED FASTQ FILES GOOD!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
	echo "AT LEAST ONE GZIPPED FASTQ FILE BAD" >> $LOG_DIR/mrmip_pb_dm_fastq.log
	exit 1
fi
echo "Gzipped fastq files have good integrity." >> $LOG_DIR/mrmip_pb_dm_fastq.log
echo "FASTQ PREP COMPLETE, READS READY FOR PEAR: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log

#MERGE READS WITH PEAR
#make pear output directory
cd $EXP_DIR
mkdir -p $PEAR_OUT_DIR
echo "Made pear output directory." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#merge reads using PEAR
cd $PEAR_IN_DIR
bash $PROGRAM_DIR/makejob_pear.sh $PEAR_OUT_DIR
echo "bsub -sla miket_sc $avoid_badnodes < pearjobs.sh"|bash
echo "Cluster merging job submitted." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#monitor status of cluster merging job periodically and do not proceed until it finishes
mergejob=$(bjobs|grep bs|head -1|awk '{print $1}')
merging=$(bjobs|grep bs|wc -l)
while [ $merging -gt 0 ]; do
  echo "Cluster merging job not complete, waiting..." >> $LOG_DIR/mrmip_pb_dm_fastq.log
  sleep 5m
  merging=$(bjobs|grep bs|wc -l)
done
echo "Cluster merging job complete." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#wait for cluster merging job report to be generated
sleep 5m

#ensure cluster merging job exited successfully
nmergings=`wc -l fastqs.txt|awk '{print $1}'`
gmergings=`bhist -l $mergejob|grep Done|wc -l`
if [ $gmergings -eq $nmergings ]; then
  echo "CLUSTER MERGING SUCCESSFUL!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
  echo "CLUSTER MERGING FAILED" >> $LOG_DIR/mrmip_pb_dm_fastq.log
  exit 1
fi
echo "Cluster merging job exited successfully." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#ensure gzipped merged fastq files have good integrity
cd $PEAR_OUT_DIR
fqcorrupt=0
for file in $(ls|grep 'fastq.gz'); do
  gzip -t $file
  fqcorrupt=`echo "$?+$fqcorrupt"|bc`
done
if [ $fqcorrupt -eq 0 ]; then
  echo "ALL GZIPPED FASTQ FILES GOOD!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
  echo "AT LEAST ONE GZIPPED FASTQ FILE BAD" >> $LOG_DIR/mrmip_pb_dm_fastq.log
  exit 1
fi
echo "Gzipped fastq files have good integrity." >> $LOG_DIR/mrmip_pb_dm_fastq.log
echo "MERGING WITH PEAR COMPLETE: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log

#MAP MERGED FASTQ BWA
#make mapping output directory
cd $EXP_DIR
mkdir -p $MAP_OUT_DIR
echo "Made mapping output directory." > $LOG_DIR/mrmip_pb_dm_fastq.log

#map reads with bwa mem
cd $MAP_IN_DIR
bash $PROGRAM_DIR/makejob_bwamem.sh $genome $gsize $MAP_OUT_DIR
echo "bsub -sla miket_sc $avoid_badnodes < bwajobs.sh"|bash
echo "Cluster mapping job submitted." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#monitor status of cluster mapping job periodically and do not proceed until it finishes
mapjob=$(bjobs|grep bs|head -1|awk '{print $1}')
mapping=$(bjobs|grep bs|wc -l)
while [ $mapping -gt 0 ]; do
  echo "Cluster mapping job not complete, waiting..." >> $LOG_DIR/mrmip_pb_dm_fastq.log
  sleep 5m
  mapping=$(bjobs|grep bs|wc -l)
done
echo "Cluster mapping job complete." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#wait for cluster mapping job report to be generated
sleep 5m

#ensure cluster mapping job exited successfully
nmappings=`wc -l fastq_files.txt|awk '{print $1}'`
gmappings=`bhist -l $mapjob|grep Done|wc -l`
if [ $gmappings -eq $nmappings ]; then
  echo "CLUSTER MAPPING SUCCESSFUL!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
  echo "CLUSTER MAPPING FAILED" >> $LOG_DIR/mrmip_pb_dm_fastq.log
  exit 1
fi
echo "Cluster mapping job exited successfully." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#ensure mapping output sam files have good integrity
cd $MAP_OUT_DIR
samcorrupt=0
for file in $(ls|grep 'fastq.gz.sam.gz'); do
  gzip -t $file
	samcorrupt=`echo "$?+$samcorrupt"|bc`
done
if [ $samcorrupt -eq 0 ]; then
  echo "ALL GZIPPED SAM FILES GOOD!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
  echo "AT LEAST ONE GZIPPED SAM FILE BAD" >> $LOG_DIR/mrmip_pb_dm_fastq.log
  exit 1
fi
echo "Gzipped sam files have good integrity." >> $LOG_DIR/mrmip_pb_dm_fastq.log
echo "CLUSTER MAPPING COMPLETE: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log

#MAPPING OUTPUT TO FINAL MIPSEQS
#generate mip sequence files from mapping output
bash $PROGRAM_DIR/makejob_mipseqs.sh $miptargets
echo "bsub -sla miket_sc $avoid_badnodes < mipseqjobs.sh"|bash
echo "Cluster mipseqs job submitted." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#monitor status of cluster mipseqs job periodically and do not proceed until it finishes
mseqjob=$(bjobs|grep bs|head -1|awk '{print $1}')
mseqing=$(bjobs|grep bs|wc -l)
while [ $mseqing -gt 0 ]; do
  echo "Cluster mipseqs job not complete, waiting..." >> $LOG_DIR/mrmip_pb_dm_fastq.log
  sleep 5m
  mseqing=$(bjobs|grep bs|wc -l)
done
echo "Cluster mipseqs job complete." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#wait for cluster mipseqs job report to be generated
sleep 5m

#ensure cluster mipseqs job exited successfully
nmipseqs=`wc -l samfiles.txt|awk '{print $1}'`
gmipseqs=`bhist -l $mseqjob|grep Done|wc -l`
if [ $gmipseqs -eq $nmipseqs ]; then
  echo "CLUSTER MIPSEQS SUCCESSFUL!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
  echo "CLUSTER MIPSEQS FAILED" >> $LOG_DIR/mrmip_pb_dm_fastq.log
  exit 1
fi
echo "Cluster mipseqs job exited successfully." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#ensure gzipped mipseqs files have good integrity
mscorrupt=0
for file in $(ls|grep 'mipseqs.gz'); do
  gzip -t $file
  mscorrupt=`echo "$?+$mscorrupt"|bc`
done
if [ $mscorrupt -eq 0 ]; then
  echo "ALL GZIPPED MIPSEQS FILES GOOD!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
  echo "AT LEAST ONE GZIPPED MIPSEQS FILE BAD" >> $LOG_DIR/mrmip_pb_dm_fastq.log
  exit 1
fi
echo "Gzipped mipseqs files have good integrity." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#count mip sequences for each sample, call final mip sequences, and generate mipcounts files for each sample
bash $PROGRAM_DIR/makejob_seqcounts.sh $miptargets
echo "bsub -sla miket_sc $avoid_badnodes < seqcountjobs.sh"|bash
echo "Cluster seqcount job submitted." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#monitor status of cluster seqcount job periodically and do not proceed until it finishes
seqctjob=$(bjobs|grep bs|head -1|awk '{print $1}')
seqcting=$(bjobs|grep bs|wc -l)
while [ $seqcting -gt 0 ]; do
  echo "Cluster seqcount job not complete, waiting..." >> $LOG_DIR/mrmip_pb_dm_fastq.log
  sleep 5m
  seqcting=$(bjobs|grep bs|wc -l)
done
echo "Cluster seqcount job complete." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#wait for cluster seqcount job report to be generated
sleep 5m

#ensure cluster seqcount job exited successfully
nscjobs=`wc -l masterlist.txt|awk '{print $1}'`
gscjobs=`bhist -l $seqctjob|grep Done|wc -l`
if [ $gscjobs -eq $nscjobs ]; then
  echo "CLUSTER SEQCOUNT SUCCESSFUL!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
  echo "CLUSTER SEQCOUNT FAILED" >> $LOG_DIR/mrmip_pb_dm_fastq.log
  exit 1
fi
echo "Cluster seqcount job exited successfully." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#ensure gzipped seqcounts files have good integrity
sccorrupt=0
for file in $(ls|grep 'seqcounts.gz'); do
  gzip -t $file
  sccorrupt=`echo "$?+$sccorrupt"|bc`
done
if [ $sccorrupt -eq 0 ]; then
  echo "ALL GZIPPED SEQCOUNTS FILES GOOD!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
  echo "AT LEAST ONE GZIPPED SEQCOUNTS FILE BAD" >> $LOG_DIR/mrmip_pb_dm_fastq.log
  exit 1
fi
echo "Gzipped seqcounts files have good integrity." >> $LOG_DIR/mrmip_pb_dm_fastq.log
echo "MIPSEQS AND SEQCOUNTS ANALYSES COMPLETE: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log

#CALL CRISPR SEQUENCE EDITS
cd $EXP_DIR
mkdir -p $OUT_DIR
echo "Made final results directory." >> $LOG_DIR/mrmip_pb_dm_fastq.log
cd $MAP_OUT_DIR
echo -e "Sample\tCRISPR\tGenotype\tIndelAlns" > ${experiment}.crisprcalls
for i in $(cut -f1 ../$barcodefile); do $PROGRAM_DIR/call_crispr_vars ${i}.dp10.af0.1.finalseqs.gz $crisprinfo >> ${experiment}.crisprcalls; done
mv ${experiment}.crisprcalls $OUT_DIR
echo "SEQUENCE GENOTYPING COMPLETE: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log

#CALL PB INTEGRATION STATUSES
echo -e "Sample\tIntStatus" > ${experiment}.intstats
for i in $(cut -f1 ../$barcodefile); do $PROGRAM_DIR/call_pb_int_status ${i}.dp10.af0.1.finalseqs.gz $ptargs $pbcode >> ${experiment}.intstats; done
mv ${experiment}.intstats $OUT_DIR
echo "PB INTEGRATION GENOTYPING COMPLETE: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log

#GET COUNTS OF MIP CAPTURE EVENTS FOR EACH GUIDE CONSTRUCT, CALL PB INTEGRATION COPY NUMBERS, AND GENERATE FILE LISTING INTEGRATED GUIDE CONSTRUCTS
echo -e "Sample\tPBCN" > ${experiment}.pbcounts
echo -e "Sample\tGuide" > ${experiment}.pbguides
for i in $(cut -f1 ../$barcodefile); do
	$PROGRAM_DIR/get_guidecounts ${i}.dp10.af0.1.finalseqs.gz $gtargs
	$PROGRAM_DIR/call_pb_cn ${i}.guidecounts >> ${experiment}.pbcounts
	$PROGRAM_DIR/get_pb_guides ${i}.guidecounts >> ${experiment}.pbguides
done
mkdir $OUT_DIR/guidecounts
mv *guidecounts $OUT_DIR/guidecounts
mv ${experiment}.pbcounts $OUT_DIR
mv ${experiment}.pbguides $OUT_DIR
echo "GUIDE CONSTRUCT GENOTYPING COMPLETE: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log

#CALL CRISPR COPY NUMBER EDITS
#combine individual mipcounts output files into three large mipcounts files (one for all for plotting, one for MGH2069, and one for GM08330)
echo -e "Sample\tContig\tCoordinate\tMip_Type\tHaplotype_1_Count\tHaplotype_2_Count" > ${experiment}.mipcounts
for i in $(cut -f1 ${EXP_DIR}/${barcodefile}); do
	{ grep -v Contig ${i}.mipcounts >> ${experiment}.mipcounts || true; }
done
echo "Combined mipcounts file generation finished." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#run the automated copy number genotyping caller for each region where copy number was interrogated
bash $PROGRAM_DIR/makejob_callcn.sh
echo "bsub -sla miket_sc $avoid_badnodes < cnjobs.sh"|bash
echo "Cluster copy number calling job submitted." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#monitor status of cluster copy number calling job periodically and do not proceed until it finishes
cnjob=$(bjobs|grep bs|head -1|awk '{print $1}')
cncalling=$(bjobs|grep bs|wc -l)
while [ $cncalling -gt 0 ]; do
  echo "Cluster copy number calling job not complete, waiting..." >> $LOG_DIR/mrmip_pb_dm_fastq.log
  sleep 5m
  cncalling=$(bjobs|grep bs|wc -l)
done
echo "Cluster copy number calling job complete." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#wait for cluster copy number calling job report to be generated
sleep 5m

#ensure cluster copy number calling job exited successfully
ncnjobs=`wc -l cncontigs.txt|awk '{print $1}'`
gcnjobs=`bhist -l $cnjob|grep Done|wc -l`
if [ $gcnjobs -eq $ncnjobs ]; then
  echo "CLUSTER COPY NUMBER CALLING SUCCESSFUL!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
  echo "CLUSTER COPY NUMBER CALLING FAILED" >> $LOG_DIR/mrmip_pb_dm_fastq.log
  exit 1
fi
echo "Cluster copy number calling job exited successfully." >> $LOG_DIR/mrmip_pb_dm_fastq.log
echo "COPY NUMBER GENOTYPING COMPLETE: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log

#move mipcounts files, copy number caller output files, and barcodekey file to final output directory
mv *mipcounts $OUT_DIR
mv *cncalls $OUT_DIR
mv *compevents $OUT_DIR
mv *simplecalls $OUT_DIR
mv $EXP_DIR/$barcodefile $OUT_DIR

#GENERATE PLOTS
#for each region where copy number was interrogated, generate plots for each individual showing read count frequencies for the two most abundant sequences at each MIP
cd $OUT_DIR
bash $PROGRAM_DIR/makejob_plot_pb.sh
echo "bsub -sla miket_sc $avoid_badnodes < plotjobs.sh"|bash
echo "Cluster plotting job submitted." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#monitor status of cluster plotting job periodically and do not proceed until it finishes
plotjob=$(bjobs|grep bs|head -1|awk '{print $1}')
plotting=$(bjobs|grep bs|wc -l)
while [ $plotting -gt 0 ]; do
	echo "Cluster plotting job not complete, waiting..." >> $LOG_DIR/mrmip_pb_dm_fastq.log
	sleep 5m
	plotting=$(bjobs|grep bs|wc -l)
done
echo "Cluster plotting job complete." >> $LOG_DIR/mrmip_pb_dm_fastq.log

#wait for cluster plotting job report to be generated
sleep 5m

#ensure cluster plotting job exited successfully
nplotted=`wc -l genes.txt|awk '{print $1}'`
gplotted=`bhist -l $plotjob|grep Done|wc -l`
if [ $gplotted -eq $nplotted ]; then
	echo "CLUSTER PLOTTING SUCCESSFUL!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
else
	echo "CLUSTER PLOTTING FAILED" >> $LOG_DIR/mrmip_pb_dm_fastq.log
	exit 1
fi
echo "Cluster plotting job exited successfully." >> $LOG_DIR/mrmip_pb_dm_fastq.log
echo "PLOTTING COMPLETE: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log

#clean up and return success
rm -rf $FASTQ_DIR
echo "PIPELINE COMPLETE: SUCCESS!" >> $LOG_DIR/mrmip_pb_dm_fastq.log
exit 0

