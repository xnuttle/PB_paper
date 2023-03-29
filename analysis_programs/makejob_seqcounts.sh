#Xander Nuttle
#makejob_seqcounts.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/makejob_seqcounts.sh miptargets_file

BARCODE_KEY=$(ls ..|grep barcodekey$)
for i in $(cut -f1 ../$BARCODE_KEY); do
	ls|grep ^${i}_|grep mipseqs.gz$ > ${i}.seqsfiles
done
NSAMPS=$(ls|grep seqsfiles$|wc -l)
ls|grep seqsfiles$ > masterlist.txt
echo "#!/usr/bin/env bash" >> seqcountjobs.sh
echo "" >> seqcountjobs.sh
echo "#BSUB -J \"seqcountjobs[1-$NSAMPS]%500\"" >> seqcountjobs.sh
echo "#BSUB -L /bin/bash" >> seqcountjobs.sh
echo "#BSUB -o /dev/null" >> seqcountjobs.sh
echo "#BSUB -e /dev/null" >> seqcountjobs.sh
echo "#BSUB -R rusage[mem=2000]" >> seqcountjobs.sh
echo "#BSUB -W 1:0" >> seqcountjobs.sh
echo "#BSUB -cwd $(pwd)" >> seqcountjobs.sh
echo "#BSUB -sp 80" >> seqcountjobs.sh
echo "#BSUB -q medium" >> seqcountjobs.sh
echo "/data/talkowski/xander/MIPs/analysis_programs/process_mipseqs.sh \`awk \"NR == \$LSB_JOBINDEX\" masterlist.txt \` $1" >> seqcountjobs.sh

