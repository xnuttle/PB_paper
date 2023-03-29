#Xander Nuttle
#makejob_pear.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/makejob_pear.sh merging_output_directory

BARCODE_KEY=$(ls ..|grep barcodekey$)
for i in $(cut -f1 ../$BARCODE_KEY); do 
	ls|grep ^${i}_|grep FS1_F >> fastqs_F.txt
	ls|grep ^${i}_|grep FS1_R >> fastqs_R.txt
done
paste fastqs_F.txt fastqs_R.txt > fastqs.txt
NFASTQS=$(paste fastqs_F.txt fastqs_R.txt|wc -l)
echo "#!/usr/bin/env bash" >> pearjobs.sh
echo "" >> pearjobs.sh
echo "#BSUB -J \"pearjobs[1-$NFASTQS]%500\"" >> pearjobs.sh
echo "#BSUB -L /bin/bash" >> pearjobs.sh
echo "#BSUB -o /dev/null" >> pearjobs.sh
echo "#BSUB -e /dev/null" >> pearjobs.sh
echo "#BSUB -R rusage[mem=2000]" >> pearjobs.sh
echo "#BSUB -W 5:0" >> pearjobs.sh
echo "#BSUB -cwd $(pwd)" >> pearjobs.sh
echo "#BSUB -sp 80" >> pearjobs.sh
echo "#BSUB -q medium" >> pearjobs.sh
echo "/data/talkowski/xander/MIPs/analysis_programs/run_pear.sh \`awk \"NR == \$LSB_JOBINDEX\" fastqs.txt \` \$LSB_JOBINDEX $1" >> pearjobs.sh

