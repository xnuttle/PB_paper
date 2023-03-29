#Xander Nuttle
#makejob_fastq_prep.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/makejob_fastq_prep.sh merging_input_directory molecular_tag_length

NSETS=$(cat samplesets.txt|wc -l)
echo "#!/usr/bin/env bash" >> prepjobs.sh
echo "" >> prepjobs.sh
echo "#BSUB -J \"prepjobs[1-$NSETS]%500\"" >> prepjobs.sh
echo "#BSUB -L /bin/bash" >> prepjobs.sh
echo "#BSUB -o /dev/null" >> prepjobs.sh
echo "#BSUB -e /dev/null" >> prepjobs.sh
echo "#BSUB -R rusage[mem=2000]" >> prepjobs.sh
echo "#BSUB -W 5:0" >> prepjobs.sh
echo "#BSUB -cwd $(pwd)" >> prepjobs.sh
echo "#BSUB -sp 80" >> prepjobs.sh
echo "#BSUB -q medium" >> prepjobs.sh
echo "/data/talkowski/xander/MIPs/analysis_programs/prep_fastqs.sh \`awk \"NR == \$LSB_JOBINDEX\" samplesets.txt \` $1 $2" >> prepjobs.sh

