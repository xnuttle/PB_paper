#Xander Nuttle
#makejob_fastq_prep_si.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/makejob_fastq_prep_si.sh merging_input_directory molecular_tag_length
#
#Differs from makejob_fastq_prep.sh in that only single indexing was used for barcoding rather than dual indexing.
#Also maximally allows 50 jobs to run at same time, as old value (500) was resulting in a lot of failures (sometimes not given exit code 1).

NSETS=$(cat samplesets.txt|wc -l)
echo "#!/usr/bin/env bash" >> prepjobs.sh
echo "" >> prepjobs.sh
echo "#BSUB -J \"prepjobs[1-$NSETS]%50\"" >> prepjobs.sh
echo "#BSUB -L /bin/bash" >> prepjobs.sh
echo "#BSUB -o /dev/null" >> prepjobs.sh
echo "#BSUB -e /dev/null" >> prepjobs.sh
echo "#BSUB -R rusage[mem=2000]" >> prepjobs.sh
echo "#BSUB -W 5:0" >> prepjobs.sh
echo "#BSUB -cwd $(pwd)" >> prepjobs.sh
echo "#BSUB -sp 80" >> prepjobs.sh
echo "#BSUB -q medium" >> prepjobs.sh
echo "/data/talkowski/xander/MIPs/analysis_programs/prep_fastqs_si.sh \`awk \"NR == \$LSB_JOBINDEX\" samplesets.txt \` $1 $2" >> prepjobs.sh

