#Xander Nuttle
#makejob_bwamem.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/makejob_bwamem.sh genome_fasta_file genome_size_file mapping_output_directory

NFASTQS=$(ls|grep fastq.gz$|wc -l)
ls|grep fastq.gz$ > fastq_files.txt
echo "#!/usr/bin/env bash" >> bwajobs.sh
echo "" >> bwajobs.sh
echo "#BSUB -J \"mapjobs[1-$NFASTQS]%500\"" >> bwajobs.sh
echo "#BSUB -L /bin/bash" >> bwajobs.sh
echo "#BSUB -o /dev/null" >> bwajobs.sh
echo "#BSUB -e /dev/null" >> bwajobs.sh
echo "#BSUB -R rusage[mem=2000]" >> bwajobs.sh
echo "#BSUB -W 1:0" >> bwajobs.sh
echo "#BSUB -cwd $(pwd)" >> bwajobs.sh
echo "#BSUB -sp 80" >> bwajobs.sh
echo "#BSUB -q medium" >> bwajobs.sh
echo "/data/talkowski/xander/MIPs/analysis_programs/map_bwamem.sh \`awk \"NR == \$LSB_JOBINDEX\" fastq_files.txt \` $1 $2 $3" >> bwajobs.sh

