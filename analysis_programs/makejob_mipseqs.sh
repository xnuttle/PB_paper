#Xander Nuttle
#makejob_mipseqs.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/makejob_mipseqs.sh miptargets_file

NSAMS=$(ls|grep sam.gz$|wc -l)
ls|grep sam.gz$ > samfiles.txt
echo "#!/usr/bin/env bash" >> mipseqjobs.sh
echo "" >> mipseqjobs.sh
echo "#BSUB -J \"mipseqjobs[1-$NSAMS]%500\"" >> mipseqjobs.sh
echo "#BSUB -L /bin/bash" >> mipseqjobs.sh
echo "#BSUB -o /dev/null" >> mipseqjobs.sh
echo "#BSUB -e /dev/null" >> mipseqjobs.sh
echo "#BSUB -R rusage[mem=2000]" >> mipseqjobs.sh
echo "#BSUB -W 1:0" >> mipseqjobs.sh
echo "#BSUB -cwd $(pwd)" >> mipseqjobs.sh
echo "#BSUB -sp 80" >> mipseqjobs.sh
echo "#BSUB -q medium" >> mipseqjobs.sh
echo "/data/talkowski/xander/MIPs/analysis_programs/get_mipseqs.sh \`awk \"NR == \$LSB_JOBINDEX\" samfiles.txt \` $1" >> mipseqjobs.sh

