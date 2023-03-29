#Xander Nuttle
#makejob_plot_pb.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/makejob_plot_pb.sh

grep -v Contig *mipcounts|cut -f2|sort|uniq > genes.txt
NGENES=$(wc -l genes.txt|awk '{print $1}')
echo "#!/usr/bin/env bash" >> plotjobs.sh
echo "" >> plotjobs.sh
echo "#BSUB -J \"plotjobs[1-$NGENES]%500\"" >> plotjobs.sh
echo "#BSUB -L /bin/bash" >> plotjobs.sh
echo "#BSUB -o /dev/null" >> plotjobs.sh
echo "#BSUB -e /dev/null" >> plotjobs.sh
echo "#BSUB -R rusage[mem=2000]" >> plotjobs.sh
echo "#BSUB -W 5:0" >> plotjobs.sh
echo "#BSUB -cwd $(pwd)" >> plotjobs.sh
echo "#BSUB -sp 80" >> plotjobs.sh
echo "#BSUB -q medium" >> plotjobs.sh
echo "/data/talkowski/xander/MIPs/analysis_programs/plot_pb.sh \`awk \"NR == \$LSB_JOBINDEX\" genes.txt \`" >> plotjobs.sh

