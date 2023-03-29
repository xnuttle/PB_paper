#Xander Nuttle
#makejob_callcn.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/makejob_callcn.sh

grep -v Contig *mipcounts|cut -f2|sort|uniq > cncontigs.txt
NCHR=$(wc -l cncontigs.txt|awk '{print $1}')
echo "#!/usr/bin/env bash" >> cnjobs.sh
echo "" >> cnjobs.sh
echo "#BSUB -J \"cnjobs[1-$NCHR]%500\"" >> cnjobs.sh
echo "#BSUB -L /bin/bash" >> cnjobs.sh
echo "#BSUB -o /dev/null" >> cnjobs.sh
echo "#BSUB -e /dev/null" >> cnjobs.sh
echo "#BSUB -R rusage[mem=2000]" >> cnjobs.sh
echo "#BSUB -W 5:0" >> cnjobs.sh
echo "#BSUB -cwd $(pwd)" >> cnjobs.sh
echo "#BSUB -sp 80" >> cnjobs.sh
echo "#BSUB -q medium" >> cnjobs.sh
echo "/data/talkowski/xander/MIPs/analysis_programs/callcn.sh \`awk \"NR == \$LSB_JOBINDEX\" cncontigs.txt \`" >> cnjobs.sh

