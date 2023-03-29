#Xander Nuttle
#get_guides_tags.sh

echo -e "Sample\tGuideInt\tGuideTag\tMipCount\tMipFrac" > PB_megapool_7I.guidetags
for i in $(cut -f1 ../final_results/indel_only_PB_megapool_7I.barcodekey); do /data/talkowski/xander/MIPs/analysis_programs/get_guides_tags ${i}.dp10.af0.1.finalseqs.gz tagged_indel_guides.txt >> PB_megapool_7I.guidetags; done

