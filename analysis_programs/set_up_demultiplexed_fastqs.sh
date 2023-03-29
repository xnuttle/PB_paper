#Xander Nuttle
#set_up_demultiplexed_fastqs.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/set_up_demultiplexed_fastqs.sh barcodekey_file datadir_file
#
#Creates soft links to demultiplexed fastq files. Also generates barcodekey files for individual samples and a file listing
#all samples using pseudonames (sample0001, sample0002, ...). The resulting links and files can then be used for extracting
#molecular tags and retaining only reads where both corresponding index reads perfectly match that sample's barcodes. The input
#"barcodekey_file" is a tab-delimited file listing sample names in column 1, index barcodes in column 2, and index2 barcodes
#in column 3. The datadir file is a text file containing the path to the directory containing original gzipped fastq files
#downloaded from Broad.

num=1
datadir=$(cat $2)
for i in $(cut -f1 $1);
	do barcode=`grep $i $1|cut -f2-3|sed 's/\t/_/g'`;
	setnum=`echo $num|awk '{printf "%04s\n", $1}'`;
	setname="sample$setnum";
	num=`echo "$num+1"|bc`; 
	fastq1=`ls $datadir|grep '\.1\.fastq'|grep $barcode`
	fastq2=`ls $datadir|grep '\.barcode_1\.fastq'|grep $barcode`
	fastq3=`ls $datadir|grep '\.barcode_2\.fastq'|grep $barcode`
	fastq4=`ls $datadir|grep '\.2\.fastq'|grep $barcode`
	ln -s $datadir/$fastq1 ${setname}.r1.fastq.gz
	ln -s $datadir/$fastq2 ${setname}.bc1.fastq.gz
	ln -s $datadir/$fastq3 ${setname}.bc2.fastq.gz
	ln -s $datadir/$fastq4 ${setname}.r2.fastq.gz
	newbarcode=`echo $barcode|sed 's/_//g'`
	echo -e "$i\t$newbarcode" > ${setname}.barcodekey
	echo $setname >> samplesets.txt
done

