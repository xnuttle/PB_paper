#Xander Nuttle
#set_up_demultiplexed_fastqs_pb9.sh
#Call: bash /data/talkowski/xander/MIPs/analysis_programs/set_up_demultiplexed_fastqs_pb9.sh barcodekey_file datadir_file
#
#Differs from set_up_demultiplexed_fastqs.sh in that only single indexing was used for barcoding rather than dual indexing
#and accommodates PB_megapool_9's having each sample sequenced in two separate lanes (i.e., having two 'read1' files, two
#'read 2' files, and two 'index read' files per sample.
#
#Creates demultiplexed fastq files with all reads from both lanes merged. Also generates barcodekey files for individual samples
#and a file listing all samples using pseudonames (sample0001, sample0002, ...). The resulting links and files can then be used
#for extracting molecular tags and retaining only reads where corresponding index reads perfectly match that sample's barcode.
#The input "barcodekey_file" is a tab-delimited file listing sample names in column 1 and index barcodes in column 2. The datadir
#file is a text file containing the path to the directory containing original gzipped fastq files downloaded from Broad.

num=1
datadir=$(cat $2)
for i in $(cut -f1 $1);
	do barcode=`grep $i $1|cut -f2`;
	setnum=`echo $num|awk '{printf "%04s\n", $1}'`;
	setname="sample$setnum";
	num=`echo "$num+1"|bc`; 
	for j in $(ls $datadir|grep '\.1\.fastq'|grep $barcode);
		do zcat $datadir/$j >> ${setname}.r1.fastq
	done
	for j in $(ls $datadir|grep '\.barcode_1\.fastq'|grep $barcode);
		do zcat $datadir/$j >> ${setname}.bc1.fastq
	done
	for j in $(ls $datadir|grep '\.2\.fastq'|grep $barcode);
		do zcat $datadir/$j >> ${setname}.r2.fastq
	done
	gzip ${setname}.r1.fastq ${setname}.bc1.fastq ${setname}.r2.fastq
	echo -e "$i\t$barcode" > ${setname}.barcodekey
	echo $setname >> samplesets.txt
done

