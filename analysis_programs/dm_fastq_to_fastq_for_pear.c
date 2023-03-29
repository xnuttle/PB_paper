//Xander Nuttle
//dm_fastq_to_fastq_for_pear.c
//Call: ./dm_fastq_to_fastq_for_pear reads1.fastq.gz reads2.fastq.gz reads3.fastq.gz reads4.fastq.gz (int)trimmed_read_length (int)molecular_tag_length (long)max_num_reads_per_output_file barcodekey_file
//
//Takes four demultiplexed gzipped fastq files together containing data for a single sample from a sequencing run (a file with all first reads,
//a file with all index1 reads, a file with all index2 reads, and a file with all second reads) and generates gzipped fastq files as output.
//Only adds reads having a barcode sequence combination perfectly matching the sample's known barcode sequence combination to the output fastq files.
//Generates 1 or more sets of gzipped fastq files for the sample. Assumes MIPs have been molecularly tagged at one location just internal to the
//extension arm to allow for counting of individual capture events.
//
//This program converts demultiplexed gzipped fastq files (containing data for a single sample from a sequencing run) into smaller gzipped fastq files, which serve
//as input to PEAR, a program to merge overlapping reads. There are four demultiplexed gzipped fastq files per sample - the first file contains the first
//sequence reads from read pairs, the second file contains index1 reads that associate the read pairs with that particular sample, the third file contains index2
//reads that associate the read pairs with that particular sample, and the fourth file contains the second sequence reads from read pairs. Sequence data for each
//read pair is included in corresponding locations in each of these files. For example, lines 1-4 of the first file contains data for the first read of the
//first read pair, lines 1-4 of the second file contains data for the index1 read corresponding to the first read pair, lines 1-4 of the third file contains
//data for the index2 read corresponding to the first read pair, and lines 1-4 of the fourth file contains data for the second read of the first read pair.
//
//This program performs several functions:
//  -converts a set of four demultiplexed large gzipped fastq files into multiple gzipped fastq files, retaining only read pairs with index reads perfectly
//   matching the sample's barcode sequences, with each read in the output gzipped fastq files tagged with the molecular tag
//   (each gzipped fastq output file will contain a maximum of a user-specified number of reads [250,000 or fewer is recommended])
//  -trims reads, retaining only a user-specified number of bases of each sequence read from the 5' ends of reads (eliminating lower quality bases at the 3'
//   ends of reads)
//
//This program does not filter reads based on quality except that read pairs with index reads not perfectly matching the sample's known barcode sequences are
//not retained. This approach allows downstream analysis to handle low-quality reads rather than throwing out reads based on not passing a stringent quality filter.
//It also prevents only one read from a read pair being retained.
//
//This program deals with a space in sequence names and accomodates dual-index barcoding. It prints reads 1 and reads 2 to separate output files, thus preparing reads
//for merging with the program PEAR (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3933873/).

#include<stdio.h>
#include<zlib.h>
#include<string.h>
#include<stdlib.h>
#define LEN 101 //maximum length of sample names and barcode sequences + 1

int main(int argc,char*argv[])
{
	//read in sample name and barcode from barcodekey file
	char sample[LEN],barcode[LEN];
	FILE*barcodekey=fopen(*(argv+8),"r");
	fscanf(barcodekey,"%s %s",sample,barcode);	
	int bc_length=strlen(barcode);
	fclose(barcodekey);
	
	//setup input files
	gzFile*in1,*in2,*in3,*in4;

	//setup output files
	char outname1[LEN],outname2[LEN];
  gzFile*outfiles[2];
	sprintf(outname1,"%s_FS1_F1.fastq.gz\0",sample);
	sprintf(outname2,"%s_FS1_R1.fastq.gz\0",sample);
	outfiles[0]=gzopen(outname1,"w");
	outfiles[1]=gzopen(outname2,"w");

	//setup variables: specify read length, index read length, and desired number of sequence reads per fastq file
	long trimmed_read_length; //(76 bp is the minimal value to cover all targeted bases (112 bp) + both hybridization arms (20 bp each)), 142 bp recommended for sequence analysis
	long tag_length;
  long reads_per_fastq=strtol(*(argv+7),NULL,10); //250,000 recommended
	int output_file_num=1;

	//read in trimmed read length value from the command line
	trimmed_read_length=strtol(*(argv+5),NULL,10);

	//read in molecular tag length value from the command line
	tag_length=strtol(*(argv+6),NULL,10);

	//setup variables
	in1=gzopen(*(argv+1),"r");
  in2=gzopen(*(argv+2),"r");
  in3=gzopen(*(argv+3),"r");
	in4=gzopen(*(argv+4),"r");
	char line[501],line2[501];
  line[500]='\0';
	line2[500]='\0';
	long reads_output=0;
	int i,k;
	char index_sequence[bc_length+1];
	index_sequence[bc_length]='\0';
	char tag_sequence[tag_length+1];
	tag_sequence[tag_length]='\0';

	//read in large gzipped fastq files line by line, trim sequences, and print output
	while(gzgets(in2,line,500))
	{
		if(reads_output>=reads_per_fastq)
		{
			gzclose(outfiles[0]);
			gzclose(outfiles[1]);
			output_file_num++;
			sprintf(outname1,"%s_FS1_F%d.fastq.gz\0",sample,output_file_num);
			sprintf(outname2,"%s_FS1_R%d.fastq.gz\0",sample,output_file_num);
			outfiles[0]=gzopen(outname1,"w");
			outfiles[1]=gzopen(outname2,"w");
			reads_output=0;
		}
		
		//process data for index1 sequence
		gzgets(in2,line,500);
		strncpy(index_sequence,line,bc_length/2); //index quality will not later explicity be taken into account, but it will in the sense that only index sequences with perfect matches to known barcodes will be used
		index_sequence[bc_length/2]='\0';
		gzgets(in2,line,500);
		gzgets(in2,line,500);		

		//process data for index2 sequence
		gzgets(in3,line,500);
		gzgets(in3,line,500);
		strncat(index_sequence,line,bc_length/2); //index quality will not later explicity be taken into account, but it will in the sense that only index sequences with perfect matches to known barcodes will be used
		gzgets(in3,line,500);
		gzgets(in3,line,500);

		//ensure barcode sequence perfectly matches a known barcode reverse complement sequence, and if so, identify which individual the read pair corresponds to
		int indiv=-1;
		if(strncmp(index_sequence,barcode,bc_length)==0)
			indiv=1;

		//process data for second sequence
    gzgets(in4,line,500);
    gzgets(in4,line2,500);
		strncpy(tag_sequence,line2,tag_length);
		if((indiv!=-1)&&(!(strchr(tag_sequence,'N')))&&(strstr(tag_sequence,"AAAAA")==NULL)&&(strstr(tag_sequence,"CCCCC")==NULL)&&(strstr(tag_sequence,"GGGGG")==NULL)&&(strstr(tag_sequence,"TTTTT")==NULL))
    {
    	line[strlen(line)-5]='\0'; //remove the newline from the string "line", as well as the '#0/3'
			if(strchr(line,' ')!=NULL)
			{
				line[strchr(line,' ')-line]='\0'; //if the string "line" contains a space, remove it and everything afterwards from the string "line"
			}
			strncat(line,"/2 MI:Z:$",9);
			strncat(line,tag_sequence,tag_length); //add molecular tag information to sequence name
			strncat(line,"\n",1);
			gzputs(outfiles[1],line);
			strncpy(line,line2+tag_length,trimmed_read_length);
			line[trimmed_read_length]='\n';
      line[trimmed_read_length+1]='\0';
      gzputs(outfiles[1],line);
    }
		gzgets(in4,line,500);
    gzgets(in4,line2,500);
		if((indiv!=-1)&&(!(strchr(tag_sequence,'N')))&&(strstr(tag_sequence,"AAAAA")==NULL)&&(strstr(tag_sequence,"CCCCC")==NULL)&&(strstr(tag_sequence,"GGGGG")==NULL)&&(strstr(tag_sequence,"TTTTT")==NULL))
    {
      gzputs(outfiles[1],line);
      strncpy(line,line2+tag_length,trimmed_read_length);
      line[trimmed_read_length]='\n';
      line[trimmed_read_length+1]='\0';
      gzputs(outfiles[1],line);
			reads_output++;
    }

		//process data for first sequence
		for(i=0;i<4;i++)
    {
      gzgets(in1,line,500);
      if((indiv!=-1)&&(!(strchr(tag_sequence,'N')))&&(strstr(tag_sequence,"AAAAA")==NULL)&&(strstr(tag_sequence,"CCCCC")==NULL)&&(strstr(tag_sequence,"GGGGG")==NULL)&&(strstr(tag_sequence,"TTTTT")==NULL))
      {
        if(i==0)
				{
					line[strlen(line)-5]='\0'; //remove the newline from the string "line", as well as the '#0/1'
      		if(strchr(line,' ')!=NULL)
					{
						line[strchr(line,' ')-line]='\0'; //if the string "line" contains a space, remove it and everything afterwards from the string "line"
					}
					strncat(line,"/1 MI:Z:$",9);
      		strncat(line,tag_sequence,tag_length); //add molecular tag information to sequence name
      		strncat(line,"\n",1);
				}
				if(i%2)
      	{
        	line[trimmed_read_length]='\n';
        	line[trimmed_read_length+1]='\0';
      	}
				gzputs(outfiles[0],line);
      }
    }
	}

	//clean up and exit
	gzclose(in1);
  gzclose(in2);
  gzclose(in3);
  gzclose(in4);
	gzclose(outfiles[0]);
	gzclose(outfiles[1]);
  return 0;
}

