//Xander Nuttle
//finalseqs_to_mipcounts.c
//Call: ./finalseqs_to_mipcounts gzipped_finalseqs_file
//
//Generates a mipcounts file for automated copy number genotyping from a finalized set of filtered MIP sequences.
//Sample name (first column of input file) should include "8330" or "2069" to specify whether the sample corresponds
//to the GM08330 background or the MGH2069 background. Otherwise, output will only include counts from 'B' type MIPs
//which target at least one SNV in each background and omit 'E' type and 'M' type MIPs which are also informative
//for copy number genotyping.
//
//The idea is to take counts from the top two most abundant sequences corresponding to each MIP target. Even though
//we do not know haplotypes, if there is a deletion or duplication, the most abundant sequences across the CNV interval
//should correspond to the same haplotype (since the other haplotype is deleted or present at one fewer copy).

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<zlib.h>
#define NLEN 200 //maximum length of names (sample, MIP, contig, mipseqs file) and of mapping coordinate converted to a string
#define LLEN 1500 //maximum length of single line of text in input finalseqs file
#define SLEN 500 //maximum length of each sequence array (and corresponding quality array)

//set up structure to store data for each sequence at a guide target
struct mipseq
{
	char mip[NLEN+1];
	char miptype;
	char crispr[NLEN+1];
	char contig[NLEN+1];
	char maploc[NLEN+1];
	char seq[SLEN+1];
	char qual[SLEN+1];
	long tagcount;
	char tagfreq[SLEN+1];
};

FILE*init_output(FILE*out,char*basename);
int countseqs(gzFile*fseqs,char*lyne,long*numseqs);
void get_seqs(gzFile*fseqs,char*lyne,long numseqs,struct mipseq*sequences);
int informative(struct mipseq*sequences,char*samp);
void print_seqs(FILE*out,char*samp,struct mipseq*sequences,long numseqs);

int main(int argc,char*argv[])
{
	//get sample name from command line
	char sample[NLEN+1];
	strncpy(sample,*(argv+1),NLEN-10); //new extension has 10 characters
	sample[strchr(sample,'.')-sample]='\0';

	//set up output file
	FILE*mipcounts=init_output(mipcounts,sample);

	//read in data for finalized MIP sequences, processing them in groups based on their associated MIP target
	char line[LLEN+1];
	long nseqs;
	struct mipseq*seqs;
	gzFile*finalseqs=gzopen(*(argv+1),"r");
	gzgets(finalseqs,line,LLEN-1); //process header line
	while(countseqs(finalseqs,line,&nseqs))
	{
		//read sequences into array of mipseq structures
		seqs=(struct mipseq*)malloc(nseqs*sizeof(struct mipseq));
		get_seqs(finalseqs,line,nseqs,seqs);

		//if MIP is informative for copy number genotyping, print data from two most abundant sequences (based on tag counts) to output
		if(informative(seqs,sample))
			print_seqs(mipcounts,sample,seqs,nseqs);
		free(seqs);
	}

	//clean up and exit
	fclose(mipcounts);
	gzclose(finalseqs);
	return 0;
}

FILE*init_output(FILE*out,char*basename)
{
	char outname[NLEN+1];
	sprintf(outname,"%s%s",basename,".mipcounts");
	out=fopen(outname,"w");
	fprintf(out,"Sample\tContig\tCoordinate\tMip_Type\tHaplotype_1_Count\tHaplotype_2_Count\n");
	return out;
}

int countseqs(gzFile*fseqs,char*lyne,long*numseqs)
{
	char mipone[NLEN+1],newmip[NLEN+1];
	long offset=0;
	(*numseqs)=0;
	while(gzgets(fseqs,lyne,LLEN-1))
	{
		offset+=strlen(lyne);
		if((*numseqs)==0)
		{
			sscanf(lyne,"%*s %s",mipone);
			(*numseqs)++;
			continue;
		}
		sscanf(lyne,"%*s %s",newmip);
		if(strncmp(newmip,mipone,NLEN)==0)
			(*numseqs)++;
		else
			break;
	}
	gzseek(fseqs,(-1*offset),SEEK_CUR);
	return ((*numseqs)>0);
}

void get_seqs(gzFile*fseqs,char*lyne,long numseqs,struct mipseq*sequences)
{
	long s;
	for(s=0;s<numseqs;s++)
	{
		gzgets(fseqs,lyne,LLEN-1);
		sscanf(lyne,"%*s %s %c %s %s %s %s %s %ld %s",sequences[s].mip,&(sequences[s].miptype),sequences[s].crispr,sequences[s].contig,sequences[s].maploc,sequences[s].seq,sequences[s].qual,&(sequences[s].tagcount),sequences[s].tagfreq);
	}
	return;
}

int informative(struct mipseq*sequences,char*samp)
{
	if((sequences[0].miptype=='B')||((strstr(samp,"8330"))&&(sequences[0].miptype=='E'))||((strstr(samp,"2069"))&&(sequences[0].miptype=='M')))
		return 1;
	else
		return 0;
}

void print_seqs(FILE*out,char*samp,struct mipseq*sequences,long numseqs)
{
	long counttwo=0;
	if(numseqs>1)
		counttwo=sequences[1].tagcount;
	fprintf(out,"%s\t%s\t%s\t%c\t%ld\t%ld\n",samp,sequences[0].contig,sequences[0].maploc,sequences[0].miptype,sequences[0].tagcount,counttwo);
	return;	
}

