//Xander Nuttle
//finalize_mipseqs.c
//Call: ./finalize_mipseqs gzipped_seqcounts_file (long)depth_cutoff (double)allele_fraction_cutoff
//
//Filters sequences in a gzipped seqcounts file based on molecular tag count depth and allele balance, outputting
//a finalized set of filtered sequences together with associated data for each sequence (quality, tag count, etc.)
//
//depth_cutoff = minimum molecular tag count; all sequences with fewer corresponding tag counts will be discarded
//allele_fraction_cutoff = minimum allele fraction; all sequences with a lower allele fraction (after discarding sequences below the depth cutoff) will be discarded

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<zlib.h>
#define NLEN 200 //maximum length of names (sample, MIP, contig, mipseqs file) and of mapping coordinate converted to a string
#define SLEN 500 //maximum length of each sequence array (and corresponding quality array)
#define LLEN 1500 //maximum length of single line of text in input seqcounts file

//set up structure to store data for each sequence at a guide target
struct mipseq
{
	char mip[NLEN+1];
	char miptype[NLEN+1];
	char crispr[NLEN+1];
	char contig[NLEN+1];
  char maploc[NLEN+1];
  char seq[SLEN+1];
  char qual[SLEN+1];
  long tagcount;
  double tagfreq;
};

gzFile* init_output(gzFile*fseqs,char*basename,long dp,double af,char*afstr);
int countseqs(gzFile*scounts,char*lyne,long*numseqs);
void get_seqs(gzFile*scounts,char*lyne,long numseqs,struct mipseq*sequences);
int compfun(const void*p1,const void*p2);
void filter_dp(struct mipseq*sequences,long numseqs,long dp);
long counttags(struct mipseq*sequences,long numseqs);
void filter_af(struct mipseq*sequences,long numseqs,double af,long count);
void print_seqs(gzFile*fseqs,char*samp,struct mipseq*sequences,long numseqs);

int main(int argc,char*argv[])
{
	//get sample name, depth cutoff, and allele fraction cutoff from command line
	char sample[NLEN+1];
	strncpy(sample,*(argv+1),NLEN-30); //leave up to 30 characters for new extension
	sample[strchr(sample,'.')-sample]='\0';
	long mindp=strtol(*(argv+2),NULL,10);
	double minaf=strtod(*(argv+3),NULL);

	//set up output file
	gzFile*finalseqs=init_output(finalseqs,sample,mindp,minaf,*(argv+3));

	//read in mipseqs and associated tag counts from seqcounts file, processing them in groups based on their associated MIP target	
	char line[LLEN+1];
	long nseqs,ntags;
	struct mipseq*seqs;
	gzFile*seqcounts=gzopen(*(argv+1),"r");	
	gzgets(seqcounts,line,LLEN-1); //process header line
	while(countseqs(seqcounts,line,&nseqs))
	{
		//read sequences into array of mipseq structures
		seqs=(struct mipseq*)malloc(nseqs*sizeof(struct mipseq));
		get_seqs(seqcounts,line,nseqs,seqs);

		//sort array of mipseq structures based on tagcount numbers
		qsort(seqs,nseqs,sizeof(struct mipseq),compfun);

		//filter out sequences based on molecular tag count depth and allele balance
		filter_dp(seqs,nseqs,mindp);		
		ntags=counttags(seqs,nseqs);
		filter_af(seqs,nseqs,minaf,ntags);

		//print remaining sequences to output file and free memory used to store data for current set of sequences
		print_seqs(finalseqs,sample,seqs,nseqs);
		free(seqs);
	}

	//clean up and exit
	gzclose(finalseqs);
	gzclose(seqcounts);
	return 0;
}

gzFile* init_output(gzFile*fseqs,char*basename,long dp,double af,char*afstr)
{
	char outname[NLEN+1];
	sprintf(outname,"%s%s%ld%s%.*lf%s",basename,".dp",dp,".af",strlen(afstr)-(strchr(afstr,'.')+1-afstr),af,".finalseqs.gz\0");
	fseqs=gzopen(outname,"w");
	gzprintf(fseqs,"Sample\tMIP\tType\tCRISPR\tContig\tCoordinate\tSequence\tQuality\tTagCount\tAlleleFraction\n");
	return fseqs;
}

int countseqs(gzFile*scounts,char*lyne,long*numseqs)
{
	char mipone[NLEN+1],newmip[NLEN+1];
	long offset=0;
	(*numseqs)=0;
	while(gzgets(scounts,lyne,LLEN-1))
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
	gzseek(scounts,(-1*offset),SEEK_CUR);
	return ((*numseqs)>0);
}

void get_seqs(gzFile*scounts,char*lyne,long numseqs,struct mipseq*sequences)
{
	long s;
	for(s=0;s<numseqs;s++)
	{
		gzgets(scounts,lyne,LLEN-1);
		sscanf(lyne,"%*s %s %s %s %s %s %s %s %ld",sequences[s].mip,sequences[s].miptype,sequences[s].crispr,sequences[s].contig,sequences[s].maploc,sequences[s].seq,sequences[s].qual,&(sequences[s].tagcount));
		sequences[s].tagfreq=0.0;
	}
	return;
}

int compfun(const void*p1,const void*p2)
{
	const struct mipseq*seq1=p1;
	const struct mipseq*seq2=p2;
	if(seq1->tagcount<seq2->tagcount)
    return 1;
  else if(seq2->tagcount<seq1->tagcount)
    return -1;
	else
		return 0;
}

void filter_dp(struct mipseq*sequences,long numseqs,long dp)
{
	long s;
	for(s=0;s<numseqs;s++)
	{
		if(sequences[s].tagcount<dp)
			sequences[s].tagcount=0;
	}
	return;
}

long counttags(struct mipseq*sequences,long numseqs)
{
	long s,count=0;
	for(s=0;s<numseqs;s++)
		count+=sequences[s].tagcount;
	return count;
}

void filter_af(struct mipseq*sequences,long numseqs,double af,long count)
{
	long s;
	for(s=0;s<numseqs;s++)
	{
		sequences[s].tagfreq=(double)sequences[s].tagcount/(double)count;
		if(sequences[s].tagfreq<af)
			sequences[s].tagcount=0;
	}
	return;
}

void print_seqs(gzFile*fseqs,char*samp,struct mipseq*sequences,long numseqs)
{
	long s;
	for(s=0;s<numseqs;s++)
	{
		if(sequences[s].tagcount>0)
			gzprintf(fseqs,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%ld\t%lf\n",samp,sequences[s].mip,sequences[s].miptype,sequences[s].crispr,sequences[s].contig,sequences[s].maploc,sequences[s].seq,sequences[s].qual,sequences[s].tagcount,sequences[s].tagfreq);
		else
			break;
	}
	return;
}

