//Xander Nuttle
//get_guides_tags.c
//Call: ./get_guides_tags gzipped_finalseqs_file guides_table
//
//This program analyzes finalized MIP sequence data to extract molecular tag information associated with each integrated
//guide RNA construct (assumes a molecularly-tagged indel guide library was used along with a MIP to capture guide sequences and
//associated molecular tags). This information can then be used for rarefaction analysis to assess diversity of integrated
//guide constructs. The guides table file is a headerless, tab-delimited text file with the first column containing names of all tagged guides
//and the second column containing guide sequences.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<zlib.h>
#include<ctype.h>
#define NLEN 200 //maximum length of names (sample, MIP, contig, mipseqs file) and of mapping coordinate converted to a string
#define SLEN 500 //maximum length of each sequence array (and corresponding quality array)
#define LLEN 1500 //maximum length of single line of text in input finalseqs file
#define TLEN 10 //length of molecular tag for guide RNA construct
#define RLEN 86 //number of reference sequence bases (excluding the guide sequence) before guide molecular tag

//set up structure to store data for each target sequence
struct mipseq
{
	char mip[NLEN+1];
	char miptype[NLEN+1];
	char crispr[NLEN+1];
	char contig[NLEN+1];
	char maploc[NLEN+1];
	char seq[SLEN+1];
	char qual[SLEN+1];
	char tagcount[NLEN+1];
	char tagfreq[NLEN+1];
	char gtag[TLEN+1];
};

//set up structure to store data for each guide
struct guide
{
	char name[NLEN+1];
	char seq[SLEN+1];
};

long count_guides(FILE*glist);
void get_guides(FILE*glist,struct guide*gyds);
int getinput(gzFile*fseqs,struct mipseq*iseq);
void process(char*samp,struct mipseq*iseq,struct guide*gyds,long numguides);
void parse_seq(char*seq,char*newseq,long*rbases);
void parse_match(char*orig,char*new,long*rb,long*o,long*n,long*r);
void parse_ins(char*orig,char*new,long*rb,long*o,long*n,long*r);
void parse_del(char*orig,char*new,long*rb,long*o,long*n,long*r);
void parse_sub(char*orig,char*new,long*rb,long*o,long*n,long*r);
void pad_tag(char*tag);

int main(int argc,char*argv[])
{
	//get sample name from command line
	char sample[NLEN+1];
	strncpy(sample,*(argv+1),NLEN);
	sample[strchr(sample,'.')-sample]='\0';

	//determine the number of guides and allocate memory to store guide information
	FILE*guidelist=fopen(*(argv+2),"r");
	long nguides=count_guides(guidelist);
	struct guide*guides;
	guides=(struct guide*)malloc(nguides*sizeof(struct guide));
	
	//read in information on guides
	get_guides(guidelist,guides);

	//read in final called MIP sequences and process them one by one
	gzFile*finalseqs=gzopen(*(argv+1),"r");
	struct mipseq inseq;
	getinput(finalseqs,&inseq); //process header line
	while(getinput(finalseqs,&inseq))
		process(sample,&inseq,guides,nguides);

	//clean up and exit
	gzclose(finalseqs);
	return 0;
}

long count_guides(FILE*glist)
{
	long numguides=0;
  char gname[NLEN+1];
  fpos_t pos;
  fgetpos(glist,&pos);
  while(fscanf(glist,"%s %*s",gname)==1)
    numguides++;
  fsetpos(glist,&pos);
  return numguides;
}

void get_guides(FILE*glist,struct guide*gyds)
{
	long g=0;
	while(fscanf(glist,"%s %s",gyds[g].name,gyds[g].seq)==2)
		g++;
	return;
}

int getinput(gzFile*fseqs,struct mipseq*iseq)
{
	char line[LLEN+1];
	int scanned=0;
	if(gzgets(fseqs,line,LLEN))
		scanned+=sscanf(line,"%*s %s %s %s %s %s %s %s %s %s",iseq->mip,iseq->miptype,iseq->crispr,iseq->contig,&(iseq->maploc),iseq->seq,iseq->qual,iseq->tagcount,iseq->tagfreq);
	return (scanned==9);
}

void process(char*samp,struct mipseq*iseq,struct guide*gyds,long numguides)
{
	char newseq[SLEN+1];
	long refbases[SLEN+1];
	long b,g,tag_reached=0;
	char base;
	char*pretag;
	if((strstr(iseq->mip,"_guide_"))&&(strstr(iseq->mip,"_MIP_0002")))
	{
		for(g=0;g<numguides;g++)
		{
			if(strncmp(iseq->contig,gyds[g].name,NLEN)==0)
				break;
		}
		parse_seq(iseq->seq,newseq,refbases);
		for(b=0;b<strlen(newseq);b++)
		{
			base=toupper(newseq[b]);
			newseq[b]=base;
		}
		pretag=strstr(newseq,"TTTTTT");
		if(pretag!=NULL)
		{
			b=(long)pretag-(long)(&(newseq[0]))+6;
			strncpy(iseq->gtag,newseq+b,TLEN);
			iseq->gtag[TLEN]='\0';
			if(strlen(iseq->gtag)<TLEN)
				pad_tag(iseq->gtag);
			if(iseq->gtag[0]!='N')
				printf("%s\t%s\t%s\t%s\t%s\n",samp,iseq->contig,iseq->gtag,iseq->tagcount,iseq->tagfreq);
		}
	}
	return;
}

void parse_seq(char*seq,char*newseq,long*rbases)
{
	long i=0,j=0,k=0;
	while(i<strlen(seq))
	{
		switch(seq[i])
		{
			case '=': i++; parse_match(seq,newseq,rbases,&i,&j,&k); break;
			case '+': i++; parse_ins(seq,newseq,rbases,&i,&j,&k); break;
			case '-': i++; parse_del(seq,newseq,rbases,&i,&j,&k); break;
			case '*': i++; parse_sub(seq,newseq,rbases,&i,&j,&k); break;
		}
	}
	newseq[j]='\0';
	return;
}

void parse_match(char*orig,char*new,long*rb,long*o,long*n,long*r)
{
	while(isalpha(orig[*o]))
	{
		new[*n]=orig[*o];
		rb[*n]=(*r);
		(*o)++;
		(*n)++;
		(*r)++;
	}
	return;
}

void parse_ins(char*orig,char*new,long*rb,long*o,long*n,long*r)
{
	while(isalpha(orig[*o]))
	{
		new[*n]=orig[*o];
		rb[*n]=(*r);
		(*o)++;
		(*n)++;
	}
	return;
}

void parse_del(char*orig,char*new,long*rb,long*o,long*n,long*r)
{
	while(isalpha(orig[*o]))
	{
		(*o)++;
		(*r)++;
	}
	return;
}

void parse_sub(char*orig,char*new,long*rb,long*o,long*n,long*r)
{
	(*o)++;
	new[*n]=orig[*o];
	rb[*n]=(*r);
	(*o)++;
	(*n)++;
	(*r)++;
	return;
}

void pad_tag(char*tag)
{
	long b;
	for(b=0;b<TLEN;b++)
	{
		if(tag[b]=='\0')
			tag[b]='N';
	}
	tag[b]='\0';
	return;
}

