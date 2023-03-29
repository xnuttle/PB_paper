//Xander Nuttle
//count_mipseqs.c
//Call: ./count_mipseqs text_file_with_names_of_gzipped_mipseqs_files(sample.seqsfiles) miptargets_file
//
//This program analyses a set of gzipped mipseqs files for a sample and outputs all distinct sequences at each MIP target
//along with all corresponding information, including the number of different molecular tags associated with that sequence.
//This information can then be used to determine which sequences should be deemed present at each MIP target site based off
//tag count and allele fraction filtering with the program finalize_mipseqs.c

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<zlib.h>
#define NLEN 200 //maximum length of names (sample, contig, mipseqs file) and of mapping coordinate converted to a string
#define SLEN 500 //maximum length of each sequence array (and corresponding quality array)
#define TLEN 8 //length of molecular tag sequences
#define LLEN 1500 //maximum length of single line of text in input mipseqs file

//set up structure to store MIP target information, including all sequences assigned to the MIP and associated tag counts
struct miptarg
{
	char mip[NLEN+1];
	char miptype[NLEN+1];
	char crispr[NLEN+1];
	struct mipseq*seqs;
};

//set up structure to store data for each distinct sequence at a MIP target
struct mipseq
{
	char contig[NLEN+1];
	char maploc[NLEN+1];
	char seq[SLEN+1];
	double qual[SLEN+1];
	long seqcount;
	long tagcount;
	struct moltag*tags;
	struct mipseq*next;
};

//set up structure to store data for each molecular tag associated with a guide target sequence
struct moltag
{
	char tag[TLEN+1];
	struct moltag*next;
};

//set up structure to store data for each input sequence
struct input
{
	char mip[NLEN+1];
	char miptype[NLEN+1];
	char crispr[NLEN+1];
	char contig[NLEN+1];
	char maploc[NLEN+1];
	char seq[SLEN+1];
	char qual[SLEN+1];
	char tag[TLEN+1];
};

long count_targs(FILE*mtargs);
void init_targs(struct miptarg*targs,FILE*mtargs);
int getinput(gzFile*mseqs,struct input*iseq);
void process(struct input*iseq,struct miptarg*targs,long numtargs);
long findtarg(char*myp,struct miptarg*targets,long ntargets);
int isnew(struct input*seqin,struct miptarg*targets,long index);
void addseq(struct input*seqin,struct miptarg*targets,long index);
void init_seqs(struct input*seqi,struct miptarg*mtargets,long indx);
void init_qual(double*curqual,char*newqual);
void update(struct input*seqin,struct miptarg*targets,long index);
void udqual(double*curqual,char*newqual);
void udmapping(char*curchr,char*curcoord,char*newchr,char*newcoord);
int newtag(char*intag,struct moltag*taglist);
void udtags(char*intag,struct moltag*taglist);
gzFile* init_output(gzFile*scounts,char*basename);
void print_data(gzFile*scounts,char*samp,struct miptarg*targs,long numtargs);
char*avgqual(double*curqual,long count,char*curseq,char*newqual);
void freeseqs(struct miptarg*targs,long numtargs);
void freetags(struct moltag*taglist);

int main(int argc,char*argv[])
{
	//get sample name
	char sample[NLEN+1];
	strncpy(sample,*(argv+1),NLEN-13); //13 characters needed for new extension
	sample[strchr(sample,'.')-sample]='\0';

	//determine the number of MIP targets and allocate memory to store MIP target information
	FILE*miptargs=fopen(*(argv+2),"r");
	long ntargs=count_targs(miptargs);
	struct miptarg*mtargs;
	mtargs=(struct miptarg*)malloc(ntargs*sizeof(struct miptarg));	

	//read in MIP names and initialize MIP target data
	init_targs(mtargs,miptargs);

	//read in names of gzipped mipseqs files to process and process gzipped mipseqs files one by one
	FILE*filelist=fopen(*(argv+1),"r");
	gzFile*mipseqs;
	char seqfile[NLEN+1];
	struct input inseq; 
	while(fscanf(filelist,"%s",seqfile)==1)
	{
		mipseqs=gzopen(seqfile,"r");
		getinput(mipseqs,&inseq); //process header line		
		while(getinput(mipseqs,&inseq))
			process(&inseq,mtargs,ntargs);	
		gzclose(mipseqs);
	}

	//set up output file and print data for each guide target
	gzFile*seqcounts=init_output(seqcounts,sample);
	print_data(seqcounts,sample,mtargs,ntargs);

	//clean up and exit
	freeseqs(mtargs,ntargs);
	free(mtargs);
	gzclose(seqcounts);
	fclose(filelist);
	fclose(miptargs);
	return 0;
}

long count_targs(FILE*mtargs)
{
	long numtargs=0;
	char mipname[NLEN+1];
	fpos_t pos;
	fscanf(mtargs,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
	fgetpos(mtargs,&pos);
	while(fscanf(mtargs,"%s %*s %*s %*s %*s %*s %*s %*s %*s %*s",mipname)==1)
	{
		numtargs++;
	}
	fsetpos(mtargs,&pos);
	return numtargs;
}

void init_targs(struct miptarg*targs,FILE*mtargs)
{
	long m=0;
	while(fscanf(mtargs,"%s %*s %*s %*s %*s %s %s %*s %*s %*s",targs[m].mip,targs[m].miptype,targs[m].crispr)==3)
	{
		targs[m].seqs=NULL;
		m++;
	}
	return;
}

int getinput(gzFile*mseqs,struct input*iseq)
{
	char line[LLEN+1];
	int scanned=0;
	if(gzgets(mseqs,line,LLEN-1))
	{
		scanned+=sscanf(line,"%*s %s %s %s %s %s %s %s %s",iseq->mip,iseq->miptype,iseq->crispr,iseq->contig,iseq->maploc,iseq->seq,iseq->qual,iseq->tag);
	}
	return (scanned==8);
}

void process(struct input*iseq,struct miptarg*targs,long numtargs)
{
	long m=findtarg(iseq->mip,targs,numtargs);
	if(isnew(iseq,targs,m))
	{
		addseq(iseq,targs,m);
	}
	else
	{
		update(iseq,targs,m);
	}
	return;
}

long findtarg(char*myp,struct miptarg*targets,long ntargets)
{
	long mnum;
	for(mnum=0;mnum<ntargets;mnum++)
	{
		if(strncmp(targets[mnum].mip,myp,NLEN)==0)
			break;
	}
	if(mnum==ntargets)
		mnum=-1;
	return mnum;
}

int isnew(struct input*seqin,struct miptarg*targets,long index)
{
	struct mipseq*current=targets[index].seqs;
	while(current!=NULL)
	{
		if(strncmp(current->seq,seqin->seq,SLEN)==0)
			return 0;
		current=current->next;
	}
	return 1;
}

void addseq(struct input*seqin,struct miptarg*targets,long index)
{
	struct mipseq*current=targets[index].seqs;
	if(current==NULL)
		init_seqs(seqin,targets,index);
	else
	{
		while(current->next!=NULL)
			current=current->next;
		current->next=(struct mipseq*)malloc(sizeof(struct mipseq));
		strncpy(current->next->contig,seqin->contig,NLEN);
		strncpy(current->next->maploc,seqin->maploc,NLEN);
		strncpy(current->next->seq,seqin->seq,SLEN);
		init_qual(current->next->qual,seqin->qual);
		current->next->seqcount=1;
		current->next->tagcount=1;
		current->next->tags=(struct moltag*)malloc(sizeof(struct moltag));
		strncpy(current->next->tags->tag,seqin->tag,TLEN);
		current->next->tags->next=NULL;
		current->next->next=NULL;
	}
	return;
}

void init_seqs(struct input*seqi,struct miptarg*mtargets,long indx)
{
	mtargets[indx].seqs=(struct mipseq*)malloc(sizeof(struct mipseq));
	struct mipseq*cur=mtargets[indx].seqs;
	strncpy(cur->contig,seqi->contig,NLEN);
	strncpy(cur->maploc,seqi->maploc,NLEN);
	strncpy(cur->seq,seqi->seq,SLEN);
	init_qual(cur->qual,seqi->qual);
	cur->seqcount=1;
	cur->tagcount=1;
	cur->tags=(struct moltag*)malloc(sizeof(struct moltag));
	strncpy(cur->tags->tag,seqi->tag,TLEN);
	cur->tags->next=NULL;
	cur->next=NULL;
	return;
}

void init_qual(double*curqual,char*newqual)
{
	long q;
	for(q=0;q<strlen(newqual);q++)
	{
		curqual[q]=(double)newqual[q];
	}
	curqual[q]=0.0;
	return;
}

void update(struct input*seqin,struct miptarg*targets,long index)
{
	struct mipseq*current=targets[index].seqs;
	while(strncmp(current->seq,seqin->seq,SLEN)!=0)
		current=current->next;
	udqual(current->qual,seqin->qual);
	udmapping(current->contig,current->maploc,seqin->contig,seqin->maploc);
	current->seqcount++;
	if(newtag(seqin->tag,current->tags))
	{
		udtags(seqin->tag,current->tags);
		current->tagcount++;
	}
	return;
}

void udqual(double*curqual,char*newqual)
{
	long q;
	for(q=0;q<strlen(newqual);q++)
	{
		curqual[q]+=(double)newqual[q];
	}
	return;
}

void udmapping(char*curchr,char*curcoord,char*newchr,char*newcoord)
{
	if(strstr(curchr,newchr)==NULL)
	{
		strncat(curchr,"/",1);
		strncat(curchr,newchr,strlen(newchr));
	}
	if(strstr(curcoord,newcoord)==NULL)
	{
		strncat(curcoord,"/",1);
		strncat(curcoord,newcoord,strlen(newcoord));
	}
	return;
}

int newtag(char*intag,struct moltag*taglist)
{
	struct moltag*cur=taglist;
	while(cur!=NULL)
  {
    if(strncmp(cur->tag,intag,TLEN)==0)
      return 0;
    cur=cur->next;
  }
  return 1;
}

void udtags(char*intag,struct moltag*taglist)
{
	struct moltag*cur=taglist;
	while(cur->next!=NULL)
		cur=cur->next;
	cur->next=(struct moltag*)malloc(sizeof(struct moltag));
	strncpy(cur->next->tag,intag,TLEN);
	cur->next->next=NULL;
	return;
}

gzFile* init_output(gzFile*scounts,char*basename)
{
	char outname[NLEN+1];
	sprintf(outname,"%s%s",basename,".seqcounts.gz\0");
	scounts=gzopen(outname,"w");
	gzprintf(scounts,"Sample\tMIP\tType\tCRISPR\tContig\tCoordinate\tSequence\tQuality\tTagCount\n");
	return scounts;
}

void print_data(gzFile*scounts,char*samp,struct miptarg*targs,long numtargs)
{
	long m;
	char finalqual[SLEN+1];
	struct mipseq*current;
	for(m=0;m<numtargs;m++)
	{
		current=targs[m].seqs;
		while(current!=NULL)
		{
			gzprintf(scounts,"%s\t%s\t%s\t%s\t",samp,targs[m].mip,targs[m].miptype,targs[m].crispr);
			gzprintf(scounts,"%s\t%s\t%s\t%s\t%ld\n",current->contig,current->maploc,current->seq,avgqual(current->qual,current->seqcount,current->seq,finalqual),current->tagcount);
			current=current->next;
		}
	}
	return;
}

char*avgqual(double*curqual,long count,char*curseq,char*newqual)
{
	long q;
	for(q=0;q<strlen(curseq);q++)
	{
		newqual[q]=(char)(curqual[q]/count);
	}
	newqual[q]='\0';
	return newqual;
}

void freeseqs(struct miptarg*targs,long numtargs)
{
	long m;
	struct mipseq*temp;
	for(m=0;m<numtargs;m++)
	{
		while(targs[m].seqs!=NULL)
		{
			freetags(targs[m].seqs->tags);
			temp=targs[m].seqs->next;
			free(targs[m].seqs);
			targs[m].seqs=temp;
		}
	}
	return;
}

void freetags(struct moltag*taglist)
{
	struct moltag*temp;
	while(taglist!=NULL)
	{
		temp=taglist->next;
		free(taglist);
		taglist=temp;
	}
	return;
}

