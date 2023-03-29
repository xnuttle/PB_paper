//Xander Nuttle
//detail_mip_targets_v5.c
//Call: ./detail_mip_targets_v5 mip_sequences_file contig_fasta_file <snvs_bedgraph_file.snvs> <crispr_sites.crispr> <(long)coord>
//
//Generates a ".miptargets" file containing information regarding a set of input MIPs and their associated targets within a contig sequence.
//
//This version of the detail_mip_targets program was created because the output from previous versions, while useful for paralog-specific
//copy number analysis, did not contain important information for the current primary application of MIPs: genotyping copy number and sequence
//variation within regions of unique sequence after attempted CRISPR editing in iPSCs corresponding to MGH2069 or GM08330. The new output
//is tailored to this purpose and is used by other programs in my automated MIP analysis pipeline.
//
//To annotate MIP targets with respect to SNVs in MGH2069 and GM08330 cell lines and/or with respect to guide RNA cut sites and/or
//prime editing sites, files ending in ".snvs" and ".crispr" can be used as optional inputs. These inputs must be sorted by chromosomal
//coordinate (column 2) or else MIP targets may not get annotated properly. If any of these optional input files are used, the program
//requires a final command line argument specifying the chromosomal coordinate corresponding to the first base of the contig sequence.

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#define SLEN 200 //size of character vectors for storing names, etc.
#define CLEN 20000000 //size of character vector for storing contig sequences (max contig size = 20 Mb)
#define MSIZE 200 //maximum distance between MIP arms in contig sequence for intervening sequence to be considered a valid MIP target

//set up structure to store MIP information
struct mip
{
	char name[SLEN+1];
	char seq[SLEN+1];
	char chr[SLEN+1];
	long start;
	long end;
	char type;
	char crispr[SLEN+1];
	char strand;
	long armlen;
	long targlen;
};

//set up structure to store contig information
struct contig
{
	char name[SLEN+1];
	char seq[CLEN+1];
	char seqrc[CLEN+1];
	long length;
};

//set up structure to store SNV information
struct snvtable
{
	long coord;
	char type;
};

//set up structure to store CRISPR target information
struct crtarg
{
	double coord;
	char name[SLEN+1];	
};

long count_mips(FILE*mseqs);
void get_mipseqs(FILE*mseqs,struct mip*mipinfo);
void get_seq(FILE*fa,struct contig*c);
void revcomp(char string[]);
long count_snvs(FILE*bgraph);
void get_snvs(FILE*bgraph,struct snvtable*snps,long cstart);
long count_crtargs(FILE*crinfo);
void get_crtargs(FILE*crinfo,struct crtarg*targs,long cstart);
void detail_mip(struct mip*mipinfo,struct contig*c,struct snvtable*snps,struct crtarg*targs,long nummips,long numsnvs,long numtargs);
void findtarg(struct mip*mipdata,long mnum,struct contig*sequence);
void classify(struct mip*mipdata,long mnum,struct snvtable*snvdata,long nsnps);
void annotate(struct mip*mipdata,long mnum,struct crtarg*targdata,long ntargs);
FILE*init_output(FILE*out,struct contig*c);
void print_data(FILE*out,struct mip*mipdata,long nummips);

int main(int argc,char*argv[])
{
	//determine the number of MIPs and allocate memory to store MIP information
	FILE*mipseqs=fopen(*(argv+1),"r");
	long nmips=count_mips(mipseqs);
	struct mip*mips;
	mips=(struct mip*)malloc(nmips*sizeof(struct mip));

	//read in MIP sequences
	get_mipseqs(mipseqs,mips);

	//allocate memory to store contig information
	struct contig*chr;
	chr=(struct contig*)malloc(sizeof(struct contig));	

	//read in contig sequence from fasta file, generate its reverse complements, and calculate its length
	FILE*fasta=fopen(*(argv+2),"r");
	get_seq(fasta,chr);

	//determine whether input includes file with locations of SNVs
	int snvinput=0;
	if((argc>3)&&(strstr(*(argv+3),"snvs")!=NULL))
		snvinput++;	

	//determine whether input includes file with locations of CRISPR targets
	int crinput=0;
	if(strstr(*(argv+argc-2),"crispr")!=NULL)
		crinput++;	

	//if input includes file(s) with SNV and/or CRISPR info, get chromosomal coordinate corresponding to first base in contig sequence from command line
	long coord;
	if(snvinput||crinput)
		coord=strtol(*(argv+argc-1),NULL,10);	

	//determine the number of SNVs, allocate memory to store SNV information, and read in SNV information from input file
	FILE*bedgraph;
	struct snvtable*snvs=NULL;
	long nsnvs=0;
	if(snvinput)
	{
		bedgraph=fopen(*(argv+3),"r");
		nsnvs=count_snvs(bedgraph);
		snvs=(struct snvtable*)malloc(nsnvs*sizeof(struct snvtable));
		get_snvs(bedgraph,snvs,coord);
	}

	//determine the number of CRISPR targets, allocate memory to store CRISPR target information, and read in CRISPR target information from input file
	FILE*crfile;
	struct crtarg*crtargs=NULL;
	long ncrispr=0;
	if(crinput)
	{
		crfile=fopen(*(argv+argc-2),"r");
		ncrispr=count_crtargs(crfile);
		crtargs=(struct crtarg*)malloc(ncrispr*sizeof(struct crtarg));
		get_crtargs(crfile,crtargs,coord);
	}

	//get data for each MIP
	detail_mip(mips,chr,snvs,crtargs,nmips,nsnvs,ncrispr);

	//set up output file
	FILE*miptargets=init_output(miptargets,chr);	

	//print MIP data to output	
	print_data(miptargets,mips,nmips);

	//clean up and exit
	free(mips);
	free(chr);
	fclose(mipseqs);
	fclose(fasta);
	if(snvinput)
	{
		fclose(bedgraph);
		free(snvs);
	}
	if(crinput)
	{
		fclose(crfile);
		free(crtargs);
	}
	return 0;
}

long count_mips(FILE*mseqs)
{
	long nummips=0;
	char mipseq[SLEN];
	fpos_t pos;
	fgetpos(mseqs,&pos);
	while(fscanf(mseqs,"%s",mipseq)==1)
		nummips++;
	fsetpos(mseqs,&pos);
	return nummips;
}

void get_mipseqs(FILE*mseqs,struct mip*mipinfo)
{
	long m=0;
	while(fscanf(mseqs,"%s",mipinfo[m].seq)==1)
		m++;
	return;
}

void get_seq(FILE*fa,struct contig*c)
{
	fscanf(fa,"%*c %s",c->name);	
	long i=0;
	char ch;
	while((ch=getc(fa))!=EOF)
	{
		if(isalpha(ch))
		{
			c->seq[i]=tolower(ch);
			i++;
		}
	}
	c->seq[i]='\0';
	strncpy(c->seqrc,c->seq,CLEN);
	revcomp(c->seqrc);
	c->length=strlen(c->seq);
	return;
}

void revcomp(char string[])
{
	long a=strlen(string);
	long b,c;
	char base;
	char newstring[a+1];
	for(b=0;b<a;b++)
	{
		newstring[b]=string[a-1-b];
	}
	for(c=0;c<a;c++)
	{
		base=newstring[c];
		switch(base)
		{
			case 'a': base='t'; break;
			case 'c': base='g'; break;
			case 'g': base='c'; break;
			case 't': base='a'; break;
			default: base='n'; break;
		}
		newstring[c]=base;
	}
	newstring[a]='\0';
	strncpy(string,newstring,CLEN);
}

long count_snvs(FILE*bgraph)
{
	long numsnvs=0;
	char snvchr[SLEN+1];
	fpos_t pos;
	fgetpos(bgraph,&pos);
	while(fscanf(bgraph,"%s %*s %*s %*s %*s",snvchr)==1)
		numsnvs++;
	fsetpos(bgraph,&pos);
	return numsnvs;
}

void get_snvs(FILE*bgraph,struct snvtable*snps,long cstart)
{
	long chrcoord,s=0;
	int mgh,gm;
	while(fscanf(bgraph,"%*s %*s %ld %d %d",&chrcoord,&mgh,&gm)==3)
	{
		snps[s].coord=chrcoord-cstart+1;
		if((mgh)&&(gm))
			snps[s].type='B';
		else if(mgh)
			snps[s].type='M';
		else
			snps[s].type='E';
		s++;
	}
	return;
}

long count_crtargs(FILE*crinfo)
{
	long numtargs=0;
	char crchr[SLEN+1];
	fpos_t pos;
	fgetpos(crinfo,&pos);
	while(fscanf(crinfo,"%s %*s %*s",crchr)==1)
		numtargs++;
	fsetpos(crinfo,&pos);
	return numtargs;
}

void get_crtargs(FILE*crinfo,struct crtarg*targs,long cstart)
{
	long t=0;
	double chrcoord;
	char targid[SLEN+1];
	while(fscanf(crinfo,"%*s %lf %s",&chrcoord,targid)==2)
	{
		targs[t].coord=chrcoord-cstart+1;
		strncpy(targs[t].name,targid,SLEN);
		t++;
	}
	return;
}

void detail_mip(struct mip*mipinfo,struct contig*c,struct snvtable*snps,struct crtarg*targs,long nummips,long numsnvs,long numtargs)
{
	long m;
	for(m=0;m<nummips;m++)
	{
		//detemine MIP name and contig targeted
		sprintf(mipinfo[m].name,"%s_MIP_%04ld",c->name,m+1);
		strncpy(mipinfo[m].chr,c->name,SLEN);
		
		//determine MIP start and end coordinates (base1), strand, length of targeting arm encountered first in sequence, and length of target sequence
		findtarg(mipinfo,m,c);

		//determine MIP type (designates whether MIP targets one or more SNVs present in MGH2069 and/or GM08330 genomes
		classify(mipinfo,m,snps,numsnvs);
		
		//determine whether MIP target sequence includes one or more CRISPR sites
		annotate(mipinfo,m,targs,numtargs);
	}
	return;
}

void findtarg(struct mip*mipdata,long mnum,struct contig*sequence)
{
	char ext[SLEN+1],lig[SLEN+1];
	long b;
	for(b=0;b<SLEN;b++)
	{
		ext[b]='\0';
		lig[b]='\0';
	}
	int found=0,rc=0;
	char*ligloc,*extloc,*searchloc;
	strncpy(ext,strrchr(mipdata[mnum].seq,'N')+1,SLEN);
	strncpy(lig,mipdata[mnum].seq,strchr(mipdata[mnum].seq,'C')-mipdata[mnum].seq);
	searchloc=sequence->seq;	
	while(!(found))
	{
		extloc=strstr(searchloc,ext);
		if(extloc==NULL)
		{
			searchloc=sequence->seqrc;
			rc++;
			continue;
		}
		searchloc=extloc+1;
		while((ligloc=strstr(searchloc,lig))!=NULL)
		{
			if((ligloc-extloc)<=MSIZE)
			{
				found++;
				break;
			}
			searchloc=ligloc+1;
		}
		searchloc=extloc+1;
	}
	if(!(rc))
	{
		mipdata[mnum].start=extloc-sequence->seq+1;
		mipdata[mnum].end=ligloc-sequence->seq+strlen(lig);
		mipdata[mnum].armlen=strlen(ext);
		mipdata[mnum].strand='+';
		mipdata[mnum].targlen=mipdata[mnum].end-mipdata[mnum].start+1-(mipdata[mnum].armlen+strlen(lig));
	}
	else
	{
		mipdata[mnum].start=sequence->length-(ligloc-sequence->seqrc+strlen(lig))+1;
		mipdata[mnum].end=sequence->length-(extloc-sequence->seqrc);
		mipdata[mnum].armlen=strlen(lig);
		mipdata[mnum].strand='-';
		mipdata[mnum].targlen=mipdata[mnum].end-mipdata[mnum].start+1-(mipdata[mnum].armlen+strlen(ext));
	}
	return;
}

void classify(struct mip*mipdata,long mnum,struct snvtable*snvdata,long nsnps)
{
	long s;
	int mgh=0,gm=0;
	mipdata[mnum].type='N';
	for(s=0;s<nsnps;s++)
	{
		if((snvdata[s].coord>=mipdata[mnum].start)&&(snvdata[s].coord<=mipdata[mnum].end))
		{
			if(snvdata[s].type=='M')
				mgh=1;
			else if(snvdata[s].type=='E')
				gm=1;
			else
			{
				mipdata[mnum].type='B';
				return;
			}
		}
		if(snvdata[s].coord>mipdata[mnum].end)
			break;
	}
	if((mgh==1)&&(gm==0))
		mipdata[mnum].type='M';
	else if((mgh==0)&&(gm==1))
		mipdata[mnum].type='E';
	else if((mgh==1)&&(gm==1))
		mipdata[mnum].type='B';
	return;
}

void annotate(struct mip*mipdata,long mnum,struct crtarg*targdata,long ntargs)
{
	long t;
	strncpy(mipdata[mnum].crispr,"none",SLEN);
	for(t=0;t<ntargs;t++)
	{
		if((targdata[t].coord>=mipdata[mnum].start)&&(targdata[t].coord<=mipdata[mnum].end))
		{
			if(strncmp(mipdata[mnum].crispr,"none",SLEN)==0)
				strncpy(mipdata[mnum].crispr,targdata[t].name,SLEN);
			else
			{
				strncat(mipdata[mnum].crispr,"/",SLEN);
				strncat(mipdata[mnum].crispr,targdata[t].name,SLEN);
			}
		}
		if(targdata[t].coord>mipdata[mnum].end)
			return;
	}
	return;
}

FILE*init_output(FILE*out,struct contig*c)
{
	char outname[SLEN+1];
	sprintf(outname,"%s%s",c->name,".miptargets");
	out=fopen(outname,"w");
	fprintf(out,"Name\tSequence\tContig\tStart\tEnd\tType\tCRISPR\tStrand\tArm1Length\tTargetLength\n");
	return out;
}

void print_data(FILE*out,struct mip*mipdata,long nummips)
{
	long m;
	for(m=0;m<nummips;m++)
	{
		fprintf(out,"%s\t%s\t%s\t%ld\t%ld\t%c\t%s\t%c\t%ld\t%ld\n",mipdata[m].name,mipdata[m].seq,mipdata[m].chr,mipdata[m].start,mipdata[m].end,mipdata[m].type,mipdata[m].crispr,mipdata[m].strand,mipdata[m].armlen,mipdata[m].targlen);
	}
	return;
}

