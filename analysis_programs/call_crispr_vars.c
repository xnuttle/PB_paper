//Xander Nuttle
//call_crispr_vars.c
//Call: ./call_crispr_vars gzipped_finalseqs_file crispr_sites_file <prime_variants_file> <coordinate_conversion_file> <(long)flank_length>
//
//This program analyzes finalized MIP sequence data together with information regarding potential CRISPR edits to call whether or not the
//sample analyzed acquired each edit. The crispr sites file "crispr_sites.crispr" should match the format of the same file used as input to
//detail_mip_targets_v5.c and should include crispr sites for all contigs analyzed. Its use here is to provide a list of CRISPR guides and
//prime edits that will be used to generate a structure to store data on their status in this particular sample. For a sample to be deemed as
//having a particular prime editing variant, the exact prime editing variant must be present at the location specified in the prime variants
//file and at least a specified number of bases (flank_length) on each side of the prime variant must be identical between the reference sequence
//and the called MIP sequence. If analyzing prime editing variants, the 4th input argument should be a file specifying chromosomes, corresponding
//contigs, and corresponding coordinates, e.g., one line might be "chrCHD8 chr14	chrcoord" where chrcoord is the base1 chr14 coordinate
//corresponding to the first base of the chrCHD8 contig. All contigs should be sequence from the '+' strand of a segment of the hg38 reference genome.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include<zlib.h>
#define NLEN 200 //size of character vectors for storing names, etc.
#define SLEN 500 //maximum length of each sequence array (and corresponding quality array)
#define LLEN 1500 //maximum length of single line of text in input finalseqs file
#define NFLK 5 //number of bases to check flanking each prime editing variant 

//set up structure to store data for CRISPR edits
struct edit
{
	char name[NLEN+1];
	char status[NLEN+1];
	long nindels;
};

//set up structure to store data for prime editing variants
struct pevariant
{
	char chr[NLEN+1];
	long coord;
	char ref[SLEN+1];
	char alt[SLEN+1];
	char name[NLEN+1];
};

//set up structure to store data for coordinate conversion
struct coordtable
{
	char contig[NLEN+1];
	char chr[NLEN+1];
	long chrcoord;
};

//set up structure to store data for each input sequence
struct input
{
	char mip[NLEN+1];
	char miptype[NLEN+1];
	char crispr[NLEN+1];
	char contig[NLEN+1];
	long maploc;
	char seq[SLEN+1];
	char qual[SLEN+1];
	char tagcount[NLEN+1];
	char tagfreq[NLEN+1];
};

long count_edits(FILE*elist);
void get_edits(FILE*elist,struct edit*edata);
long count_vars(FILE*vars);
void get_vars(FILE*vars,struct pevariant*pvars);
long count_contigs(FILE*cconvert);
void get_coords(FILE*cconvert,struct coordtable*cdata);
int getinput(gzFile*fseqs,struct input*iseq);
void process(struct input*iseq,struct edit*edata,long numedits,struct pevariant*pvars,long numvars,long flank,struct coordtable*cdata,long numcontigs);
long findedit(char*ename,struct edit*editdata,long ncrispr);
long findvar(char*vname,struct pevariant*variants,long nvariants);
long findcontig(char*cname,struct coordtable*ctigs,long nctigs);
void clearseqs(char*rseq,char*aseq,char*rseg,char*aseg,char*rlflk,char*alflk,char*rrflk,char*arflk);
void clearcoords(long*coords);
void parse_cs(char*cs,char*refaln,char*altaln,long coord,long*coordaln);
void alnncpy(char*dest1,char*dest2,char*src1,char*src2,long num1,long num2);
void print_data(char*samp,struct edit*edata,long numedits);

int main(int argc,char*argv[])
{
	//get sample name
	char sample[NLEN+1];
	strncpy(sample,*(argv+1),NLEN);
	sample[strchr(sample,'.')-sample]='\0';

	//determine the number of genotyped CRISPR edits, allocate memory to store editing information, read in CRISPR edits from input file, and initialize editing information
	FILE*editlist=fopen(*(argv+2),"r");
	long nedits=count_edits(editlist);
	struct edit*edits;
	edits=(struct edit*)malloc(nedits*sizeof(struct edit));
	get_edits(editlist,edits);

	//determine whether input includes file with information about prime editing variants; if so, check for input on flank length and read it in if present
	int primevars=0;
	long nflk=NFLK;
	if(argc>3)
	{
		primevars++;
		if(argc>5)
			nflk=strtol(*(argv+5),NULL,10);
	}

	//determine the number of prime editing variants, allocate memory to store variant information, and read in variant information from input file
	//also read in data for conversion from contig coordinates to chromosomal coordinates
	FILE*variants,*ctable;
	long nvars=0,ncontigs=0;
	struct pevariant*pevars=NULL;
	struct coordtable*contigs=NULL;
	if(primevars)
	{
		variants=fopen(*(argv+3),"r");
		nvars=count_vars(variants);
		pevars=(struct pevariant*)malloc(nvars*sizeof(struct pevariant));
		get_vars(variants,pevars);
		ctable=fopen(*(argv+4),"r");
		ncontigs=count_contigs(ctable);
		contigs=(struct coordtable*)malloc(ncontigs*sizeof(struct coordtable));
		get_coords(ctable,contigs);
	}

	//read in final called MIP sequences and process them one by one
	gzFile*finalseqs=gzopen(*(argv+1),"r");
	struct input inseq;
	getinput(finalseqs,&inseq); //process header line
	while(getinput(finalseqs,&inseq))
		process(&inseq,edits,nedits,pevars,nvars,nflk,contigs,ncontigs);
	
	//print status for each genotype CRISPR edit
	print_data(sample,edits,nedits);
	
	//clean up and exit
	free(edits);
	gzclose(finalseqs);
	fclose(editlist);
	if(primevars)
	{
		free(pevars);
		fclose(variants);
	}
	return 0;
}

long count_edits(FILE*elist)
{
	long numedits=0;
	char crchr[SLEN+1];
	fpos_t pos;
	fgetpos(elist,&pos);
	while(fscanf(elist,"%s %*s %*s",crchr)==1)
		numedits++;
	fsetpos(elist,&pos);
	return numedits;
}

void get_edits(FILE*elist,struct edit*edata)
{
	long e=0;
	while(fscanf(elist,"%*s %*s %s",edata[e].name)==1)
	{
		strncpy(edata[e].status,"uncallable\0",NLEN);
		edata[e].nindels=0;
		e++;
	}
	return;
}

long count_vars(FILE*vars)
{
	long numvars=0;
	char varname[NLEN+1];
	fpos_t pos;
	fscanf(vars,"%*s %*s %*s %*s %*s");
	fgetpos(vars,&pos);
	while(fscanf(vars,"%*s %*s %*s %*s %s",varname)==1)
		numvars++;
	fsetpos(vars,&pos);
	return numvars;
}

void get_vars(FILE*vars,struct pevariant*pvars)
{
	long v=0;
	while(fscanf(vars,"%s %ld %s %s %s",pvars[v].chr,&(pvars[v].coord),pvars[v].ref,pvars[v].alt,pvars[v].name)==5)
		v++;
	return;
}

long count_contigs(FILE*cconvert)
{
	long numcontigs=0;
	char cname[NLEN+1];
	fpos_t pos;
	fgetpos(cconvert,&pos);
	while(fscanf(cconvert,"%s %*s %*s",cname)==1)
		numcontigs++;
	fsetpos(cconvert,&pos);
	return numcontigs;
}

void get_coords(FILE*cconvert,struct coordtable*cdata)
{
	long c=0;
	while(fscanf(cconvert,"%s %s %ld",cdata[c].contig,cdata[c].chr,&(cdata[c].chrcoord))==3)
		c++;
	return;
}

int getinput(gzFile*fseqs,struct input*iseq)
{
	char line[LLEN+1];
	int scanned=0;
	if(gzgets(fseqs,line,LLEN))
		scanned+=sscanf(line,"%*s %s %s %s %s %ld %s %s %s %s",iseq->mip,iseq->miptype,iseq->crispr,iseq->contig,&(iseq->maploc),iseq->seq,iseq->qual,iseq->tagcount,iseq->tagfreq);
	return (scanned==9);
}

void process(struct input*iseq,struct edit*edata,long numedits,struct pevariant*pvars,long numvars,long flank,struct coordtable*cdata,long numcontigs)
{
	char*cr;
	char edit[NLEN+1];
	long e,v,c,chrloc,j;
	char refseq[SLEN+1]="",altseq[SLEN+1]="",refseg[SLEN+1]="",altseg[SLEN+1]="",reflflk[SLEN+1]="",altlflk[SLEN+1]="",refrflk[SLEN+1]="",altrflk[SLEN+1]="";
	long refcoord[SLEN+1];
	cr=strtok(iseq->crispr,"/");
	do
	{
		if((strncmp(cr,"none",NLEN)==0)||(strncmp(cr,"PB-",3)==0)||(strncmp(cr,"plasmid",NLEN)==0))
		{
			cr=strtok(NULL,"/");
			continue;
		}
		e=findedit(cr,edata,numedits);
		if(e==-1)
		{
			cr=strtok(NULL,"/");
			continue;
		}
		v=findvar(cr,pvars,numvars);
		if(v==-1)
		{
			if(strpbrk(iseq->seq,"+-")!=NULL)
			{
				strncpy(edata[e].status,"has_indel\0",NLEN);
				edata[e].nindels++;
			}
			else if(strncmp(edata[e].status,"uncallable",10)==0)
				strncpy(edata[e].status,"no_indel\0",NLEN);
			cr=strtok(NULL,"/");
			continue;
		}
		c=findcontig(iseq->contig,cdata,numcontigs);
		chrloc=(iseq->maploc)+cdata[c].chrcoord-1;
		clearseqs(refseq,altseq,refseg,altseg,reflflk,altlflk,refrflk,altrflk);
		clearcoords(refcoord);
		parse_cs(iseq->seq,refseq,altseq,chrloc,refcoord);
		for(j=0;j<SLEN;j++)
		{
			if(refcoord[j]==pvars[v].coord)
			{
				alnncpy(refseg,altseg,refseq+j,altseq+j,strlen(pvars[v].ref),strlen(pvars[v].alt));
				while((j-flank)<0)
					flank--;
				alnncpy(reflflk,altlflk,refseq+j-flank,altseq+j-flank,flank,flank);				
			}
			if(refcoord[j]==(pvars[v].coord+strlen(pvars[v].ref)))
			{
				alnncpy(refrflk,altrflk,refseq+j,altseq+j,flank,flank);				
				break;
			}
		}	
		if((strncmp(refseg,pvars[v].ref,SLEN)==0)&&(strncmp(altseg,pvars[v].alt,SLEN)==0)&&(strncmp(reflflk,altlflk,flank)==0)&&(strncmp(refrflk,altrflk,flank)==0))
			strncpy(edata[e].status,"has_prime_edit\0",NLEN);
		else if(strncmp(edata[e].status,"uncallable",10)==0)
			strncpy(edata[e].status,"no_prime_edit\0",NLEN);
		cr=strtok(NULL,"/");
	}
	while(cr!=NULL);
	return;
}

long findedit(char*ename,struct edit*editdata,long ncrispr)
{
	long ed;
	for(ed=0;ed<ncrispr;ed++)
	{
		if(strncmp(editdata[ed].name,ename,NLEN)==0)
			return ed;
	}
	return -1;
}

long findvar(char*vname,struct pevariant*variants,long nvariants)
{
	long vnum;
	for(vnum=0;vnum<nvariants;vnum++)
	{
		if(strncmp(variants[vnum].name,vname,NLEN)==0)
			return vnum;
	}
	return -1;
}

long findcontig(char*cname,struct coordtable*ctigs,long nctigs)
{
	long cnum;
	for(cnum=0;cnum<nctigs;cnum++)
	{
		if(strncmp(ctigs[cnum].contig,cname,NLEN)==0)
			return cnum;
	}
	return -1;
}

void clearseqs(char*rseq,char*aseq,char*rseg,char*aseg,char*rlflk,char*alflk,char*rrflk,char*arflk)
{
	long i;
	for(i=0;i<SLEN;i++)
	{
		rseq[i]='\0';
		aseq[i]='\0';
		rseg[i]='\0';
		aseg[i]='\0';
		rlflk[i]='\0';
		alflk[i]='\0';
		rrflk[i]='\0';
		arflk[i]='\0';
	}
	return;
}

void clearcoords(long*coords)
{
	long i;
	for(i=0;i<SLEN;i++)
		coords[i]=0;
	return;
}

void parse_cs(char*cs,char*refaln,char*altaln,long coord,long*coordaln)
{
	long i=0,a=0;
	char alntype;
	while(cs[i]!='\0')
	{
		if(!(isalpha(cs[i])))
		{
			alntype=cs[i];
			i++;
		}
		if(alntype=='=')
		{
			refaln[a]=cs[i];
			altaln[a]=cs[i];
			coordaln[a]=coord;
			coord++;
		}
		else if(alntype=='+')
		{
			refaln[a]='-';
			altaln[a]=toupper(cs[i]);
			coordaln[a]=-1;
		}
		else if(alntype=='-')
		{
			refaln[a]=toupper(cs[i]);
			altaln[a]='-';
			coordaln[a]=coord;
			coord++;
		}
		else if(alntype=='*')
		{
			refaln[a]=toupper(cs[i]);
			i++;
			altaln[a]=toupper(cs[i]);
			coordaln[a]=coord;
			coord++;
		}
		i++;
		a++;
	}
	refaln[a]='\0';
	altaln[a]='\0';
	return;
}

void alnncpy(char*dest1,char*dest2,char*src1,char*src2,long num1,long num2)
{
	long c1=0,c2=0,d1=0,d2=0,s1=0,s2=0;
	while((c1<num1)||(c2<num2))
	{
		if((src1[s1]=='\0')||(src2[s2]=='\0'))
			break;
		if(src1[s1]!='-')
		{
			dest1[d1]=src1[s1];
			d1++;
			c1++;
		}
		if(src2[s2]!='-')
		{
			dest2[d2]=src2[s2];
			d2++;
			c2++;
		}
		s1++;
		s2++;
	}
	dest1[d1]='\0';
	dest2[d2]='\0';
	return;
}

void print_data(char*samp,struct edit*edata,long numedits)
{
	long e;
	for(e=0;e<numedits;e++)
		printf("%s\t%s\t%s\t%ld\n",samp,edata[e].name,edata[e].status,edata[e].nindels);
	return;
}

