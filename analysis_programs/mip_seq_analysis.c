//Xander Nuttle
//mip_seq_analysis.c
//Call: ./mip_seq_analysis gzipped_sam_file miptargets_file <(double)mapping_location_wiggle_room>
//
//This program analyzes reads in a gzipped sam mapping output file generated in a MIP experiment.
//It assigns each read to a MIP target of interest, annotates sequence variation in an easily parsed
//format, and outputs the annotated sequence and quality string alongside the sample name, information
//corresponding to the assigned MIP target, and the molecular tag. The gzipped output file should then
//be analyzed together with other output files from the same sample to count the number of capture events
//corresponding to each observed sequence at each MIP target of interest.  
//
//The default value for the mapping location wiggle room (4.5) generally works well, though in some cases
//it may be appropriate to adjust this value via a third optional command line argument. For example, if
//you have multiple MIPs with targets shifted by a few bases or less, setting this value to zero would
//allow you to keep sequences corresponding to these nearby MIP targets separate for further analysis.

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<zlib.h>
#define NLEN 200 //size of character vectors for storing names, etc.
#define SLEN 500 //size of character vectors for storing sequence and quality strings 
#define TLEN 8 //length of molecular tag sequences
#define LLEN 1501 //maximum length of single line of text in mapping output (gzipped sam) file + 1
#define MWIG 4.5 //default mapping location wiggle room

//set up structure to store MIP target information
struct miptarg
{
	char name[NLEN+1];
	char contig[NLEN+1];
	long start;
	char type;
	char crispr[NLEN+1];
	long tstart;
	long armlen;
	long tlength;
};

//set up structure to store read information
struct readdata
{
	char contig[NLEN+1];
	long maploc;
	char cigar[SLEN+1];
	char seq[SLEN+1];
	char qual[SLEN+1];
	char md[SLEN+1];
	char tag[TLEN+1];
};

gzFile* init_output(gzFile*mseqs,char*basename);
long count_targs(FILE*mtargs);
void get_targ_info(FILE*mtargs,struct miptarg*targs);
int getread(gzFile*samgz,struct readdata*reed);
void parseread(struct readdata*reed,struct miptarg*targs,long numtargs,gzFile*mseqs,char*samp,double wigg);
long findtarg(char*chr,long coord,struct miptarg*miptargets,long numtargets,double wig);
void parse_aln(char*cigar,char*md,char*readseq,char*readqual,long maploc,long targloc,long targlength,char*csseq,char*csqual);
void makecs(char aln,long nbases,long*t_index,long*r_index,char*mdtag,char*rseq,char*rqual,char*seqcs,char*qualcs,long targlen);
void parse_mapped(long num_bases,long*tindex,long*rindex,char*mdstring,char*read_seq,char*read_qual,char*cs_seq,char*cs_qual,long findex);
void parse_ins(long num_bases,long*tindex,long*rindex,char*read_seq,char*read_qual,char*cs_seq,char*cs_qual,long findex);
void parse_del(long num_bases,long*tindex,long*rindex,char*mdstring,char*read_qual,char*cs_seq,char*cs_qual,long findex);
void parse_clipped(long num_bases,long*rindex,char*read_seq,char*read_qual);
char*lowercase(char*string);

int main(int argc,char*argv[])
{
	//get sample name
	char sample[NLEN+1];
	strncpy(sample,*(argv+1),NLEN-11);
	sample[strchr(sample,'.')-sample]='\0';

	//set up output file
	gzFile*mipseqs=init_output(mipseqs,sample);

	//determine the number of MIP targets and allocate memory to store MIP information
	FILE*miptargs=fopen(*(argv+2),"r");
	long ntargs=count_targs(miptargs);
	struct miptarg*targets;
	targets=(struct miptarg*)malloc(ntargs*sizeof(struct miptarg));
	
	//read in information on MIP targets and link MIP targets to guide RNAs
	get_targ_info(miptargs,targets);

	//get value of mapping location wiggle room from command line
	double wiggle=MWIG;
	if(argc==4)
		wiggle=strtod(*(argv+3),NULL);

	//open gzipped sam file and process reads one by one
	gzFile*sam=gzopen(*(argv+1),"r");
	struct readdata read;
	while(getread(sam,&read))
		parseread(&read,targets,ntargs,mipseqs,sample,wiggle);

	//clean up and exit
	free(targets);
	gzclose(mipseqs);
	gzclose(sam);
	fclose(miptargs);
	return 0;
}

gzFile* init_output(gzFile*mseqs,char*basename)
{
	char outname[NLEN+1];
	sprintf(outname,"%s%s",basename,".mipseqs.gz\0");
	mseqs=gzopen(outname,"w");
	gzprintf(mseqs,"Sample\tMIP\tType\tCRISPR\tContig\tCoordinate\tSequence\tQuality\tTag\n");
	return mseqs;
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

void get_targ_info(FILE*mtargs,struct miptarg*targs)
{
	long m=0;
	while(fscanf(mtargs,"%s %*s %s %ld %*s %c %s %*s %ld %ld",targs[m].name,targs[m].contig,&(targs[m].start),&(targs[m].type),targs[m].crispr,&(targs[m].armlen),&(targs[m].tlength))==7)
	{
		targs[m].tstart=targs[m].start+targs[m].armlen;
		m++;
	}	
	return;
}

int getread(gzFile*samgz,struct readdata*reed)
{
	char line[LLEN],finalmd[SLEN+1];;
	int parsed=0;
	if(gzgets(samgz,line,LLEN-1))
	{
		parsed+=sscanf(line,"%*s %*s %s %ld %*s %s %*s %*s %*s %s %s %*s %s",reed->contig,&(reed->maploc),reed->cigar,reed->seq,reed->qual,reed->md);
		sprintf(finalmd,"%s",reed->md+5);
		strcpy(reed->md,finalmd);		
		parsed+=(strncpy(reed->tag,strstr(line,"MI:Z:$")+6,TLEN)!=NULL);
		reed->tag[TLEN]='\0';
	}
	return (parsed==7);
}

void parseread(struct readdata*reed,struct miptarg*targs,long numtargs,gzFile*mseqs,char*samp,double wigg)
{
	long m=findtarg(reed->contig,reed->maploc,targs,numtargs,wigg);
	char finalseq[SLEN+1],finalqual[SLEN+1];
	long i;
	for(i=0;i<SLEN;i++)
	{
		finalseq[i]='\0';
		finalqual[i]='\0';
	}
	finalseq[i]='\0';
	finalqual[i]='\0';
	if(m>=0)
	{
		parse_aln(reed->cigar,reed->md,reed->seq,reed->qual,reed->maploc,targs[m].tstart,targs[m].tlength,finalseq,finalqual);
		gzprintf(mseqs,"%s\t%s\t%c\t%s\t%s\t%ld\t%s\t%s\t%s\n",samp,targs[m].name,targs[m].type,targs[m].crispr,reed->contig,targs[m].tstart,finalseq,finalqual,reed->tag);
	}
	return;
}

long findtarg(char*chr,long coord,struct miptarg*miptargets,long numtargets,double wig)
{
	long mnum;
	for(mnum=0;mnum<numtargets;mnum++)
	{
		if((strncmp(miptargets[mnum].contig,chr,NLEN)==0)&&((coord-wig)<=miptargets[mnum].start)&&((coord+wig)>=miptargets[mnum].start))
			break;
	}
	if(mnum==numtargets)
		mnum=-1;
	return mnum;
}

void parse_aln(char*cigar,char*md,char*readseq,char*readqual,long maploc,long targloc,long targlength,char*csseq,char*csqual)
{
	long targ_index=maploc-targloc; //target index gives relation of read base to targeted contig bases
	long read_index=0; //read index gives position in read sequence
	long cs_index=0; //cs index gives position in output cs-formatted sequence and quality strings
	long numbases; //number of bases mapped, soft clipped, inserted, or deleted
	char alignment; //alignment type 'M', 'I', 'D', or 'S'
	char newcigar[SLEN+1];
	while(sscanf(cigar,"%ld %c",&numbases,&alignment)==2)
	{
		makecs(alignment,numbases,&targ_index,&read_index,md,readseq,readqual,csseq,csqual,targlength);
		sprintf(newcigar,"%s",strchr(cigar,alignment)+1);
		strcpy(cigar,newcigar);
	}
	return;
}

void makecs(char aln,long nbases,long*t_index,long*r_index,char*mdtag,char*rseq,char*rqual,char*seqcs,char*qualcs,long targlen)
{
	switch(aln)
	{
		case 'M': parse_mapped(nbases,t_index,r_index,mdtag,rseq,rqual,seqcs,qualcs,targlen-1); break;
		case 'I':	parse_ins(nbases,t_index,r_index,rseq,rqual,seqcs,qualcs,targlen-1); break;
		case 'D':	parse_del(nbases,t_index,r_index,mdtag,rqual,seqcs,qualcs,targlen-1); break;
		case 'S':	parse_clipped(nbases,r_index,rseq,rqual); break;
	}
	return;
}

void parse_mapped(long num_bases,long*tindex,long*rindex,char*mdstring,char*read_seq,char*read_qual,char*cs_seq,char*cs_qual,long findex)
{
	long b,nmatch,base_read,newtract=1;
	char rbase;
	char mapseq[4],newmd[SLEN+1];
	for(b=0;b<num_bases;b++)
	{
		base_read=0;
		if(isdigit(mdstring[0]))
		{
			if(sscanf(mdstring,"%ld %s",&nmatch,newmd)==1)
				strncpy(newmd,"\0",1);
			if(nmatch>0)
			{
				if((*tindex>=0)&&(*tindex<=findex))
				{
					if(newtract)
					{
						strncat(cs_seq,"=",1);
						strncat(cs_qual,"\"",1);
						newtract=0;
					}
					strncat(cs_seq,read_seq+(*rindex),1);
					strncat(cs_qual,read_qual+(*rindex),1);
				}
				base_read=1;
				nmatch--;
				sprintf(mdstring,"%ld%s",nmatch,newmd);
			}
			if(nmatch==0)
				strcpy(mdstring,newmd);
		}		
		if((isalpha(mdstring[0]))&&(!(base_read)))
		{
			sscanf(mdstring,"%c %s",&rbase,newmd);
			if((*tindex>=0)&&(*tindex<=findex))
			{
				sprintf(mapseq,"*%c%c",tolower(rbase),tolower(read_seq[*rindex]));
				strncat(cs_seq,mapseq,3);
				sprintf(mapseq,"\"%c%c",read_qual[*rindex],read_qual[*rindex]);
				strncat(cs_qual,mapseq,3);
			}
			newtract=1;
			strcpy(mdstring,newmd);
		}
		(*tindex)++; //increment index corresponding to target bases
		(*rindex)++; //increment index corresponding to read bases
	}
	return;
}

void parse_ins(long num_bases,long*tindex,long*rindex,char*read_seq,char*read_qual,char*cs_seq,char*cs_qual,long findex)
{
	char insseq[num_bases+1],insqual[num_bases+1];
	strncpy(insseq,read_seq+(*rindex),num_bases);
	strncpy(insqual,read_qual+(*rindex),num_bases);
	insseq[num_bases]='\0';
	insqual[num_bases]='\0';
	if((*tindex>=0)&&(*tindex<=findex))
	{
		strncat(cs_seq,"+",1);
		strncat(cs_seq,lowercase(insseq),strlen(insseq));
		strncat(cs_qual,"\"",1);
		strncat(cs_qual,insqual,strlen(insqual));		
	}
	(*rindex)+=num_bases; //increment index corresponding to read bases only
	return;
}

void parse_del(long num_bases,long*tindex,long*rindex,char*mdstring,char*read_qual,char*cs_seq,char*cs_qual,long findex)
{
	char delseq[num_bases+1],delqual[num_bases+1],newmd[SLEN+1];
	strncpy(delseq,mdstring+1,num_bases);
	delseq[num_bases]='\0';
	delqual[num_bases]='\0';
	long b;
	for(b=0;b<num_bases;b++)
	{
		if((*tindex<0)||(*tindex>findex))
		{
			delseq[b]='\t';
			delqual[b]='\t';
		}
		else
			delqual[b]=(read_qual[*rindex-1]+read_qual[*rindex])/2;
		(*tindex)++; //increment index corresponding to target bases only
	}	
	if((sscanf(delseq,"%s",delseq)==1)&&(sscanf(delqual,"%s",delqual)==1))
	{
		strncat(cs_seq,"-",1);
		strncat(cs_seq,lowercase(delseq),strlen(delseq));
		strncat(cs_qual,"\"",1);
		strncat(cs_qual,delqual,strlen(delqual));
	}
	sprintf(newmd,"%s",mdstring+strlen(delseq)+1);
	strcpy(mdstring,newmd);
	return;	
}

void parse_clipped(long num_bases,long*rindex,char*read_seq,char*read_qual)
{
	long b;
	for(b=0;b<num_bases;b++)
		(*rindex)++;
	return;
}

char*lowercase(char*string)
{
	long i;
	for(i=0;i<strlen(string);i++)
		string[i]=tolower(string[i]);
	return string;
}

//parse_aln, makecs, parse_mapped, parse_ins, parse_del, parse_clipped, and lowercase are functions to parse through CIGAR string and MD tag
//corresponding to a single read, generating a new sequence string in cs format and a corresponding new quality string; these strings can be
//easily parsed to quickly analyze the alignment of the read to the reference sequence it mapped to for SNV mutations and indels
//
//see https://github.com/lh3/minimap2#cs for details of the cs format; this program annotates the alignment using the cs long format
//
//the CIGAR string and MD tag must be parsed simultaneously with the sequence and quality strings because they inform each other's interpretation, see http://davetang.org/muse/2011/01/28/perl-and-sam/
//
//this code has been tested on several gzipped sam inputs and reproducibly outputs annotated sequences equivalent to those produced using the
//parse_cigar_and_md function from my old MIP sequence annotation program used to analyze SRGAP2 MIP sequence data

