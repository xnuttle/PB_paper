//Xander Nuttle
//get_guidecounts.c
//Call: ./get_guidecounts gzipped_finalseqs_file miptargets_file
//
//Generates a ".guidecounts" file containing counts of MIP capture events for each guide constuct for a single sample,
//taking the gzipped finalseqs file for that sample and a miptargets file as inputs. The miptargets file should only
//include MIPs targeting integrated guide construct sequences. The output ".guidecounts" file enables assessment of
//which guide constructs integrated into the genome of the sample analyzed.

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<zlib.h>
#define NLEN 200 //size of character vectors for storing names, etc.
#define LLEN 1500 //maximum length of single line of text in input finalseqs file

//set up structure to store MIP target information
struct miptarg
{
	char name[NLEN+1];
	char contig[NLEN+1];
};

long count_targs(FILE*mtargs);
void get_targ_info(FILE*mtargs,struct miptarg*targs,long*numconsts,long*nummips);
FILE*init_output(FILE*out,char*basename,long nummips);
void init_guides(struct miptarg*targs,long numconsts,long nummips,long numtargs,char gyds[numconsts][NLEN+1],long cownts[numconsts][nummips]);
long findguide(char*mipname,long numconsts,char gyds[numconsts][NLEN+1]);
void print_data(FILE*out,long numconsts,long nummips,char gyds[numconsts][NLEN+1],long cownts[numconsts][nummips]);

int main(int argc,char*argv[])
{
	//get sample name
	char sample[NLEN+1];
	strncpy(sample,*(argv+1),NLEN-13);
	sample[strchr(sample,'.')-sample]='\0';

	//determine the number of MIP targets and allocate memory to store MIP information
	FILE*miptargs=fopen(*(argv+2),"r");
	long ntargs=count_targs(miptargs);
	struct miptarg*targets;
	targets=(struct miptarg*)malloc(ntargs*sizeof(struct miptarg));

	//read in information on MIP targets; also, determine the number of guide constructs and the maximum number of MIPs targeting any guide construct
	long nconsts=1,nmips=1;	
	get_targ_info(miptargs,targets,&nconsts,&nmips);

	//set up output file
	FILE*guidecounts=init_output(guidecounts,sample,nmips);

	//allocate memory to store guide construct information and initialize this info
	char guides[nconsts][NLEN+1];
	long counts[nconsts][nmips];
	init_guides(targets,nconsts,nmips,ntargs,guides,counts);

	//read in final called MIP sequences and process them one by one
	gzFile*finalseqs=gzopen(*(argv+1),"r");
	char line[LLEN+1],mip[NLEN+1];
	long g,m,tagcount;
	sscanf(line,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s"); //process header line
	while(gzgets(finalseqs,line,LLEN))
	{
		sscanf(line,"%*s %s %*s %*s %*s %*s %*s %*s %ld %*s",mip,&tagcount);
		g=findguide(mip,nconsts,guides);
		if(g>=0)
		{
			for(m=0;m<nmips;m++)
			{
				sprintf(line,"%s_MIP_%04ld",guides[g],m+1);
				if(strncmp(mip,line,NLEN)==0)
				{
					counts[g][m]+=tagcount;
					break;
				}
			}
		}
	}	

	//print counts for each guide construct
	print_data(guidecounts,nconsts,nmips,guides,counts);

	//clean up and exit
	free(targets);
	gzclose(finalseqs);
	fclose(miptargs);
	fclose(guidecounts);
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
		numtargs++;
	fsetpos(mtargs,&pos);
	return numtargs;
}

void get_targ_info(FILE*mtargs,struct miptarg*targs,long*numconsts,long*nummips)
{
	long t=0,cmips=0;
	while(fscanf(mtargs,"%s %*s %s %*s %*s %*s %*s %*s %*s %*s",targs[t].name,targs[t].contig)==2)
	{
		if((t>0)&&(strncmp(targs[t].contig,targs[t-1].contig,NLEN)!=0))
		{
			if(cmips>(*nummips))
				(*nummips)=cmips;
			(*numconsts)++;
			cmips=1;
		}
		else
			cmips++;
		t++;
	}
	return;
}

FILE*init_output(FILE*out,char*basename,long nummips)
{
	char outname[NLEN+1];
	long m;
	sprintf(outname,"%s%s",basename,".guidecounts");
	out=fopen(outname,"w");
	fprintf(out,"Guide");
	for(m=0;m<nummips;m++)
		fprintf(out,"\tMIP%ld_Count",m+1);
	fprintf(out,"\n");
	return out;
}

void init_guides(struct miptarg*targs,long numconsts,long nummips,long numtargs,char gyds[numconsts][NLEN+1],long cownts[numconsts][nummips])
{
	long g,m,t;
	for(g=0;g<numconsts;g++)
	{
		for(m=0;m<nummips;m++)
			cownts[g][m]=-1;
	}
	g=-1;
	for(t=0;t<numtargs;t++)
	{
		if(strstr(targs[t].name,"0001")!=NULL)
		{
			g++;
			strncpy(gyds[g],targs[t].contig,NLEN);
			m=0;
			cownts[g][m]=0;
		}
		else
		{
			m++;
			cownts[g][m]=0;
		}
	}
	return;
}

long findguide(char*mipname,long numconsts,char gyds[numconsts][NLEN+1])
{
	long guide;
	for(guide=0;guide<numconsts;guide++)
	{
		if(strstr(mipname,gyds[guide])!=NULL)
			return guide;
	}
	return -1;
}

void print_data(FILE*out,long numconsts,long nummips,char gyds[numconsts][NLEN+1],long cownts[numconsts][nummips])
{
	long guide,mip;
	for(guide=0;guide<numconsts;guide++)
	{
		fprintf(out,"%s",gyds[guide]);
		for(mip=0;mip<nummips;mip++)
			fprintf(out,"\t%ld",cownts[guide][mip]);
		fprintf(out,"\n");
	}
	return;
}

