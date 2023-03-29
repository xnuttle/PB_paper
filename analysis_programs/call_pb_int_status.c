//Xander Nuttle
//call_pb_int_status.c
//Call: ./call_pb_int_status gzipped_finalseqs_file plasmid_miptargets_file pb_annotation_codes_file
//
//This program analyzes finalized MIP sequence data to determine whether a given sample has one or more PB integrations
//in its genome, and if so, the type of PB integrations present (e.g., indel guide PB integrations or prime editing guide PB
//integrations). The plasmid miptargets file must contain only those MIPs targeting plasmid sequences, with annotations in the
//CRISPR column detailing whether each MIP targets sequence inside a PB transposon (and if so, which type) or outside the PB
//transposon but present in the plasmid backbone or the PBase encoding plasmid. The pb annotation codes file "pb_annotations.pbcode"
//lists different genotypes regarding PB integration status in the first column and strings of 0s and 1s in the 2nd column indicating
//among all plasmid-targeting MIPs, which must have captured and which must not have captured for the caller to assign the
//corresponding genotype to the sample being analyzed.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include<zlib.h>
#define NLEN 200 //size of character vectors for storing names, etc.
#define LLEN 1500 //maximum length of single line of text in input finalseqs file

//set up structure to store MIP target information
struct miptarg
{
	char name[NLEN+1];
	int captured;
};

//set up structure to store PB integration genotypes and MIP presence/absence patterns
struct genotype
{
	char status[NLEN+1];
	char code[NLEN+1];
};

long count_targs(FILE*mtargs);
void get_targ_info(FILE*mtargs,struct miptarg*targs);
long count_genotypes(FILE*gcodes);
void get_geno_info(FILE*gcodes,struct genotype*gtypes);
long findmip(char*mipname,struct miptarg*targs,long numtargs);
long findgeno(char*code,struct genotype*gtypes,long numstats);

int main(int argc,char*argv[])
{
	//get sample name
	char sample[NLEN+1];
	strncpy(sample,*(argv+1),NLEN);
	sample[strchr(sample,'.')-sample]='\0';
	
	//determine the number of MIP targets and allocate memory to store MIP information
	FILE*miptargs=fopen(*(argv+2),"r");
	long ntargs=count_targs(miptargs);
	struct miptarg*targets;
	targets=(struct miptarg*)malloc(ntargs*sizeof(struct miptarg));

	//read in information on MIP targets
	get_targ_info(miptargs,targets);

	//determine the number of PB integration genotypes and allocate memory to store genotype information
	FILE*pbcode=fopen(*(argv+3),"r");
	long nstats=count_genotypes(pbcode);
	struct genotype*genotypes;
	genotypes=(struct genotype*)malloc(nstats*sizeof(struct genotype));

	//read in information on PB integration genotypes
	get_geno_info(pbcode,genotypes);

	//read in final called MIP sequences and process them one by one
	gzFile*finalseqs=gzopen(*(argv+1),"r");
	char line[LLEN+1],mip[NLEN+1];
	long m;
	while(gzgets(finalseqs,line,LLEN))
	{
		sscanf(line,"%*s %s %*s %*s %*s %*s %*s %*s %*s %*s",mip);
		m=findmip(mip,targets,ntargs);
		if(m>=0)
			targets[m].captured=1;
	}

	//create string of 0s and 1s corresponding to whether or not each MIP captured (i.e., whether or not its target sequence is present)
	char mipcode[NLEN+1]="",temp[NLEN+1]="";
	for(m=0;m<ntargs;m++)
	{
		sprintf(temp,"%d",targets[m].captured);
		strncat(mipcode,temp,NLEN);
	}

	//compare observed pattern of MIP target presence/absence with patterns expected under different genotypes; call genotype if a match is found
	char genotype[NLEN+1]="";
	long g=findgeno(mipcode,genotypes,nstats);
	if(g>=0)
		strncpy(genotype,genotypes[g].status,NLEN);
	else
		strncpy(genotype,"random",NLEN);

	//print genotype call
	printf("%s\t%s\n",sample,genotype);

	//clean up and exit
	free(targets);
	free(genotypes);
	gzclose(finalseqs);
	fclose(miptargs);
	fclose(pbcode);
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

void get_targ_info(FILE*mtargs,struct miptarg*targs)
{
	long m=0;
	while(fscanf(mtargs,"%s %*s %*s %*s %*s %*s %*s %*s %*s %*s",targs[m].name)==1)
	{
		targs[m].captured=0;
		m++;
	}
	return;
}

long count_genotypes(FILE*gcodes)
{
	long numstats=0;
	char codename[NLEN+1];
	fpos_t pos;
	fgetpos(gcodes,&pos);
	while(fscanf(gcodes,"%s %*s",codename)==1)
		numstats++;
	fsetpos(gcodes,&pos);
	return numstats;
}

void get_geno_info(FILE*gcodes,struct genotype*gtypes)
{
	long g=0;
	while(fscanf(gcodes,"%s %s",gtypes[g].status,gtypes[g].code)==2)
		g++;
	return;
}

long findmip(char*mipname,struct miptarg*targs,long numtargs)
{
	long mipp;
	for(mipp=0;mipp<numtargs;mipp++)
	{
		if(strncmp(mipname,targs[mipp].name,NLEN)==0)
			return mipp;
	}
	return -1;
}

long findgeno(char*code,struct genotype*gtypes,long numstats)
{
	long geno;
	for(geno=0;geno<numstats;geno++)
	{
		if(strncmp(code,gtypes[geno].code,NLEN)==0)
			return geno;
	}
	return -1;
}

