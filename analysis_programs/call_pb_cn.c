//Xander Nuttle
//call_pb_cn.c
//Call: ./call_pb_cn guidecounts_file
//
//Takes a file containing counts of molecular tags corresponding to U6 and H1 MIP capture events for each
//guide RNA and each guide RNA pair and prints to standard output the sample name followed by the number
//of piggyBac (PB) insertions for that sample, inferred from the U6 MIP capture event counts.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#define CMIN 20 //calling minimum; if the highest U6 count is below this value, PB copy number will be called as uncertain but nonzero (-1) or zero
#define CNZO 5 //count minimum to be considered nonzero; if the highest U6 count is below this value, PB copy number will be called as zero
#define FRAC 0.1 //fraction to multiply highest U6 count by in PB copy number calling; other counts below FRAC*highest_U6_count will be considered as zero
#define CCUT 5 //count cutoff in PB copy number calling; counts (other than the highest U6 count) below this cutoff will be considered as zero 

void count_guides(FILE*gcounts,long*nguides);
int compfun(const void*p1,const void*p2);


int main(int argc,char*argv[])
{
	//get sample name
	char sample[101];
	strncpy(sample,*(argv+1),88);
	sample[strrchr(sample,'.')-sample]='\0';

	//determine number of guides and/or guide pairs to be analyzed
	FILE*guidecounts=fopen(*(argv+1),"r");
  long numguides=-1;
  count_guides(guidecounts,&numguides);

	//initialize array of U6 MIP capture event counts
	long counts[numguides];
	long g;
	for(g=0;g<numguides;g++)
		fscanf(guidecounts,"%*s %ld %*s",&(counts[g]));

	//sort array of U6 MIP capture event counts
	qsort(counts,numguides,sizeof(long),compfun);

	//call PB copy number
	long pbcn=1;
	double mincount=FRAC*counts[0];
	if(counts[0]<CNZO)
		pbcn=0;
	else if(counts[0]<CMIN)
		pbcn=-1;
	else
	{
		for(g=1;g<numguides;g++)
		{
			if((counts[g]>=CCUT)&&(counts[g]>=mincount))
				pbcn++;
			else
				break;
		}
	}

	//print sample and corresponding called PB copy number to standard output
	printf("%s\t%ld\n",sample,pbcn);	
	
	//clean up and exit
	return 0;
}

void count_guides(FILE*gcounts,long*nguides)
{
  char dummy[51];
  fpos_t pos;
  fgetpos(gcounts,&pos);
  while(fscanf(gcounts,"%s %*s %*s",dummy)==1)
    (*nguides)++;
  fsetpos(gcounts,&pos);
	fscanf(gcounts,"%s %*s %*s",dummy);
  return;
}

int compfun(const void*p1,const void*p2)
{
	const long*n1=p1;
	const long*n2=p2;
	if(*n1<*n2)
		return 1;
	else
		return -1;
}


