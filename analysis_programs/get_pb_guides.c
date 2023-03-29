//Xander Nuttle
//get_pb_guides.c
//Call: ./get_pb_guides guidecounts_file
//
//Takes a file containing counts of molecular tags corresponding to U6 and H1 MIP capture events for each
//guide RNA and each guide RNA pair and prints to standard output the sample name followed by the names of
//all piggyBac (PB) insertions for that sample, inferred from the U6 MIP capture event counts.

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#define CMIN 20 //calling minimum; if the highest U6 count is below this value, PB copy number will be called as uncertain but nonzero (-1) or zero
#define CNZO 5 //count minimum to be considered nonzero; if the highest U6 count is below this value, PB copy number will be called as zero
#define FRAC 0.1 //fraction to multiply highest U6 count by in PB copy number calling; other counts below FRAC*highest_U6_count will be considered as zero
#define CCUT 5 //count cutoff in PB copy number calling; counts (other than the highest U6 count) below this cutoff will be considered as zero 
#define NONE "none" //placeholder string for samples determined to have zero PB integrations and this zero guide RNAs/guide RNA pairs integrated
#define LOWC "lowcounts" //placeholder string for samples called as having an uncertain number of PB integrations due to low counts

//set up guidecount structure
struct gcounts
{
	char gname[51];
	long count;
};

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

	//initialize guidecount structure
	struct gcounts*guides;	
	guides=(struct gcounts*)malloc(numguides*sizeof(struct gcounts));
	long g;
	for(g=0;g<numguides;g++)
		fscanf(guidecounts,"%s %ld %*s",guides[g].gname,&(guides[g].count));

	//sort array of U6 MIP capture event counts
	qsort(guides,numguides,sizeof(struct gcounts),compfun);

	//call PB copy number
	long pbcn=1;
	double mincount=FRAC*guides[0].count;
	if(guides[0].count<CNZO)
		printf("%s\t%s\n",sample,NONE);		
	else if(guides[0].count<CMIN)
		printf("%s\t%s\n",sample,LOWC);		
	else
	{
		for(g=0;g<numguides;g++)
		{
			if((guides[g].count>=CCUT)&&(guides[g].count>=mincount))
				printf("%s\t%s\n",sample,guides[g].gname);
			else
				break;
		}
	}

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
	const struct gcounts*n1=p1;
	const struct gcounts*n2=p2;
	if((n1->count)<(n2->count))
		return 1;
	else
		return -1;
}

