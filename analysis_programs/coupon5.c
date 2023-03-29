//Xander Nuttle
//coupon5.c
//Call: ./coupon5 (long)number_of_coupon_draws_per_simulation_D (long)number_of_coupons_N (long)number_of_simulations_Z
//
//This program simulates outcomes from the generalized coupon's collector problem, namely, to obtain simulated counts of each coupon
//collected after D coupon draws. For PB experiments, these counts can be interpreted as simulated counts of each intrgrated gRNA construct
//after D integrations. This program specifically determines the average counts of the most abundant collected coupon to the least abundant
//collected coupon across Z simulations, such that these average counts can be compared to PB integration data.

#include<stdio.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<stdlib.h>
#include<string.h>

const gsl_rng* setup_rng(const gsl_rng*random);
int compfun(const void*p1,const void*p2);

int main(int argc,char*argv[])
{
	//get desired numbers of draws per simulation (D), coupons (N), and simulations (Z) from the command line
	long ndraws=strtol(*(argv+1),NULL,10);
	long ncoupons=strtol(*(argv+2),NULL,10);
	long nsims=strtol(*(argv+3),NULL,10);	

	//set up random number generator
	const gsl_rng*r=setup_rng(r);

	//set up arrays to store cumulative counts of most to least abundant coupons, counts for a single simulation, and average counts
	long*counts=(long*)malloc(ncoupons*sizeof(long));
	unsigned int*simcounts=(unsigned int*)malloc(ncoupons*sizeof(unsigned int));
	double*avgcounts=(double*)malloc(ncoupons*sizeof(double));

	//set up array to store probabilities of drawing each coupon during simulations
	double*probs=(double*)malloc(ncoupons*sizeof(double));
	long n;
	for(n=0;n<ncoupons;n++)
	{
		probs[n]=(double)(1.0/ncoupons);
		counts[n]=0;
		simcounts[n]=0;
		avgcounts[n]=0.0;
	}	

	//run simulations and update counts
	long z;
	for(z=0;z<nsims;z++)
	{
		gsl_ran_multinomial(r,ncoupons,ndraws,probs,simcounts);
		qsort(simcounts,ncoupons,sizeof(unsigned int),compfun);
		for(n=0;n<ncoupons;n++)
			counts[n]+=(long)simcounts[n];
	}

	//determines the average counts of the most abundant collected coupon to the least abundant collected coupon
	for(n=0;n<ncoupons;n++)
		avgcounts[n]=(double)(counts[n])/(double)nsims;

	//report results 
	printf("Rank\tAvgcount\n");
	for(n=0;n<ncoupons;n++)
		printf("%ld\t%lf\n",n+1,avgcounts[n]);

	//clean up and exit
	free(counts);
	free(simcounts);
	free(avgcounts);
	free(probs);
	return 0;
}

const gsl_rng* setup_rng(const gsl_rng*random)
{
	unsigned long seed;
	FILE*seedin;
	random=gsl_rng_alloc(gsl_rng_mt19937);
	seedin=popen("date +%N","r");
	fscanf(seedin,"%lu",&seed);
	pclose(seedin);
	gsl_rng_set(random,seed);
	return random;
}

int compfun(const void*p1,const void*p2)
{
	const unsigned int*n1=p1;
	const unsigned int*n2=p2;
	if(*n1<*n2)
		return 1;
	else
		return -1;
}

