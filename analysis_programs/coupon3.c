//Xander Nuttle
//coupon3.c
//Call: ./coupon3 (long)number_of_sets_to_complete_M (long)number_of_coupons_N (long)number_of_simulations_Z
//
//This program implements a Markov chain in order to simulate outcomes from the generalized coupon's collector problem,
//namely, to simulate the number of coupons drawn to obtain M complete sets of N coupons. For PB experiments, this
//same number can be interpreted as the number of PB integrations needed to obtain at least M integrations of each of
//N gRNA constructs.

#include<stdio.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<stdlib.h>
#include<string.h>

const gsl_rng* setup_rng(const gsl_rng*random);
long simcc(long numsets,long numcoupons,const gsl_rng*random);
void init_state(long*svec,long n_sets,long n_coupons);
void setpaths(long n_sets,long n_coupons,long st[n_sets+1],long newsts[n_sets+1][n_sets+1],double p[n_sets+1]);
void update_state(long n_sets,long st[n_sets+1],long newsts[n_sets+1][n_sets+1],unsigned int ivec[n_sets+1]);
double avgsims(long n_sims,long simsvec[n_sims]);

int main(int argc,char*argv[])
{
	//get desired numbers of sets (M), coupons (N), and simulations (Z) from the command line
	long nsets=strtol(*(argv+1),NULL,10);
	long ncoupons=strtol(*(argv+2),NULL,10);
	long nsims=strtol(*(argv+3),NULL,10);

	//set up random number generator
	const gsl_rng*r=setup_rng(r);	

	//set up array to store simulation outcomes and run simulations
	long*sims=(long*)malloc(nsims*sizeof(long));
	long z;
	for(z=0;z<nsims;z++)
		sims[z]=simcc(nsets,ncoupons,r);

	//calculate average of simulation outsomes and report results
	double avgdraws=avgsims(nsims,sims);	
	printf("number of sets: %ld\nnumber of coupons: %ld\nnumber of simulations: %ld\naverage number of draws needed to complete sets: %lf\n",nsets,ncoupons,nsims,avgdraws);

	//clean up and exit
	free(sims);
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

long simcc(long numsets,long numcoupons,const gsl_rng*random)
{
	long ndraws=0;
	long state[numsets+1];
	long newstates[numsets+1][numsets+1];
	double probs[numsets+1];
	unsigned int newst_indx[numsets+1];
	init_state(state,numsets,numcoupons);
	while(state[numsets]<numcoupons) //all sets completed when last set has all coupons
	{
		setpaths(numsets,numcoupons,state,newstates,probs);
		gsl_ran_multinomial(random,numsets+1,1,probs,newst_indx);
		update_state(numsets,state,newstates,newst_indx);
		ndraws++;
	}
	return ndraws;
}

void init_state(long*svec,long n_sets,long n_coupons)
{
	long m;
	for(m=0;m<(n_sets+1);m++)
		svec[m]=0;
	svec[0]=n_coupons;
	return;
}

void setpaths(long n_sets,long n_coupons,long st[n_sets+1],long newsts[n_sets+1][n_sets+1],double p[n_sets+1])
{
	long s,m;
	s=0; //new state 0: remain in same state
	for(m=0;m<(n_sets+1);m++)
		newsts[s][m]=st[m];
	p[s]=(double)st[n_sets]/n_coupons;
	for(s=1;s<(n_sets+1);s++) //new states 1-M: transition to a new state with set s incremented by 1
	{
		for(m=0;m<(n_sets+1);m++)
			newsts[s][m]=st[m];
		newsts[s][s]++;
		p[s]=(double)(st[s-1]-st[s])/n_coupons;
	}
	return;
}

void update_state(long n_sets,long st[n_sets+1],long newsts[n_sets+1][n_sets+1],unsigned int ivec[n_sets+1])
{
	long s,m;
	for(s=0;s<(n_sets+1);s++)
	{
		if(ivec[s]==1)
		{
			for(m=0;m<(n_sets+1);m++)
				st[m]=newsts[s][m];
			break;
		}
	}
	return;
}

double avgsims(long n_sims,long simsvec[n_sims])
{
	long long draws=0;
	long s;
	for(s=0;s<n_sims;s++)		
		draws+=simsvec[s];
  return (double)draws/(double)n_sims;
}

