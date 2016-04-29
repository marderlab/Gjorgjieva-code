#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
//#include </n/sw/gsl-1.15_gcc-4.7.2/include/gsl/gsl_randist.h>
//#include </n/sw/gsl-1.15_gcc-4.7.2/include/gsl/gsl_rng.h>
#include <time.h>
#include <string>
using namespace std;


int main ()
{
  double mean = 0.0; 
  double g_syn = 0.6;
  double tau_c = 1.0;
  double gA = 0;
  
  char cmd[999];

  int ii,jj,kk,tt,ll;
  
  //ll = system ("g++ model_7conduct_V.c -o m.out -lgsl");   
  //while (g_syn < 0.5001) //run through g_syn
    {
      mean = 0.0;
      while (mean < 10.0001)  //run through mean current, in units of nA
	{
	  gA = 0.0;	  
	  while (gA < 100.0001)
	    {
	      printf("mean=%g gA=%g gsyn=%g \n",mean, gA,g_syn);
	      
	      //sprintf(cmd, "sbatch run_script.txt %f %f %f %f %f %f %f", mean, sigma, gNa[ii], ratio1[jj], ratio2[kk], ratio3[tt], tau_c);

	       //NOTE: if num_run=0; then use model with regular synapses ./m.out
	      //if num_run=1; then Grashow synapse model ./m2.out
	      sprintf(cmd, "./m.out %f 0 227.052 2.7781 3.3469 %f 121.124 75.649 1.3194 0.1631 %f 0 0 0",mean, gA, g_syn);
	      ll = system(cmd);
	      
	      sprintf(cmd, "./m2.out %f 0 227.052 2.7781 3.3469 %f 121.124 75.649 1.3194 0.1631 %f 0 0 1",mean, gA, g_syn);	      
	      ll = system(cmd);
	      
	      gA += 10;
	    }
	  mean += 0.5;      
	} //end of loop for mean
      //if (g_syn < 0.010000001)
	{
	  //g_syn += 0.001;
	}
	//else
	{
	  g_syn += 0.1;
	}
    } //end of gA (g_syn) loop
  
  return 0;	      
}
