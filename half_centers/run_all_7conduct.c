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
  double g_syn = 0.0;
  
  char cmd[999];

  int ii,jj,kk,tt,ll;
  
  //ll = system ("g++ model_7conduct_V.c -o m.out -lgsl");   
  
  double tau_c = 1.0;
  
  while (mean < 10.0001)  //run through mean current, in units of nA
    {
      g_syn = 0.0;
      while (g_syn < 0.020001) //run through g_syn
	{
	  printf("mean %g g_syn %g \n", mean, g_syn);
	  
	  //sprintf(cmd, "sbatch run_script.txt %f %f %f %f %f %f %f", mean, sigma, gNa[ii], ratio1[jj], ratio2[kk], ratio3[tt], tau_c);
	  
	  sprintf(cmd, "./m.out %f 0 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 158.6874 2.0942 2.3907 21.32181 84.7673 52.9122 0.8532 0.0561 %f %f 0 1 2",mean, g_syn, g_syn);
	 
	  ll = system(cmd);
	  g_syn += 0.001;		  		
	} //end of g_syn loop
      mean += 0.1;      
    } //end of loop for mean    
  
  return 0;	      
}
