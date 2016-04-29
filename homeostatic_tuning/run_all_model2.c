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
  double tau_c = 1.0;
  double g_CaT = 5.5;
  double g_CaS = 0;
  double g_KCa = 0;

  double eps = 0.000001;
  
  char cmd[999];

  int ii,jj,kk,tt,ll;
  
  //ll = system ("g++ model_7conduct_V.c -o m.out -lgsl");   
  while (g_CaT < 7+eps) //run through g_syn
    {
      g_CaS = 0.0;
      while (g_CaS < 7+eps)  //run through mean current, in units of nA
	{
	  g_KCa = 0.0;	  
	  while (g_KCa < 1500+eps)
	    {
	      printf("g_CaT=%g, g_CaS=%g, g_KCa=%g, gsyn=%g \n",g_CaT,g_CaS,g_KCa,g_syn);
	      
	      //sprintf(cmd, "sbatch run_script.txt %f %f %f %f %f %f %f", mean, sigma, gNa[ii], ratio1[jj], ratio2[kk], ratio3[tt], tau_c);

	       //NOTE: if num_run=0; then use model with regular synapses ./m.out
	      //if num_run=1; then Grashow synapse model ./m2.out
	      /*
	      sprintf(cmd, "./m.out %f 0 227.052 2.7781 3.3469 %f 121.124 75.649 1.3194 0.1631 %f 0 0 0",mean, gA, g_syn);
	      ll = system(cmd);
	      */
	      
	      sprintf(cmd, "./m2.out 0 0 227.052 %f %f 30.4321 %f 75.649 1.3194 0.1631 %f 0 0 1",g_CaT,g_CaS,g_KCa,g_syn);	      
	      ll = system(cmd);

	      g_KCa += 100; 
	      
	    } //end ofg_KCa
	  g_CaS += 0.5;      
	} //end of loop for g_CaS
      g_CaT += 0.5;
    } //end of g_CaT
  
  return 0;	      
}
