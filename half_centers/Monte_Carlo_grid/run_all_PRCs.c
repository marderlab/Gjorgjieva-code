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

  int cell1 = 312;
  int cell2 = 1125;
  
  char cmd[999];
  int ll;
  //ll = system ("g++ half_centers_PRC.c -o hc.out -lgsl -lgslcblas");
  ll = system ("g++ -x c++ -I/usr/local/include -L/usr/local/lib half_centers_PRC.c -o hc.out -lgsl -lgslcblas");

  
  /*For random number generation:*/
  gsl_rng * r;
  r = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(r, time(NULL));
  /* done */

  double rand_num = 0;  

  /*
  double gmax_vec[5]; //strength of synaptic conductance
  gmax_vec[0] = 1; gmax_vec[1] = 10; gmax_vec[2] = 50; gmax_vec[3] = 100; gmax_vec[4] = 500; */
  double gmax_vec[3]; //strength of synaptic conductance
  gmax_vec[0] = 0.01; gmax_vec[1] = 0.05; gmax_vec[2] = 0.2; 
  double Tdur_vec[5]; //duration as percent of period
  Tdur_vec[0] = 0.05; Tdur_vec[1] = 0.1; Tdur_vec[2] = 0.2; Tdur_vec[3] = 0.3; Tdur_vec[4] = 0.4; 
  double Vsyn_vec[2];
  Vsyn_vec[0] = -78; Vsyn_vec[1] = 0;
 
  double g_syn1 = 0.07;
  double g_syn2 = 0.09;

   
  double phase;
  int ff;
  
  double gmax, Tdur, Vsyn;
  for (int gg=0;gg<3;gg++) //different amplitude of syn input
    {
      gmax = gmax_vec[gg];
      for (int tt=0;tt<5;tt++) //different duration of syn input
	{
	  Tdur = Tdur_vec[tt];
	  {
	    for (int vv=0;vv<2;vv++) //exc or inh input
	      {
		Vsyn = Vsyn_vec[vv];

		for (int ff=0; ff<19; ff++) //different phases for PRC
		  {
		    phase = (((double)ff+1.0)/20.0);

		    //for (int cell=0;cell<2;cell++) //which cell to provide input to...
		    int cell = 2;
		      {
			//for gsyn non-0
			printf("%d %d %f %f %f %f %f %d %f %f %f %f %f \n",cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur,cell,phase,2040.3,2134.5,179.2,179.2);
			sprintf(cmd, "./hc.out %d %d %f %f %f %f %f %d %f %f %f %f %f",cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur,cell,phase,2040.3,2134.5,179.2,179.2);

			// for gsyn = 0
			//printf("%d %d %f %f %f %f %f %d %f %f %f %f %f \n",cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur,cell,phase,2045.6,2071.3,154.6,143.0);
			//sprintf(cmd, "./hc.out %d %d %f %f %f %f %f %d %f %f %f %f %f",cell1,cell2,g_syn1,g_syn2,Vsyn,gmax,Tdur,cell,phase,2045.6,2071.3,154.6,143.0);
			ll = system(cmd);
		      }
		  }
	      }
	  }
	}
    }
  
  return 0;	      
}
