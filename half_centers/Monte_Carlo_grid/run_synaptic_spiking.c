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

  //specify these... for the pair
  int cell1[17];
  int cell2[17];              
                                                           
  cell1[0] = 76; cell2[0] = 1301;
  cell1[1] = 76; cell2[1] = 1692;
  cell1[2] = 99; cell2[2] = 1080;
  cell1[3] = 99; cell2[3] = 1376;
  cell1[4] = 99; cell2[4] = 1589;
  cell1[5] = 883; cell2[5] = 1376;
  cell1[6] = 883; cell2[6] = 1526;
  cell1[7] = 883; cell2[7] = 1589;
  cell1[8] = 1301; cell2[8] = 1376;
  cell1[9] = 1301; cell2[9] = 1526;
  cell1[10] = 1301; cell2[10] = 1589;
  cell1[11] = 1376; cell2[11] = 1692;
  cell1[12] = 1376; cell2[12] = 1846;
  cell1[13] = 1526; cell2[13] = 1692;
  cell1[14] = 1526; cell2[14] = 1846;
  cell1[15] = 1589; cell2[15] = 1692;
  cell1[16] = 1589; cell2[16] = 1846;
  //cell1[17] = 1753; cell2[17] = 1919;
  //cell1[18] = 1864; cell2[18] = 1873; 

  
  char cmd[999];
  int ll;
  //ll = system ("g++ -x c++ -I/usr/local/include -L/usr/local/lib half_center_7conduct_V_sp.c -o h.out -lgsl -lgslcblas");
  ll = system ("g++ half_center_7conduct_V_sp.c -o h.out -lgsl -lgslcblas");
  
  /*For random number generation:*/
  gsl_rng * r;
  r = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(r, time(NULL));
  /* done */

  double rand_num = 0;  
  
  double V_half = 45; 
  double tau_syn = 100; //ms
  double E_syn = -78; //mV
  
  double del_g = 0.01;
  double g_syn1;
  double g_syn2;

  for (int tt=0;tt<17;tt++)
    {
      //printf("%d %d \n",cell1[tt],cell2[tt]);
      //although examine a grid of synapses, make then not exactly idential but within 5% of the chosen value - don't do it
      for (int g1=1; g1<14; g1++)
	{
	  g_syn1 = del_g*g1;// * (1 + 0.05*gsl_rng_uniform(r));
	  for (int g2=1; g2<14; g2++)
	    {
	      g_syn2 = del_g*g2;// * (1 + 0.05*gsl_rng_uniform(r));
	      
	      //record the parameters
	      //fprintf(spfile,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n",g_Na1,g_Na2,g_CaT1,g_CaT2,g_CaS1,g_CaS2,g_A1,g_A2,g_KCa1,g_KCa2,g_K1,g_K2,g_H1,g_H2,g_syn1,g_syn2,V_half,tau_syn,E_syn);
	      
	      sprintf(cmd, "qsub -q all.q -cwd -o /dev/null -j y -ckpt reloc pair_file.txt %d %d %f %f %f %f %f",cell1[tt]-1,cell2[tt]-1,g_syn1,g_syn2,V_half,tau_syn,E_syn);
	      //sprintf(cmd, "./h.out %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",g_Na1,g_Na2,g_CaT1,g_CaT2,g_CaS1,g_CaS2,g_A1,g_A2,g_KCa1,g_KCa2,g_K1,g_K2,g_H1,g_H2,g_syn1,g_syn2,V_half,tau_syn,E_syn);
	      ll = system(cmd);
	    }      
	}
    }
  
  //fclose(spfile);
  return 0;	      
}
