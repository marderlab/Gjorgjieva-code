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

  int N = 42; //numer of pairs to test

  //specify these... for the pair
  int cell1[N];
  int cell2[N];              
                                                           
  cell1[0] = 96;   cell1[1] = 76;   cell1[2] = 163;  cell1[3] = 19;   cell1[4] = 133;
  cell1[5] = 297;  cell1[6] = 64;   cell1[7] = 276;  cell1[8] = 363;  cell1[9] = 241;
  cell1[10] = 146; cell1[11] = 358; cell1[12] = 470; cell1[13] = 157; cell1[14] = 377;
  cell1[15] = 50;  cell1[16] = 314; cell1[17] = 466; cell1[18] = 137; cell1[19] = 96;
  cell1[20] = 287; cell1[21] = 285; cell1[22] = 247; cell1[23] = 221; cell1[24] = 198;
  cell1[25] = 181; cell1[26] = 181; cell1[27] = 77;  cell1[28] = 77;  cell1[29] = 77; 
  cell1[30] = 61;  cell1[31] = 61;  cell1[32] = 61;  cell1[33] = 52;  cell1[34] = 327; 
  cell1[35] = 283; cell1[36] = 283; cell1[37] = 357; cell1[38] = 145; cell1[39] = 163;
  cell1[40] = 430; cell1[41] = 339;
  
  cell2[0] = 1116;  cell2[1] = 1780;  cell2[2] = 1968;  cell2[3] = 200;   cell2[4] = 1179;
  cell2[5] = 924;   cell2[6] = 1346;  cell2[7] = 1249;  cell2[8] = 424;   cell2[9] = 875;
  cell2[10] = 392;  cell2[11] = 1265; cell2[12] = 872;  cell2[13] = 1550; cell2[14] = 766;
  cell2[15] = 1936; cell2[16] = 773;  cell2[17] = 509;  cell2[18] = 1821; cell2[19] = 727;
  cell2[20] = 1989; cell2[21] = 1509; cell2[22] = 1445; cell2[23] = 1939; cell2[24] = 1092;
  cell2[25] = 1653; cell2[26] = 1219; cell2[27] = 1918; cell2[28] = 1876; cell2[29] = 1796;
  cell2[30] = 973;  cell2[31] = 567;  cell2[32] = 320;  cell2[33] = 995;  cell2[34] = 644;
  cell2[35] = 661;  cell2[36] = 471;  cell2[37] = 1658; cell2[38] = 389;  cell2[39] = 1940;
  cell2[40] = 583;  cell2[41] = 1902;
  
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
  
  double del_g = 0.02;
  double g_syn1;
  double g_syn2;
  
  for (int tt=18;tt<24;tt++)
    {
      //printf("%d %d \n",cell1[tt],cell2[tt]);
      //although examine a grid of synapses, make then not exactly idential but within 5% of the chosen value - don't do it
      for (int g1=1; g1<25; g1++)
	{
	  g_syn1 = del_g*g1;// * (1 + 0.05*gsl_rng_uniform(r));
	  for (int g2=1; g2<25; g2++)
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
