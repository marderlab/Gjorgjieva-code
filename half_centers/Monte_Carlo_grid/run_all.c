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
  double g_Na_min = 800; double g_Na_max = 1200; //uS/nF 
  double g_CaT_min = 0; double g_CaT_max = 6;
  double g_CaS_min = 0; double g_CaS_max = 12;
  double g_A_min = 20; double g_A_max = 130;
  double g_KCa_min = 20; double g_KCa_max = 500;
  double g_Kd_min = 90; double g_Kd_max = 120;
  double g_H_min = 1; double g_H_max = 5;

  double g_Na, g_CaT, g_CaS, g_A, g_KCa, g_Kd, g_H;
  char cmd[999];

  int ii,jj,kk,tt,ll;
  
  //ll = system ("g++ -x c++ -I/usr/local/include -L/usr/local/lib model_7conduct_V_sp.c -o m.out -lgsl -lgslcblas");
  ll = system ("g++ model_7conduct_V_sp.c -o m.out -lgsl -lgslcblas");

  /*For random number generation:*/
  gsl_rng * r;
  r = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(r, time(NULL));
  /* done */

  double rand_num = 0;  
  int Num_sims = 2000;

  //open file to record model parameters
  char sp_file [999];
  sprintf (sp_file, "conductances_2_%d.dat",Num_sims);
  FILE *spfile;
  spfile = fopen(sp_file,"w");
  
  for (ii=0; ii<Num_sims; ii++)
    {
      //generate random conductances for the model
      rand_num = gsl_rng_uniform(r);
      g_Na = rand_num*(g_Na_max-g_Na_min) + g_Na_min;
      rand_num = gsl_rng_uniform(r);
      g_CaT = rand_num*(g_CaT_max-g_CaT_min) + g_CaT_min;
      rand_num = gsl_rng_uniform(r);
      g_CaS = rand_num*(g_CaS_max-g_CaS_min) + g_CaS_min;
      rand_num = gsl_rng_uniform(r);
      g_A = rand_num*(g_A_max-g_A_min) + g_A_min;
      rand_num = gsl_rng_uniform(r);
      g_KCa = rand_num*(g_KCa_max-g_KCa_min) + g_KCa_min;
      rand_num = gsl_rng_uniform(r);
      g_Kd = rand_num*(g_Kd_max-g_Kd_min) + g_Kd_min;
      rand_num = gsl_rng_uniform(r);
      g_H = rand_num*(g_H_max-g_H_min) + g_H_min;

      //record the conductances
      fprintf(spfile,"%d %g %g %g %g %g %g %g \n",ii,g_Na,g_CaT,g_CaS,g_A,g_KCa,g_Kd,g_H);
		 
      sprintf(cmd, "qsub -q all.q -cwd -o /dev/null -j y -ckpt reloc run_file.txt %f %f %f %f %f %f %f",g_Na,g_CaT,g_CaS,g_A,g_KCa,g_Kd,g_H);
      //sprintf(cmd, "./m.out %f %f %f %f %f %f %f",g_Na,g_CaT,g_CaS,g_A,g_KCa,g_Kd,g_H); 
      ll = system(cmd);
      
    }

  fclose(spfile);
  return 0;	      
}
