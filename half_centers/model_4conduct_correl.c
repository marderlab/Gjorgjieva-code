/*Code to run the model neuron simulation*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
//#include </n/sw/gsl-1.15_gcc-4.7.2/include/gsl/gsl_randist.h>
//#include </n/sw/gsl-1.15_gcc-4.7.2/include/gsl/gsl_rng.h>
#include <time.h>

/*Define some parameters */
#define total_time 1000000.0   //the time of the presynaptic input in ms
#define dt 0.01            //time step in ms 
//#define N 10000

// Parameters of the model: global
//p values
int p_Na = 3;
int p_A = 3;
int p_Kd = 4;
int p_H = 1;
//reversal potentials
double E_Na = 50; //mV
double E_A = -80; //mV
double E_Kd = -80; //mV
double E_H = -20; //mV
double E_L = -50; //mV
double C = 1; //nF 
//synaptic
double E_syn = -80; //mV
double tau_syn = 100; //ms
double V_half = -45; //mV
double V_slope = 5; //mV


//int N = (int) total_time/dt;
//double stim[N];

double get_power(double x, int y)
{

  double prod = 1.0;
  for (int kk=0;kk<y;kk++)
    {
      prod = prod*x;
    }
  return prod;
}

int main(int argc, char *argv[])
{
  double mean = atof(argv[1]);
  double sigma = atof(argv[2]);
  double g_Na = atof(argv[3]); //uS/nF 
  //for cell 1
  double ratio1 = atof(argv[4]); //ratio between g_Na/g_Kd
  double ratio2 = atof(argv[5]); //ratio between g_Kd/g_A
  double ratio3 = atof(argv[6]); //ratio between g_A/g_H
  //for cell 2
  double ratio1_cell2 = atof(argv[7]); //ratio between g_Na/g_Kd
  double ratio2_cell2 = atof(argv[8]); //ratio between g_Kd/g_A
  double ratio3_cell2 = atof(argv[9]); //ratio between g_A/g_H
  //synaptic
  double g_syn1 =  atof(argv[10]);
  double g_syn2 =  atof(argv[11]);
  //
  double corr_coef = atof(argv[12]); //correlation coefficient of the input
  double tau_c = atof(argv[13]); //correlation time const
  int num_run = atof(argv[14]);

  //don't use, use voltage derivative
  double V_th = -15; //mV, threshold for detecting spikes

  /*open files where spikes will be recorded*/
  char sp_file [999];
  sprintf (sp_file, "spikes1_mean_%g_sig_%g_gNa_%g_rat1_%g_%g_rat2_%g_%g_rat3_%g_%g_gsyn_%g_%g_corr_%g_tc_%g_dt_%g_num_%d.dat", mean, sigma, g_Na, ratio1, ratio1_cell2, ratio2, ratio2_cell2, ratio3, ratio3_cell2, g_syn1, g_syn2, corr_coef,tau_c, dt, num_run);
  FILE *spfile;
  spfile = fopen(sp_file,"w");

  char sp_file2 [999];
  sprintf (sp_file2, "spikes2_mean_%g_sig_%g_gNa_%g_rat1_%g_%g_rat2_%g_%g_rat3_%g_%g_gsyn_%g_%g_corr_%g_tc_%g_dt_%g_num_%d.dat", mean, sigma, g_Na, ratio1, ratio1_cell2, ratio2, ratio2_cell2, ratio3, ratio3_cell2, g_syn1, g_syn2,  corr_coef,tau_c, dt, num_run);
  FILE *spfile2;
  spfile2 = fopen(sp_file2,"w");

  /*
  char V_file [999];
  sprintf (V_file, "V1_mean_%g_sig_%g_gNa_%g_rat1_%g_%g_rat2_%g_%g_rat3_%g_%g_gsyn_%g_%g_corr_%g_tc_%g_dt_%g.dat", mean, sigma, g_Na, ratio1, ratio1_cell2, ratio2, ratio2_cell2, ratio3, ratio3_cell2, g_syn1, g_syn2, corr_coef, tau_c, dt);
  FILE *Vfile;
  Vfile = fopen(V_file,"w");
  
  char V_file2 [999];
  sprintf (V_file2, "V2_mean_%g_sig_%g_gNa_%g_rat1_%g_%g_rat2_%g_%g_rat3_%g_%g_gsyn_%g_%g_corr_%g_tc_%g_dt_%g.dat", mean, sigma, g_Na, ratio1, ratio1_cell2, ratio2, ratio2_cell2, ratio3, ratio3_cell2,g_syn1, g_syn2,  corr_coef, tau_c, dt);
  FILE *Vfile2;
  Vfile2 = fopen(V_file2,"w");
  */

   /*For random number generation:*/
  gsl_rng * r;
  r = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(r, time(NULL));
  /*done*/

  int ii, jj, nn;

/* Defining the model
     gate = [ m h ][4]   gating variables for the four conductances
     g = [ g_Na g_K g_A g_H g_leak ]    maximal conductances
     E = [ E_Na E_K E_A E_H E_leak ]    Nernst (reversal) potentials 
     gmh = [ conductances as a function of time ]
*/

  double g_L = 0.01; //uS/nF  
  //for cell 1
  double g_Kd = g_Na*ratio1;  //uS/nF 
  double g_A = g_Na*ratio2;  //uS/nF 
  double g_H = g_Na*ratio3; //uS/nF  
  //for cell 2
  double g_Kd_cell2 = g_Na*ratio1_cell2;  //uS/nF 
  double g_A_cell2 = g_Na*ratio2_cell2;  //uS/nF 
  double g_H_cell2 = g_Na*ratio3_cell2; //uS/nF

  //printf("%g %g %g %g \n", g_Na, g_Kd, g_A, g_H);

  /*Initial conditions*/
  double V[2];
  V[0] = -70; //mV
  V[1] = -70;
  double g[4][2];
  double gate[2][4][2];
  double tau[2][4][2];
  double gate_inf[2][4][2];
  double E[4];
  int p[4];
  double gnmh[4][2];
  for (nn=0;nn<2;nn++)
    {
      for (ii=0;ii<4;ii++)
	{
	  gate[0][ii][nn] = 0.02; //m
	  gate[1][ii][nn] = 0.7; //h
	  gnmh[ii][nn] = 0.0;
	} 
      g[0][nn] = g_Na; 
    }
   E[0] = E_Na;
   E[1] = E_Kd;
   E[2] = E_A;
   E[3] = E_H;
   p[0] = p_Na;
   p[1] = p_Kd;
   p[2] = p_A;
   p[3] = p_H;   
   
   g[1][0] = g_Kd;   g[1][1] = g_Kd_cell2;
   g[2][0] = g_A;    g[2][1] = g_A_cell2;
   g[3][0] = g_H;    g[3][1] = g_H_cell2;

  double current_time = 0.0;
  double s_c = mean;
  double s_1 = mean;
  double s_2 = mean;
  double s[2]; 
  double Svar[2]; Svar[0] = 0; Svar[1] = 0;
  double S_inf[2]; S_inf[0] = 0; S_inf[1] = 0;
  double g_syn[2]; g_syn[0] = g_syn1; g_syn[1] = g_syn2;
  int counter = 0;

  int voltage_high1 = 0;
  double spike1 = 0.0; //to impose an absolute refractory period when reading spikes
  int voltage_high2 = 0;
  double spike2 = 0.0; //to impose an absolute refractory period when reading spikes
  double rand_g1 = 0.0;
  double rand_g2 = 0.0;
  double rand_g3 = 0.0;

  double ds = 0.0;
  //  double I_ion[2];   I_ion[0] = 0.0;     I_ion[1] = 0.0;
 
  double sum_g[2], sum_Eg[2], tau_V[2], V_inf[2];
  while (current_time < total_time)
   {       
     current_time += dt;
     counter ++;
     
     //common input
     rand_g1 = gsl_ran_gaussian (r,1.0);
     //s = mean + sigma/sqrt(dt)*rand_g; //if white noise
     //if correlated noise with correlation time constant tau_c
     ds = (mean-s_c)*(1-exp(-dt/tau_c)) + sigma/sqrt(dt)*sqrt(1-exp(-2*dt/tau_c))*rand_g1; 
     s_c += ds;
     //independent input cell 1
     rand_g2 = gsl_ran_gaussian (r,1.0);
     ds = (mean-s_1)*(1-exp(-dt/tau_c)) + sigma/sqrt(dt)*sqrt(1-exp(-2*dt/tau_c))*rand_g2; 
     s_1 += ds;
     //independent input cell 2
     rand_g3 = gsl_ran_gaussian (r,1.0);
     ds = (mean-s_2)*(1-exp(-dt/tau_c)) + sigma/sqrt(dt)*sqrt(1-exp(-2*dt/tau_c))*rand_g3; 
     s_2 += ds;

     s[0] = s_1*(1-corr_coef) + s_c*corr_coef; //input to cell 1
     s[1] = s_2*(1-corr_coef) + s_c*corr_coef; //input to cell 2

     /*
       Model equations:
       C dv/dt = s - I_ion
       I_ion = sum( gate.*(V-E) )
       d gate/dt = -1/tau(v).*(gate-gate_0(v))
     */

     for (nn=0;nn<2;nn++)
       {
	 tau[0][0][nn] = 1.32 - 1.26/(1.0+exp((V[nn]+120.0)/(-25.0))); //m, Na
	 tau[1][0][nn] = 0.67/(1+exp((V[nn]+62.9)/(-10.0)))*(1.5+1.0/(1.0+exp((V[nn]+34.9)/3.6))); //h, Na
	 tau[0][1][nn] = 7.2 - 6.4/(1+exp((V[nn]+28.3)/(-19.2))); //m, Kd
	 tau[1][1][nn] = 0.0; //h, Kd
	 tau[0][2][nn] = 11.6 - 10.4/(1.0+exp((V[nn]+32.9)/(-15.2))); //m, A
	 tau[1][2][nn] = 38.6 - 29.2/(1.0+exp((V[nn]+38.9)/(-26.5))); //h, A
	 tau[0][3][nn] = 272.0 + 1499.0/(1.0+exp((V[nn]+42.2)/(-8.73))); //m, H
	 tau[1][3][nn] = 0.0; //h, H
	 gate_inf[0][0][nn] = 1.0/(1.0+exp((V[nn]+25.5)/(-5.29))); //m, Na
	 gate_inf[1][0][nn] = 1.0/(1.0+exp((V[nn]+48.9)/5.18)); //h, Na
	 gate_inf[0][1][nn] = 1.0/(1.0+exp((V[nn]+12.3)/(-11.8))); //m, Kd
	 gate_inf[1][1][nn] = 0.0; //h, Kd
	 gate_inf[0][2][nn] = 1.0/(1.0+exp((V[nn]+27.2)/(-8.7))); //m, A
	 gate_inf[1][2][nn] = 1.0/(1.0+exp((V[nn]+56.9)/4.9)); //h, A
	 gate_inf[0][3][nn] = 1.0/(1.0+exp((V[nn]+70.0)/6.0)); //m, H
	 gate_inf[1][3][nn] = 0.0; //h, H
       }
     //V_pre = V[0]
     if (V[0]>V_half)
       {
	 S_inf[1] = tanh((V[0]-V_half)/V_slope);
       }
     else
       {
	 S_inf[1] = 0;
       }
     //V_pre = V[1]
     if (V[1]>V_half)
       {
	 S_inf[0] = tanh((V[1]-V_half)/V_slope);
       }
     else
       {
	 S_inf[0] = 0;
       }			

     for (nn=0;nn<2;nn++)
       {       
	 //evolve Svar
	 Svar[nn] = Svar[nn] + (1-exp(-dt/(tau_syn*(1-S_inf[nn]))))*(S_inf[nn] - Svar[nn]);

	 //evolve the conductances, first order Euler
	 for (ii=0;ii<4;ii++)
	   {
	     gate[0][ii][nn] = gate[0][ii][nn] + (1-exp(-dt/tau[0][ii][nn]))*(gate_inf[0][ii][nn] - gate[0][ii][nn]);
	     gate[1][ii][nn] = gate[1][ii][nn] + (1-exp(-dt/tau[1][ii][nn]))*(gate_inf[1][ii][nn] - gate[1][ii][nn]);
	   }
	 //reset
	 gate[1][1][nn] = 1.0; //h, Kd
	 gate[1][3][nn] = 1.0; //h, H

	 sum_g[nn] = g_L + g_syn[nn]*Svar[nn];
	 sum_Eg[nn] = g_L*E_L + g_syn[nn]*Svar[nn]*E_syn;
	 for (ii=0;ii<4;ii++)
	   {
	     gnmh[ii][nn] = g[ii][nn]*get_power(gate[0][ii][nn],p[ii])*gate[1][ii][nn];
	     sum_g[nn] += gnmh[ii][nn];
	     sum_Eg[nn] += gnmh[ii][nn]*E[ii];
	   }
	 /*
	 I_ion[nn] = g_L*(V[nn]-E_L); //sum of all ionic currents
	 for (ii=0;ii<4; ii++)
	   {
	     I_ion[nn] += gnmh[ii][nn]*(V[nn]-E[ii]);
	   }
	 */
	 tau_V[nn]  = C/sum_g[nn];  //Membrane time constant.
	 
	 // V_inf is steady-state voltage.
	 V_inf[nn] = (sum_Eg[nn] + s[nn])*tau_V[nn];
	 
	 //evolve the voltage using exponential Euler
	 V[nn] = V_inf[nn] + (V[nn]-V_inf[nn])*exp(-dt/tau_V[nn]);
       }
     
     /*
       if (fmod(current_time,0.1)<dt)
       {
	 //printf("%g \n",current_time);
	 fprintf(Vfile,"%f \n", V[0]);
	 fprintf(Vfile2,"%f \n", V[1]);
       }
     */
     
     //find out if there's a spike and record it
     if (V[0] > V_th  & current_time>5000)
       {
	 //if a high threshold was detected, record the spike
	 if (voltage_high1 == 0 && current_time-spike1>5) //5 is tau_abs
	   {
	     fprintf(spfile,"%f \n", current_time);
	     spike1 = current_time;
	   }
	 
	 //if the high threshold was continuous (after spike has occured)
	 //then wait for it to cross the threshold again before recording
	 voltage_high1 ++;
       }
     else
       {
	 if (voltage_high1 > 0)
	   {
	     voltage_high1 = 0; //reset
	   }
       }

     //same for cell 2
     if (V[1] > V_th  & current_time>5000)
       {
	 //if a high threshold was detected, record the spike
	 if (voltage_high2 == 0 && current_time-spike2>5) //5 is tau_abs
	   {
	     fprintf(spfile2,"%f \n", current_time);
	     spike2 = current_time;
	   }
	 
	 //if the high threshold was continuous (after spike has occured)
	 //then wait for it to cross the threshold again before recording
	 voltage_high2 ++;
       }
     else
       {
	 if (voltage_high2 > 0)
	   {
	     voltage_high2 = 0; //reset
	   }
       }
    
   }
     
  fclose(spfile);
  fclose(spfile2);
  //fclose(Vfile);
  //fclose(Vfile2);
 
 return 0;
}
