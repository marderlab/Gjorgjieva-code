/*Code to run the STG model neuron simulation*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
//#include </n/sw/gsl-1.15_gcc-4.7.2/include/gsl/gsl_randist.h>
//#include </n/sw/gsl-1.15_gcc-4.7.2/include/gsl/gsl_rng.h>
#include <time.h>

/*Define some parameters */
#define total_time 2005000.0   //the time of the presynaptic input in ms
#define dt 0.01            //time step in ms 
#define R_F 8.6174e-005 //needed for the computation of the revsersal for Ca
#define T 10

// Parameters of the model: global
//p values
int p_Na = 3;
int p_A = 3;
int p_K = 4;
int p_H = 1;
int p_CaT = 3;
int p_CaS = 3;
int p_KCa = 4;

double get_power(double x, int y)
{

  double prod = 1.0;
  for (int kk=0;kk<y;kk++)
    {
      prod = prod*x;
    }
  return prod;
}
double Xinf(double V0, double A, double B)
{

  double X = 1.0/(1.0+exp((V0+A)/B));
  return X;
}

double m_KCainf(double V0, double Ca0)
{
  double X = (Ca0/(Ca0+3.0))/(1.0+exp((V0+28.3)/-12.6));
  return X;
}

double tauX(double V0, double A, double B, double D, double E)
{
  double X = A - B/(1+exp((V0+D)/E));
  return X;
}

double tauhNa(double V0)
{
  double X = (0.67/(1.0+exp((V0+62.9)/-10.0)))*(1.5 + 1.0/(1.0+exp((V0+34.9)/3.6)));
  return X;
}

double taumCaS(double V0)
{
  double X = 1.4 + (7.0/((exp((V0+27.0)/10.0))+(exp((V0+70.0)/-13.0))));
  return X;
}

double tauhCaS(double V0)
{
  double X = 60.0 + (150.0/((exp((V0+55)/9.0))+(exp((V0+65.0)/-16.0))));
  return X;
}


int main(int argc, char *argv[])
{

  double g_Na = atof(argv[1]); //uS/nF 
  double g_CaT = atof(argv[2]);
  double g_CaS = atof(argv[3]);
  double g_A = atof(argv[4]);
  double g_KCa = atof(argv[5]);
  double g_K = atof(argv[6]);
  double g_H = atof(argv[7]);
  double g_L = 0.01; //uS/nF  

  //don't use, use voltage derivative
  double V_th = -15; //mV, threshold for detecting spikes

  /*open files where spikes will be recorded*/
  char sp_file [999];
  sprintf (sp_file, "sp_gNa_%g_gCaT_%g_gCaS_%g_gA_%g_gKCa_%g_gK_%g_gH_%g_gL_%g.dat",g_Na,g_CaT,g_CaS,g_A,g_KCa,g_K,g_H,g_L);
  FILE *spfile;
  spfile = fopen(sp_file,"w");
  
  char V_file [999];
  sprintf (V_file, "V_gNa_%g_gCaT_%g_gCaS_%g_gA_%g_gKCa_%g_gK_%g_gH_%g_gL_%g.dat",g_Na,g_CaT,g_CaS,g_A,g_KCa,g_K,g_H,g_L);
  FILE *Vfile;
  Vfile = fopen(V_file,"w");

   /*For random number generation:*/
  gsl_rng * r;
  r = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(r, time(NULL));
  /*done*/

  int ii, jj;

  /* Defining the model
     gate = [ m h ][4]   gating variables for the four conductances
     g = [ g_Na g_CaT g_CaS g_A g_KCa g_K g_H g_leak ]    maximal conductances
     E = [ E_Na E_CaT E_CaS E_A E_KCa E_K E_H E_leak ]    Nernst (reversal) potentials 
     gmh = [ conductances as a function of time ]
  */
  
 //reversal potentials
  double E_Na = 50; //mV or 30?
  double E_A = -80; //mV
  double E_K = -80; //mV
  double E_H = -20; //mV
  double E_Ca = 120; //mV default, but will vary
  double E_L = -50; //mV
  double C = 1; //nF 
  double tau_Ca = 20; //ms
  double I_Ca = 0;
  //synaptic
  double E_syn = -80; //mV
  double tau_syn = 100; //ms
  //Grashow V_half=-45; V_slope = 5;
  //otherwise I used:  V_half = 40; V_slope = -2.0;
  double V_half = -45.0;  //mV
  double V_slope = 5.0;  //mV

  /*Initial conditions*/
  double V = -50;
  double Ca = 0.2;
  double Ca_inf = Ca;
  double g[7];
  double gate[2][7];
  double tau[2][7];
  double gate_inf[2][7];
  double E[7];
  int p[7];
  double gnmh[7];

  gate_inf[0][0] = Xinf(V,25.5,-5.29); //m, Na
  gate_inf[1][0] = Xinf(V,48.9,5.18); //h, Na
  gate_inf[0][1] = Xinf(V,27.1,-7.2); //m, CaT
  gate_inf[1][1] = Xinf(V,32.1,5.5); //h, CaT
  gate_inf[0][2] = Xinf(V,33.0,-8.1); //m, CaS
  gate_inf[1][2] = Xinf(V,60.0,6.2); //h, CaT
  gate_inf[0][3] = Xinf(V,27.2,-8.7); //m, A
  gate_inf[1][3] = Xinf(V,56.9,4.9); //h, A
  gate_inf[0][4] = m_KCainf(V,Ca); //m, kCa 
  gate_inf[1][4] = 0.0; //h, kCa
  gate_inf[0][5] = Xinf(V,12.3,-11.8); //m, Kd
  gate_inf[1][5] = 0.0; //h, Kd	
  gate_inf[0][6] = Xinf(V,70.0,6.0); //m, H
  gate_inf[1][6] = 0.0; //h, H
  for (ii=0;ii<7;ii++)
    {
      gate[0][ii] =  gate_inf[0][ii];
      gate[1][ii] =  gate_inf[1][ii];
      gnmh[ii] = 0.0;
    } 

  E[0] = E_Na; 
  E[1] = E_Ca; 
  E[2] = E_Ca;
  E[3] = E_A;
  E[4] = E_K; 
  E[5] = E_K;
  E[6] = E_H;
  //
  p[0] = p_Na;
  p[1] = p_CaT;
  p[2] = p_CaS;
  p[3] = p_A;
  p[4] = p_KCa;
  p[5] = p_K;
  p[6] = p_H; 
  
  g[0] = g_Na;  
  g[1] = g_CaT; 
  g[2] = g_CaS; 
  g[3] = g_A;  
  g[4] = g_KCa;
  g[5] = g_K;  
  g[6] = g_H;   

  double current_time = 0.0;
  int counter = 0;

  int voltage_high = 0;
  double spike = 0.0; //to impose an absolute refractory period when reading spikes

  double sum_g = 0;
  double sum_Eg = 0;
  double tau_V = 0;
  double V_inf = 0;
  
  while (current_time < total_time)
   {       
     current_time += dt;
     counter ++;
     //printf("%f \n",current_time);

     /*
       Model equations:
       C dv/dt = s - I_ion
       I_ion = sum( gate.*(V-E) )
       d gate/dt = -1/tau(v).*(gate-gate_0(v))
     */

     E_Ca = 500.0*R_F*(T + 273.15)*log(3000.0/Ca);
     E[1] = E_Ca;
     E[2] = E_Ca;
 	 
     Ca_inf = 0.05 - 0.94*I_Ca;  
     // integrate Ca dynamics 
     Ca = Ca + (1-exp(-dt/tau_Ca))*(Ca_inf - Ca);

     tau[0][0] = tauX(V,1.32,1.26,120.0,-25); //m, Na
     tau[1][0] = tauhNa(V); //h, Na
     tau[0][1] = tauX(V,21.7,21.3,68.1,-20.5); //m, CaT
     tau[1][1] = tauX(V,105.0,89.8,55.0,-16.9); //h, CaT
     tau[0][2] = taumCaS(V); //m, CaS
     tau[1][2] = tauhCaS(V); //h, CaS
     tau[0][3] = tauX(V,11.6,10.4,32.9,-15.2); //m, A
     tau[1][3] = tauX(V,38.6,29.2,38.9,-26.5); //h, A
     tau[0][4] = tauX(V,90.3,75.1,46.0,-22.7); //m, KCa
     tau[1][4] = 0.0; //h, kCa
     tau[0][5] = tauX(V,7.2,6.4,28.3,-19.2); //m, Kd
     tau[1][5] = 0.0; //h, Kd	
     tau[0][6] = tauX(V,272.0,-1499.0,42.2,-8.73); //m, H
     tau[1][6] = 0.0; //h, H
	 
     gate_inf[0][0] = Xinf(V,25.5,-5.29); //m, Na
     gate_inf[1][0] = Xinf(V,48.9,5.18); //h, Na
     gate_inf[0][1] = Xinf(V,27.1,-7.2); //m, CaT
     gate_inf[1][1] = Xinf(V,32.1,5.5); //h, CaT
     gate_inf[0][2] = Xinf(V,33.0,-8.1); //m, CaS
     gate_inf[1][2] = Xinf(V,60.0,6.2); //h, CaT
     gate_inf[0][3] = Xinf(V,27.2,-8.7); //m, A
     gate_inf[1][3] = Xinf(V,56.9,4.9); //h, A
     gate_inf[0][4] = m_KCainf(V,Ca); //m, kCa 
     gate_inf[1][4] = 0.0; //h, kCa
     gate_inf[0][5] = Xinf(V,12.3,-11.8); //m, Kd
     gate_inf[1][5] = 0.0; //h, Kd	
     gate_inf[0][6] = Xinf(V,70.0,6.0); //m, H
     gate_inf[1][6] = 0.0; //h, H

     //evolve the conductances, first order Euler
     for (ii=0;ii<7;ii++)
       {
	 gate[0][ii] = gate[0][ii] + (1-exp(-dt/tau[0][ii]))*(gate_inf[0][ii] - gate[0][ii]);
	 gate[1][ii] = gate[1][ii] + (1-exp(-dt/tau[1][ii]))*(gate_inf[1][ii] - gate[1][ii]);
       }
     //reset
     gate[1][4] = 1.0; //h, KCa
     gate[1][5] = 1.0; //h, Kd
     gate[1][6] = 1.0; //h, H
     
     sum_g = g_L;
     sum_Eg = g_L*E_L;
     for (ii=0;ii<7;ii++)
       {
	 gnmh[ii] = g[ii]*get_power(gate[0][ii],p[ii])*gate[1][ii];
	 sum_g += gnmh[ii];
	 sum_Eg += gnmh[ii]*E[ii];
       }
     I_Ca = (gnmh[1] + gnmh[2])*(V-E_Ca);
     
     tau_V  = C/sum_g;  //Membrane time constant.
	 
     // V_inf is steady-state voltage.
     V_inf = sum_Eg/sum_g;
	 
     //evolve the voltage using exponential Euler
     V = V_inf + (V-V_inf)*exp(-dt/tau_V);
     
     //find out if there's a spike and record it
     if (V > V_th  & current_time>5000)
       {
	 //if a high threshold was detected, record the spike
	 if (voltage_high == 0 && current_time-spike>5) //5 is tau_abs
	   {
	     fprintf(spfile,"%f \n", current_time);
	     spike = current_time;
	   }
	     
	 //if the high threshold was continuous (after spike has occured)
	 //then wait for it to cross the threshold again before recording
	 voltage_high ++;
       }
     else
       {
	 if (voltage_high > 0)
	   {
	     voltage_high = 0; //reset
	   }
       }
     
     //only record a small subset of the voltages
     if (current_time > (total_time-5000)  && fmod(current_time,0.1)<dt)
       {
	 fprintf(Vfile,"%f \n",V);
       }
   }
     
  fclose(spfile);
  fclose(Vfile);
 
 return 0;
}
