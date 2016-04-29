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
#define total_time 10000.0   //the time of the presynaptic input in ms
#define dt 0.001            //time step in ms 
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

void read_stimulus(FILE *stimfile, double g_vec[8][2000])
{ 
  /*Read in stimulus times from a file provided.*/
 double x = 1.0;

 for (int ii=0;ii<2000;ii++)
   {
     for (int jj=0;jj<8;jj++)
       {
	 fscanf(stimfile,"%lf \t",&x);
	 g_vec[jj][ii] = x;
       }
   }
 
 return;
}

int main(int argc, char *argv[])
{
  double mean = 0;//atof(argv[1]);
  double sigma = 0;//atof(argv[2]);
  //for cell 1
  int in1 = atof(argv[1]); //uS/nF 
  double g_L1 = 0.01; //uS/nF  
  //for cell 2
  int in2 = atof(argv[2]);
  double g_L2 = 0.01; //uS/nF
  //synaptic
  double g_syn1 =  atof(argv[3]);
  double g_syn2 =  atof(argv[4]);
  //
  //double corr_coef = atof(argv[21]); //correlation coefficient of the input
  //double tau_c = atof(argv[22]); //correlation time const
  //properties of the synapse
  double V_half =  45;  //atof(argv[5]); //half-activation, determines if synapses are graded, or spike-mediated
  double tau_syn = 100; //atof(argv[6]); //synaptic time constant
  double E_syn = -78;   //atof(argv[7]); //reversal potential, determines if synapses are inhibitory or excitatory
  
  double g_Na1, g_Na2, g_CaT1, g_CaT2, g_CaS1, g_CaS2, g_A1, g_A2, g_KCa1, g_KCa2, g_K1, g_K2, g_H1, g_H2;
  //load the conductances for each cell
  double g_vec[8][2000];
  char stim_file [999];
  FILE *stimfile;
  sprintf (stim_file, "conductances_1_2000.dat");
  stimfile = fopen(stim_file,"r");
  read_stimulus(stimfile,g_vec);
  fclose(stimfile);

  g_Na1 = g_vec[1][in1]; g_CaT1 = g_vec[2][in1]; g_CaS1 = g_vec[3][in1];
  g_A1 = g_vec[4][in1]; g_KCa1 = g_vec[5][in1]; g_K1 = g_vec[6][in1]; g_H1 = g_vec[7][in1];
  
  g_Na2 = g_vec[1][in2]; g_CaT2 = g_vec[2][in2]; g_CaS2 = g_vec[3][in2];
  g_A2 = g_vec[4][in2]; g_KCa2 = g_vec[5][in2]; g_K2 = g_vec[6][in2]; g_H2 = g_vec[7][in2];
  
  //don't use, use voltage derivative
  double V_th = -15; //mV, threshold for detecting spikes

  
  /*open files where voltage will be recorded*/
  
  char V_file [999];
  sprintf (V_file, "V_in1_%d_in2_%d_gsyn_%g_%g.dat",in1,in2,g_syn1,g_syn2);

  FILE *Vfile;
  Vfile = fopen(V_file,"w");

   /*For random number generation:*/
  gsl_rng * r;
  r = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(r, time(NULL));
  /*done*/

  int ii, jj, nn;

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
  double E_Ca[2]; E_Ca[0] = 120; E_Ca[1] = 120; //mV default, but will vary
  double E_L = -50; //mV
  double C = 1; //nF 
  double tau_Ca = 20; //ms
  double I_Ca[2]; I_Ca[0] = 0; I_Ca[1] = 0;
  //synaptic
  //double E_syn = -80; //mV
  //double tau_syn = 100; //ms
  //Grashow V_half=-45; V_slope = 5;
  //otherwise I used: 
  //double V_half = 45.0;  //mV
  double V_slope = -2.0;  //mV

  /*Initial conditions*/
  double V[2];
  V[0] = -50; //mV
  V[1] = -50;
  double Ca[2]; Ca[0] = 0.2; Ca[1] = 0.2;
  double Ca_inf[2]; Ca_inf[0] = Ca[0]; Ca_inf[1] = Ca[1];
  double g[7][2];
  double gate[2][7][2];
  double tau[2][7][2];
  double gate_inf[2][7][2];
  double E[7][2];
  int p[7];
  double gnmh[7][2];

  for (nn=0;nn<2;nn++)
    {
      gate_inf[0][0][nn] = Xinf(V[nn],25.5,-5.29); //m, Na
      gate_inf[1][0][nn] = Xinf(V[nn],48.9,5.18); //h, Na
      gate_inf[0][1][nn] = Xinf(V[nn],27.1,-7.2); //m, CaT
      gate_inf[1][1][nn] = Xinf(V[nn],32.1,5.5); //h, CaT
      gate_inf[0][2][nn] = Xinf(V[nn],33.0,-8.1); //m, CaS
      gate_inf[1][2][nn] = Xinf(V[nn],60.0,6.2); //h, CaT
      gate_inf[0][3][nn] = Xinf(V[nn],27.2,-8.7); //m, A
      gate_inf[1][3][nn] = Xinf(V[nn],56.9,4.9); //h, A
      gate_inf[0][4][nn] = m_KCainf(V[nn],Ca[nn]); //m, kCa 
      gate_inf[1][4][nn] = 0.0; //h, kCa
      gate_inf[0][5][nn] = Xinf(V[nn],12.3,-11.8); //m, Kd
      gate_inf[1][5][nn] = 0.0; //h, Kd	
      gate_inf[0][6][nn] = Xinf(V[nn],70.0,6.0); //m, H
      gate_inf[1][6][nn] = 0.0; //h, H
      for (ii=0;ii<7;ii++)
	{
	  gate[0][ii][nn] =  gate_inf[0][ii][nn];
	  gate[1][ii][nn] =  gate_inf[1][ii][nn];
	  gnmh[ii][nn] = 0.0;
	} 
    }
  

  E[0][0] = E_Na; E[0][1] = E_Na;
  E[1][0] = E_Ca[0]; E[1][1] = E_Ca[1]; 
  E[2][0] = E_Ca[0]; E[2][1] = E_Ca[1];
  E[3][0] = E_A;  E[3][1] = E_A;
  E[4][0] = E_K;  E[4][1] = E_K; 
  E[5][0] = E_K;  E[5][1] = E_K;
  E[6][0] = E_H;  E[6][1] = E_H;
  //
  p[0] = p_Na;
  p[1] = p_CaT;
  p[2] = p_CaS;
  p[3] = p_A;
  p[4] = p_KCa;
  p[5] = p_K;
  p[6] = p_H; 
  
  g[0][0] = g_Na1;   g[0][1] = g_Na2;
  g[1][0] = g_CaT1;  g[1][1] = g_CaT2;
  g[2][0] = g_CaS1;  g[2][1] = g_CaS2;
  g[3][0] = g_A1;    g[3][1] = g_A2;
  g[4][0] = g_KCa1;  g[4][1] = g_KCa2;
  g[5][0] = g_K1;    g[5][1] = g_K2;
  g[6][0] = g_H1;    g[6][1] = g_H2;

  double g_L[2]; g_L[0] = g_L1; g_L[1] = g_L2;
  double current_time = 0.0;
  double s_c = mean;
  double s_1 = mean;
  double s_2 = mean;
  double s[2]; s[0] = 0; s[1] = 0;
  double Svar[2];
  Svar[1] = Xinf(V[0],V_half,V_slope);
  Svar[0] = Xinf(V[1],V_half,V_slope); 
  double S_inf[2];
  S_inf[0] = Svar[0]; 
  S_inf[1] = Svar[1];
  double g_syn[2]; g_syn[0] = g_syn1; g_syn[1] = g_syn2;
  int counter = 0;

  int voltage_high[2]; voltage_high[0] = 0; voltage_high[1] = 0;
  double rand_g1 = 0.0;
  double rand_g2 = 0.0;
  double rand_g3 = 0.0;

  double ds = 0.0;
 
  double sum_g[2], sum_Eg[2], tau_V[2], V_inf[2];
  while (current_time < total_time)
   {       
     current_time += dt;
     counter ++;
     //printf("%f \n",current_time);
     //s = mean + sigma/sqrt(dt)*rand_g; //if white noise

     /*
     //common input
     rand_g1 = gsl_ran_gaussian (r,1.0);
     //if correlated noise with correlation time constant tau_c
     ds = (mean-s_c)*(1-exp(-dt/tau_c)) + sigma*sqrt(1-exp(-2*dt/tau_c))*rand_g1; 
     s_c += ds;
     //independent input cell 1
     rand_g2 = gsl_ran_gaussian (r,1.0);
     ds = (mean-s_1)*(1-exp(-dt/tau_c)) + sigma*sqrt(1-exp(-2*dt/tau_c))*rand_g2; 
     s_1 += ds;
     //independent input cell 2
     rand_g3 = gsl_ran_gaussian (r,1.0);
     ds = (mean-s_2)*(1-exp(-dt/tau_c)) + sigma*sqrt(1-exp(-2*dt/tau_c))*rand_g3; 
     s_2 += ds;

     s[0] = s_1*(1-corr_coef) + s_c*corr_coef; //input to cell 1
     s[1] = s_2*(1-corr_coef) + s_c*corr_coef; //input to cell 2
     */
     
     /*
       Model equations:
       C dv/dt = s - I_ion
       I_ion = sum( gate.*(V-E) )
       d gate/dt = -1/tau(v).*(gate-gate_0(v))
     */

    
     for (nn=0;nn<2;nn++)
       {
	 //next line: variable reversal potential for Ca
	 E_Ca[nn] = 500.0*R_F*(T + 273.15)*log(3000.0/Ca[nn]);
	 E[1][nn] = E_Ca[nn];
	 E[2][nn] = E_Ca[nn];
 	 
	 Ca_inf[nn] = 0.05 - 0.94*I_Ca[nn];  
	 // integrate Ca dynamics 
	 Ca[nn] = Ca[nn] + (1-exp(-dt/tau_Ca))*(Ca_inf[nn] - Ca[nn]);
	 
	 tau[0][0][nn] = tauX(V[nn],1.32,1.26,120.0,-25); //m, Na
	 tau[1][0][nn] = tauhNa(V[nn]); //h, Na
	 tau[0][1][nn] = tauX(V[nn],21.7,21.3,68.1,-20.5); //m, CaT
	 tau[1][1][nn] = tauX(V[nn],105.0,89.8,55.0,-16.9); //h, CaT
	 tau[0][2][nn] = taumCaS(V[nn]); //m, CaS
	 tau[1][2][nn] = tauhCaS(V[nn]); //h, CaS
	 tau[0][3][nn] = tauX(V[nn],11.6,10.4,32.9,-15.2); //m, A
	 tau[1][3][nn] = tauX(V[nn],38.6,29.2,38.9,-26.5); //h, A
	 tau[0][4][nn] = tauX(V[nn],90.3,75.1,46.0,-22.7); //m, KCa
	 tau[1][4][nn] = 0.0; //h, kCa
	 tau[0][5][nn] = tauX(V[nn],7.2,6.4,28.3,-19.2); //m, Kd
	 tau[1][5][nn] = 0.0; //h, Kd	
	 tau[0][6][nn] = tauX(V[nn],272.0,-1499.0,42.2,-8.73); //m, H
	 tau[1][6][nn] = 0.0; //h, H

	 gate_inf[0][0][nn] = Xinf(V[nn],25.5,-5.29); //m, Na
	 gate_inf[1][0][nn] = Xinf(V[nn],48.9,5.18); //h, Na
	 gate_inf[0][1][nn] = Xinf(V[nn],27.1,-7.2); //m, CaT
	 gate_inf[1][1][nn] = Xinf(V[nn],32.1,5.5); //h, CaT
	 gate_inf[0][2][nn] = Xinf(V[nn],33.0,-8.1); //m, CaS
	 gate_inf[1][2][nn] = Xinf(V[nn],60.0,6.2); //h, CaT
	 gate_inf[0][3][nn] = Xinf(V[nn],27.2,-8.7); //m, A
	 gate_inf[1][3][nn] = Xinf(V[nn],56.9,4.9); //h, A
	 gate_inf[0][4][nn] = m_KCainf(V[nn],Ca[nn]); //m, kCa 
	 gate_inf[1][4][nn] = 0.0; //h, kCa
	 gate_inf[0][5][nn] = Xinf(V[nn],12.3,-11.8); //m, Kd
	 gate_inf[1][5][nn] = 0.0; //h, Kd	
	 gate_inf[0][6][nn] = Xinf(V[nn],70.0,6.0); //m, H
	 gate_inf[1][6][nn] = 0.0; //h, H
       }
     //graded synapses
      
     S_inf[1] = Xinf(V[0],V_half,V_slope);
     S_inf[0] = Xinf(V[1],V_half,V_slope);
      
      
     //old synapses, non-differentiable
     /*
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
     */
     for (nn=0;nn<2;nn++)
       {       
	 //evolve Svar
	 Svar[nn] = Svar[nn] + (1-exp(-dt/tau_syn))*(S_inf[nn] - Svar[nn]);
	 //old from Grashow papers: "weird" function
	 //Svar[nn] = Svar[nn] + (1-exp(-dt/(tau_syn*(1-S_inf[nn]))))*(S_inf[nn] - Svar[nn]);

	 //evolve the conductances, first order Euler
	 for (ii=0;ii<7;ii++)
	   {
	     gate[0][ii][nn] = gate[0][ii][nn] + (1-exp(-dt/tau[0][ii][nn]))*(gate_inf[0][ii][nn] - gate[0][ii][nn]);
	     gate[1][ii][nn] = gate[1][ii][nn] + (1-exp(-dt/tau[1][ii][nn]))*(gate_inf[1][ii][nn] - gate[1][ii][nn]);
	   }
	 //reset
	 gate[1][4][nn] = 1.0; //h, KCa
	 gate[1][5][nn] = 1.0; //h, Kd
	 gate[1][6][nn] = 1.0; //h, H
     
	 sum_g[nn] = g_L[nn] + g_syn[nn]*Svar[nn];
	 sum_Eg[nn] = g_L[nn]*E_L + g_syn[nn]*Svar[nn]*E_syn;
	 for (ii=0;ii<7;ii++)
	   {
	     gnmh[ii][nn] = g[ii][nn]*get_power(gate[0][ii][nn],p[ii])*gate[1][ii][nn];
	     sum_g[nn] += gnmh[ii][nn];
	     sum_Eg[nn] += gnmh[ii][nn]*E[ii][nn];
	   }
	 I_Ca[nn] = (gnmh[1][nn] + gnmh[2][nn])*(V[nn]-E_Ca[nn]);

	 tau_V[nn]  = C/sum_g[nn];  //Membrane time constant.
	 
	 // V_inf is steady-state voltage.
	 V_inf[nn] = (sum_Eg[nn] + s[nn])/sum_g[nn];
	 
	 //evolve the voltage using exponential Euler
	 V[nn] = V_inf[nn] + (V[nn]-V_inf[nn])*exp(-dt/tau_V[nn]);

	     
	 //only record a small subset of the voltages
	 if (current_time > (total_time-5000)  && fmod(current_time,0.1)<dt)
	   {
	     fprintf(Vfile,"%d %f \n", nn,V[nn]);
	   }
       } //end of nn loop
   } //end of time loop
     
  fclose(Vfile);
 
 return 0;
}
