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
#define rang 0.2   //?
#define dt 0.1            //time step in ms 
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

  /*For random number generation:*/
  gsl_rng * r;
  r = gsl_rng_alloc (gsl_rng_default);
  gsl_rng_set(r, time(NULL));
  /*done*/

  //properties of the common input, for no it's not there
  double mean = atof(argv[1]);
  double sigma = atof(argv[2]);
  double tau_c = atof(argv[3]);
  double corr_coef = atof(argv[4]);
  int Ntot = atof(argv[5]);

  double taug = 100;
  double taum = 10*taug;
  double tstop = 1200*taug; //changed: JG Nov 11, 2015 from 600*taug
  double ms_per_samp = tstop/2000.0;

  //example: bursting cell
  double tau_i[2][7];
  double tau_Na = 15.0;  //Na
  double tau_CaT = 0.18;  //CaT
  double tau_CaS = 0.22;  //CaS
  double tau_A = 2.0;   //A
  double tau_KCa = 8.0;   //KCa
  double tau_Kd = 5.0;   //Kd
  double tau_h = 0.08;  //h
  double tau_syni[2];
  //tau_syni[0] = 0.001; tau_syni[1] = 0.001;
  tau_syni[0] = 0.01; tau_syni[1] = 0.01;
  //double target = 10.0;
  double target = 30.0;

  //DECREASING TARGET DECREASES gsyn
  //DECREASING tau_isyn DECREASES gsyn
  
  /*---------------------------------------------------------------------------------------------------------------
    so it turns out that target 30 and large tau_isyn 0.01 generates 3 sp alternating bursts with DC=5 (g_syn around 0.5-0.6)
    DC5_tau0.01_target30_RUN1-2.fig
    
    Decreasing target to 7 (with tau_isyn tp 0.01) generates sometimes oscillating weights (smaller g_syn 0.09-0.11)
    and alternating 2 sp bursts spikes in V1 and V2
    DC5_tau0.01_target7_RUN1-3.fig

    Decrease target further for 1 sp?
    if target=5, still 2 sp, but not alternating sometimes
    The further we decrease target, the less Na there is, and so the spikes are very short.
    
    Decreasing tau_isyn to 0.001 but target at 30 generates stable weights (smaller g_syn 0.07-0.08)
    and alternating single spikes in V1 and V2
    DC5_tau0.001_target30_RUN1-3.fig
    
    Reducing the target to 10 (with tau_isyn to 0.001), generates stable weights, 
    but cells are either 1-2 sp with irregular bursts (smaller g_syn, 0.013-0.020)
    DC5_tau0.001_target10_RUN1-3.fig
    
    Maybe keep on reducing the target?
    Reducing the target to 7 (with tau_isyn to 0.001), generates stable weights, 
    but cells are 2 sp with irregular bursts (smaller g_syn, 0.007-0.009)
    DC5_tau0.001_target7_RUN1-3.fig
    ---------------------------------------------------------------------------------------------------------------

    DC=0, and tau_isyn=0.001 with target 7 produces 3 sp burts that don't really alternate (small g_syn 0.009-0.017)
    - but not always, sometimes one of the cells stops if g_syn is very different between the two cells 
    DC0_tau0.001_target7_RUN1-6.fig

    Decrease to 5 more unreliable. 
    DC0_tau0.001_target5_RUN1-3.fig

    But increase to 10 more reliable? But not always. Once silences cell. 
    DC0_tau0.001_target10_RUN1-4.fig
    CHECKED THAT WHEN gsyn=0, ALSO BURSTING

    But increase to 12 more reliable? But not always alternating. CHECK MORE
    DC0_tau0.001_target12_RUN1-4.fig
    CHECKED THAT WHEN gsyn=0, ALSO BURSTING

    Increase to 15 is too much, start to get 2 sp bursts.
    DC0_tau0.001_target15_RUN1-3.fig

    Increase target to 30 (tau_isyn=0.001), 2 sp alternating bursts (g_syn 0.05, doesn't vary much)
    Unbounded growth of the weights, even if sometimes bounded, one set of conductances grow huge!

    but when large tau_isyn (0.01) produces large oscillations in the conductances for any target
    ---------------------------------------------------------------------------------------------------------------  */

  //initial conditions for the conductances
  /*
  double g_Na = rang*gsl_ran_flat (r,0.0,1.0); //Na 0.1033;
  double g_CaT = rang*gsl_ran_flat (r,0.0,1.0); //CaT  0.1359;
  double g_CaS = rang*gsl_ran_flat (r,0.0,1.0); //CaS 0.1082;
  double g_A = rang*gsl_ran_flat (r,0.0,1.0); //A 0.0307;
  double g_KCa = rang*gsl_ran_flat (r,0.0,1.0); //KCa 0.1907;
  double g_Kd = rang*gsl_ran_flat (r,0.0,1.0); //Kd 0.1405;
  double g_h = rang*gsl_ran_flat (r,0.0,1.0); //h 0.0073;
  double g_leak = 0.001 + 0.199*gsl_ran_flat (r,0.0,1.0); //leak 0.1790;
  */
  double gbar_vec[2][7]; double g_syn[2];
  double m[2][7]; double msyn[2];

  for (int num_run=1;num_run<Ntot+1;num_run++)
    {
      //for now, make synaptic conductance be 0 and treat the two cells as independent
      //double gsyn = 0.0;
 
      /*open files where conductances will be recorded*/
      char sp_file1 [999];
      sprintf (sp_file1, "g1_mean_%g_sig_%g_tau_%g_coef_%g_target_%g_taug_%g_taum_%g_num_%d.dat",mean,sigma,tau_c,corr_coef,target,taug,taum,num_run);
      FILE *spfile1;
      spfile1 = fopen(sp_file1,"w");
      //
      char sp_file2 [999];
      sprintf (sp_file2, "g2_mean_%g_sig_%g_tau_%g_coef_%g_target_%g_taug_%g_taum_%g_num_%d.dat",mean,sigma,tau_c,corr_coef,target,taug,taum,num_run);
      FILE *spfile2;
      spfile2 = fopen(sp_file2,"w");

      int ii, jj, nn;
  
      //reversal potentials
      double E_Na = 50; //mV or 30 - to reproduce some of Tim's figure 3, need to make this 30 mV
      double E_A = -80; //mV
      double E_K = -80; //mV
      double E_H = -20; //mV
      double E_Ca[2]; E_Ca[0] = 120; E_Ca[1] = 120; //mV default, but will vary
      double E_L = -50; //mV
      double C = 1; //nF 
      double tau_Ca = 20; //ms
      double I_Ca[2]; I_Ca[0] = 0; I_Ca[1] = 0;
      //synaptic
      double E_syn = -80; //mV
      double tau_syn = 100; //ms
      //Grashow V_half=-45; V_slope = 5;
      //otherwise I used:  V_half = 40; V_slope = -2.0;
      double V_half = -45.0;  //mV
      double V_slope = 5.0;  //mV

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

	  tau_i[nn][0] = tau_Na;  gbar_vec[nn][0] = rang*gsl_ran_flat (r,0.0,1.0); 
	  tau_i[nn][1] = tau_CaT; gbar_vec[nn][1] = rang*gsl_ran_flat (r,0.0,1.0);
	  tau_i[nn][2] = tau_CaS; gbar_vec[nn][2] = rang*gsl_ran_flat (r,0.0,1.0);
	  tau_i[nn][3] = tau_A;   gbar_vec[nn][3] = rang*gsl_ran_flat (r,0.0,1.0);   
	  tau_i[nn][4] = tau_KCa; gbar_vec[nn][4] = rang*gsl_ran_flat (r,0.0,1.0); 
	  tau_i[nn][5] = tau_Kd;  gbar_vec[nn][5] = rang*gsl_ran_flat (r,0.0,1.0); 
	  tau_i[nn][6] = tau_h;   gbar_vec[nn][6] = rang*gsl_ran_flat (r,0.0,1.0);  

	  E[0][nn] = E_Na;      m[nn][0] = gbar_vec[nn][0]; //Na
	  E[1][nn] = E_Ca[nn];  m[nn][1] = gbar_vec[nn][1];  //CaT
	  E[2][nn] = E_Ca[nn];  m[nn][2] = gbar_vec[nn][2]; //CaS
	  E[3][nn] = E_A;       m[nn][3] = gbar_vec[nn][3]; //A
	  E[4][nn] = E_K;       m[nn][4] = gbar_vec[nn][4]; //KCa
	  E[5][nn] = E_K;       m[nn][5] = gbar_vec[nn][5]; //Kd
	  E[6][nn] = E_H;       m[nn][6] = gbar_vec[nn][6]; //h

	  g_syn[nn] = 0.001*gsl_ran_flat (r,0.0,1.0); 
	  msyn[nn] = g_syn[nn];
	}
 
      //
      p[0] = p_Na;
      p[1] = p_CaT;
      p[2] = p_CaS;
      p[3] = p_A;
      p[4] = p_KCa;
      p[5] = p_K;
      p[6] = p_H; 


      double g_L[2]; g_L[0] = 0.001 + 0.199*gsl_ran_flat (r,0.0,1.0); g_L[1] = 0.001 + 0.199*gsl_ran_flat (r,0.0,1.0);
      double current_time = 0.0;
      double s_c = mean;
      double s_1 = mean;
      double s_2 = mean;
      double s[2]; 
      double Svar[2];
      Svar[1] = Xinf(V[0],45.0,-2.0);
      Svar[0] = Xinf(V[1],45.0,-2.0); 
      double S_inf[2];
      S_inf[0] = Svar[0]; 
      S_inf[1] = Svar[1];
      int counter = 0;

      int voltage_high1 = 0;
      double spike1 = 0.0; //to impose an absolute refractory period when reading spikes
      int voltage_high2 = 0;
      double spike2 = 0.0; //to impose an absolute refractory period when reading spikes
      double rand_g1 = 0.0;
      double rand_g2 = 0.0;
      double rand_g3 = 0.0;

      double ds = 0.0;
 
      double sum_g[2], sum_Eg[2], tau_V[2], V_inf[2];
      for (int gg=0;gg<7;gg++)
	{
	  printf("%f \t", gbar_vec[0][gg]);
	  fprintf(spfile1,"%f \t", gbar_vec[0][gg]);
	  fprintf(spfile2,"%f \t", gbar_vec[1][gg]);
	}
      fprintf(spfile1,"%f \n",g_syn[0]);
      fprintf(spfile2,"%f \n",g_syn[1]);
      printf("\n");
  
      while (current_time < tstop)
	{       
	  current_time += dt;
	  counter ++;
	  //printf("%f \n",current_time);
	  //s = mean + sigma/sqrt(dt)*rand_g; //if white noise

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
	  /*
	    S_inf[1] = Xinf(V[0],V_half,V_slope);
	    S_inf[0] = Xinf(V[1],V_half,V_slope);
	  */
      
	  //old synapses, non-differentiable
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
	 
	      //evolve Svar: new
	      //Svar[nn] = Svar[nn] + (1-exp(-dt/tau_syn))*(S_inf[nn] - Svar[nn]);
	      //old from Grashow papers: "weird" function
	      Svar[nn] = Svar[nn] + (1-exp(-dt/(tau_syn*(1-S_inf[nn]))))*(S_inf[nn] - Svar[nn]);

	      //evolve the conductances, first order Euler
	      for (ii=0;ii<7;ii++)
		{
		  gate[0][ii][nn] += (1-exp(-dt/tau[0][ii][nn]))*(gate_inf[0][ii][nn] - gate[0][ii][nn]);
		  gate[1][ii][nn] += (1-exp(-dt/tau[1][ii][nn]))*(gate_inf[1][ii][nn] - gate[1][ii][nn]);
		}
	      //reset because these conductances don't have an h gate
	      gate[1][4][nn] = 1.0; //h, KCa
	      gate[1][5][nn] = 1.0; //h, Kd
	      gate[1][6][nn] = 1.0; //h, H
     
	      sum_g[nn] = g_L[nn] + g_syn[nn]*Svar[nn];
	      sum_Eg[nn] = g_L[nn]*E_L + g_syn[nn]*Svar[nn]*E_syn;
	      for (ii=0;ii<7;ii++)
		{
		  //reset the maximal conductances:
		  g[ii][nn] = gbar_vec[nn][ii];
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

	      //now, evolve/grow the synaptic conductances
	      for (int gg=0;gg<7;gg++)
		{
		  gbar_vec[nn][gg] += (1-exp(-dt/taug))*(m[nn][gg] -  gbar_vec[nn][gg]);
		  //gbar_vec[nn][gg] += dt/taug*(m[nn][gg] - gbar_vec[nn][gg]);
		  if (gbar_vec[nn][gg] < 0)
		    {
		      gbar_vec[nn][gg] = 0;
		      //and don't change m
		    }
		  else //change m
		    {
		      m[nn][gg] += (1-exp(-dt/(-taum/tau_i[nn][gg])))*(Ca[nn] -  target);
		      //m[nn][gg] += dt*(-tau_i[nn][gg]/taum)*(Ca[nn] - target);
		    }
		}
	 
	      //also evolve the synaptic conductances:
	      g_syn[nn] += (1-exp(-dt/taug))*(msyn[nn] -  g_syn[nn]);
	      msyn[nn] += (1-exp(-dt/(-taum/tau_syni[nn])))*(Ca[nn] -  target);
	    }

     
	  if (fmod(current_time,ms_per_samp)<dt)
	    {
	      for (int gg=0;gg<7;gg++)
		{	    
		  fprintf(spfile1,"%f \t", gbar_vec[0][gg]);
		  fprintf(spfile2,"%f \t", gbar_vec[1][gg]);
		}
	      fprintf(spfile1,"%f \t %f \n",g_L[nn],g_syn[0]);
	      fprintf(spfile2,"%f \t %f \n",g_L[nn],g_syn[1]);
	    }
	}
     
      fclose(spfile1);
      fclose(spfile2);

      for (nn=0;nn<2;nn++)
	{
	  for (ii=0;ii<7;ii++)
	    {
	      printf("%f \t",gbar_vec[nn][ii]);
	    }
	  printf("%f \t %f \n",g_L[nn],g_syn[nn]);
	}
    }
							       
 
 return 0;
}
