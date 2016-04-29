/*Code to run the STG model neuron simulation.
NO CA DYNAMICS - DOES NOT CONVERGE
Adapt the maximal intrinsic conductances of some currents of each cell in a half-center 
oscillator made up of STG neurons to counteract the change in the synaptic conductance.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
//#include </n/sw/gsl-1.15_gcc-4.7.2/include/gsl/gsl_randist.h>
//#include </n/sw/gsl-1.15_gcc-4.7.2/include/gsl/gsl_rng.h>
#include <time.h>

/*Define some parameters */
#define total_time 2000.0   //the time of the presynaptic input in ms
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

double tau_hNa(double V0)
{
  double X = (0.67/(1.0+exp((V0+62.9)/-10.0)))*(1.5 + 1.0/(1.0+exp((V0+34.9)/3.6)));
  return X;
}

double taum_CaS(double V0)
{
  double X = 1.4 + (7.0/((exp((V0+27.0)/10.0))+(exp((V0+70.0)/-13.0))));
  return X;
}

double tauh_CaS(double V0)
{
  double X = 60.0 + (150.0/((exp((V0+55)/9.0))+(exp((V0+65.0)/-16.0))));
  return X;
}

double dXinf(double V0, double A, double B)
{
  double X = -exp((A + V0)/B)/(B*(exp((A + V0)/B) + 1)*(exp((A + V0)/B) + 1));
  return X;
}

double dmKCainf(double V0, double Ca0)
{
  double X = (5.0*Ca0)/(63.0*exp((5.0*V0)/63.0 + 283.0/126.0)*(Ca0 + 3.0)*(1.0/exp((5.0*V0)/63.0 + 283.0/126.0) + 1)*(1.0/exp((5.0*V0)/63.0 + 283.0/126.0) + 1));
  return X;
}

void compare_tau(double a[2], double b[2], double tau[2], double tauf[2], double taus[2], double tauu[2])
{
  for (int ii=0;ii<2; ii++)
    {
      if (tau[ii]<=tauf[ii])
	{
	  a[ii] = 1;
	  b[ii] = 1;
	}
      if (tauf[ii]<tau[ii] && tau[ii]<=taus[ii])
	{
	  a[ii] = (log(taus[ii])-log(tau[ii]))/(log(taus[ii])-log(tauf[ii]));
	  b[ii] = 1;
	}
      if (taus[ii]<tau[ii] && tau[ii]<=tauu[ii])
	{
	  a[ii] = 0;
	  b[ii] = (log(tauu[ii])-log(tau[ii]))/(log(tauu[ii])-log(taus[ii]));
	}
      if (tauu[ii]<tau[ii])
	{
	  a[ii] = 0;
	  b[ii] = 0;
	}
    }
      return;
}

int find_coef(double Vv[300001], double A)
{
  int X = 0;
  while (Vv[X]!=A && X<300001)
    {
      X++;
    }
  return X;
}

void get_inverse(double A[4][4],double Ainv[4][4])
{

  double a11 = A[0][0]; double a12 = A[0][1]; double a13 = A[0][2]; double a14 = A[0][3];
  double a21 = A[1][0]; double a22 = A[1][1]; double a23 = A[1][2]; double a24 = A[1][3];
  double a31 = A[2][0]; double a32 = A[2][1]; double a33 = A[2][2]; double a34 = A[2][3];
  double a41 = A[3][0]; double a42 = A[3][1]; double a43 = A[3][2]; double a44 = A[3][3];
  
  Ainv[0][0] = (a22*a33*a44 - a22*a34*a43 - a23*a32*a44 + a23*a34*a42 + a24*a32*a43 - a24*a33*a42)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);
  
  Ainv[0][1] = -(a12*a33*a44 - a12*a34*a43 - a13*a32*a44 + a13*a34*a42 + a14*a32*a43 - a14*a33*a42)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);
  
  Ainv[0][2] = (a12*a23*a44 - a12*a24*a43 - a13*a22*a44 + a13*a24*a42 + a14*a22*a43 - a14*a23*a42)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);
  
  Ainv[0][3] = -(a12*a23*a34 - a12*a24*a33 - a13*a22*a34 + a13*a24*a32 + a14*a22*a33 - a14*a23*a32)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[1][0] = -(a21*a33*a44 - a21*a34*a43 - a23*a31*a44 + a23*a34*a41 + a24*a31*a43 - a24*a33*a41)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[1][1] = (a11*a33*a44 - a11*a34*a43 - a13*a31*a44 + a13*a34*a41 + a14*a31*a43 - a14*a33*a41)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[1][2] = -(a11*a23*a44 - a11*a24*a43 - a13*a21*a44 + a13*a24*a41 + a14*a21*a43 - a14*a23*a41)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[1][3] = (a11*a23*a34 - a11*a24*a33 - a13*a21*a34 + a13*a24*a31 + a14*a21*a33 - a14*a23*a31)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[2][0] = (a21*a32*a44 - a21*a34*a42 - a22*a31*a44 + a22*a34*a41 + a24*a31*a42 - a24*a32*a41)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[2][1] = -(a11*a32*a44 - a11*a34*a42 - a12*a31*a44 + a12*a34*a41 + a14*a31*a42 - a14*a32*a41)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[2][2] = (a11*a22*a44 - a11*a24*a42 - a12*a21*a44 + a12*a24*a41 + a14*a21*a42 - a14*a22*a41)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[2][3] = -(a11*a22*a34 - a11*a24*a32 - a12*a21*a34 + a12*a24*a31 + a14*a21*a32 - a14*a22*a31)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[3][0] = -(a21*a32*a43 - a21*a33*a42 - a22*a31*a43 + a22*a33*a41 + a23*a31*a42 - a23*a32*a41)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[3][1] = (a11*a32*a43 - a11*a33*a42 - a12*a31*a43 + a12*a33*a41 + a13*a31*a42 - a13*a32*a41)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[3][2] = -(a11*a22*a43 - a11*a23*a42 - a12*a21*a43 + a12*a23*a41 + a13*a21*a42 - a13*a22*a41)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  Ainv[3][3] = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31)/(a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41);

  return;
}
      

int main(int argc, char *argv[])
{
  double mean = atof(argv[1]);
  double sigma = atof(argv[2]);
  //intrinsic
  double gNa = atof(argv[3]); //uS/nF 
  double gCaT = atof(argv[4]);
  double gCaS = atof(argv[5]);
  double gA = atof(argv[6]);
  double gKCa = atof(argv[7]);
  double gKd = atof(argv[8]);
  double gH = atof(argv[9]);
  double gleak = atof(argv[10]); //uS/nF    
  //synaptic
  double gsyn =  atof(argv[11]);
  double new_gsyn = atof(argv[12]);
  double corr_coef = 0;//atof(argv[12]); //correlation coefficient of the input
  double tau_c = 0;//atof(argv[13]); //correlation time const
  int num_run = 0;//atof(argv[14]);
  
  double Iapp = mean;
  //reversal potentials
  double VNa = 50; //mV or 30?
  double VA = -80; //mV
  double VK = -80; //mV
  double VH = -20; //mV
  double VCa = 120; //mV 
  double Vleak = -50; //mV
  double C = 1; //nF 
  double tauCa = 20; //ms
  //synaptic
  double Vsyn = -80; //mV
  double tausyn = 100; //ms
  double V_half = 45;
  double V_slope = -2;
  //Initiate V
  //Compute the different time constants - all vectors, at each value of V
  double V0[300001];

  double Ca = 10.0;
  double Istatic0[300001];
  
  int nn = 0; //for determining Vosc
  double Vosc = 0;
  double Vzero[3]; Vzero[0]=0; Vzero[1]=0; Vzero[2]=0;
  for (int ii=0;ii<300001;ii++)
    {
      V0[ii] = -100+ii*0.001;
       //Compute the static current (I/V curve)
      Istatic0[ii] = -gNa*get_power(Xinf(V0[ii],25.5,-5.29),3)*Xinf(V0[ii],48.9,5.18)*(V0[ii]-VNa) - gCaT*get_power(Xinf(V0[ii],27.1,-7.2),3)*Xinf(V0[ii],32.1,5.5)*(V0[ii]-VCa) - gCaS*get_power(Xinf(V0[ii],33,-8.1),3)*Xinf(V0[ii],60,6.2)*(V0[ii]-VCa) - gA*get_power(Xinf(V0[ii],27.2,-8.7),3)*Xinf(V0[ii],56.9,4.9)*(V0[ii]-VK) -gKCa*get_power(m_KCainf(V0[ii],Ca),4)*(V0[ii]-VK) - gKd*get_power(Xinf(V0[ii],12.3,-11.8),4)*(V0[ii]-VK) - gH*Xinf(V0[ii],70.0,6.0)*(V0[ii]-VH) - gleak*(V0[ii]-Vleak) - gsyn*Xinf(V0[ii],40,-2)*(V0[ii]-Vsyn) + Iapp;

      /* Look for the zeros of the I/V curve. If there are 3 zeros, V_osc is the
	 more depolarized (the third). If there is only one zero, it corresponds
	 to V_osc */
      if (ii>0)
	{
	  if ((Istatic0[ii-1]>0 && Istatic0[ii]<=0) || (Istatic0[ii-1]<0 && Istatic0[ii]>=0))
	    {
	      Vzero[nn] = V0[ii-1];
	      nn++;
	    }
	}
    }
      
  if (Vzero[1]!=0)
    {
      Vosc = Vzero[2];
    }
  else
    {
      Vosc = Vzero[0];
    }
  double VTC = -50; //fzero(@(V)Jx(V),-50);
  /*************************************************************************************************/
  //Initialize V
  double V[2]; V[0] = VTC; V[1] = Vosc;

  double Istatic[2];
  int X11 = find_coef(V0,VTC);
  int X12 = find_coef(V0,Vosc); 
  Istatic[0] = Istatic0[X11]; Istatic[1] = Istatic0[X12];
  //printf("%f %f %f %f \n",VTC,V0[X11],Vosc,V0[X12]);  

  //Compute the different time constants - all vectors, at each value of V
  double taumNa[2], tauhNa[2], taumCaT[2], tauhCaT[2], taumCaS[2], tauhCaS[2], taumA[2], tauhA[2], taumKCa[2], taumKd[2], taumH[2], taumsyn[2];
  double mNainf[2], hNainf[2], mCaTinf[2], hCaTinf[2], mCaSinf[2], hCaSinf[2], mAinf[2], hAinf[2], mKCainf[2], mKdinf[2], mHinf[2], msyninf[2];
  double BAL_mNa[2], BAL_hNa[2], BAL_mCaT[2], BAL_hCaT[2], BAL_mCaS[2], BAL_hCaS[2], BAL_mA[2], BAL_hA[2], BAL_mKCa[2], BAL_mKd[2], BAL_mH[2], BAL_msyn[2], newBAL_msyn[2];
  
  //Choose the fast, slow and ultraslow timescales
  //taumNa -Fastest (regenerative)
  //taumKd -Fastest restorative
  //tauhCaS -Slowest
  
  // Initialize a (w_fs) and b (w_su)
  double amNa[2], ahNa[2], amCaT[2], ahCaT[2], amCaS[2], ahCaS[2], amA[2], ahA[2], amKCa[2], amKd[2], amH[2], amsyn[2];
  double bmNa[2], bhNa[2], bmCaT[2], bhCaT[2], bmCaS[2], bhCaS[2], bmA[2], bhA[2], bmKCa[2], bmKd[2], bmH[2], bmsyn[2];
  
  double ICa[2];
  
  for (int ii=0;ii<2;ii++)
    {      
      /* 1. Aggregate timescales
	 This part computes the contributions of the different variables in the three timescales (a = w_fs and b = w_su). 
	 fast: a=1 and b=1, slow: a=0 and b=1, ultraslow: a=0 and b=0.
	 Between fast and slow: 0<a<1 and b=1, between slow and ultraslow: a=0 and 0<b<1. */
      taumNa[ii] = tauX(V[ii],1.32,1.26,120,-25);
      tauhNa[ii] = tau_hNa(V[ii]);
      taumCaT[ii] = tauX(V[ii],21.7,21.3,68.1,-20.5);
      tauhCaT[ii] = tauX(V[ii],105,89.8,55,-16.9);
      taumCaS[ii] = taum_CaS(V[ii]);
      tauhCaS[ii] = tauh_CaS(V[ii]);
      taumA[ii] = tauX(V[ii],11.6,10.4,32.9,-15.2);
      tauhA[ii] = tauX(V[ii],38.6,29.2,38.9,-26.5);
      taumKCa[ii] = tauX(V[ii],90.3,75.1,46,-22.7);
      taumKd[ii] = tauX(V[ii],7.2,6.4,28.3,-19.2);
      taumH[ii] = tauX(V[ii],272,-1499,42.2,-8.73);
      taumsyn[ii] = tausyn;

      // 2. Compute V_TC and Ca_TC
      mNainf[ii] = Xinf(V[ii],25.5,-5.29);
      hNainf[ii] = Xinf(V[ii],48.9,5.18);
      mCaTinf[ii] = Xinf(V[ii],27.1,-7.2);
      hCaTinf[ii] = Xinf(V[ii],32.1,5.5);
      mCaSinf[ii] = Xinf(V[ii],33,-8.1);
      hCaSinf[ii] = Xinf(V[ii],60,6.2);
      mAinf[ii] = Xinf(V[ii],27.2,-8.7);
      hAinf[ii] = Xinf(V[ii],56.9,4.9);
      mKdinf[ii] = Xinf(V[ii],12.3,-11.8);
      mHinf[ii] = Xinf(V[ii],70.0,6.0);
      ICa[ii] = (gCaT*get_power(mCaTinf[ii],3)*hCaTinf[ii] + gCaS*get_power(mCaSinf[ii],3)*hCaSinf[ii])*(V[ii]-VCa);
      mKCainf[ii] = m_KCainf(V[ii],Ca); 
      msyninf[ii] = Xinf(V[ii],40,-2); //treat V1=V2 for the two neurons

      //Compute the static current (I/V curve)
      //no need for this HERE
      //because it's computed above and used to find Vosc
      //Istatic[ii] = -gNa*get_power(mNainf[ii],3)*hNainf[ii]*(V[ii]-VNa) - gCaT*get_power(mCaTinf[ii],3)*hCaTinf[ii]*(V[ii]-VCa) - gCaS*get_power(mCaSinf[ii],3)*hCaSinf[ii]*(V[ii]-VCa) - gA*get_power(mAinf[ii],3)*hAinf[ii]*(V[ii]-VK) -gKCa*get_power(mKCainf[ii],4)*(V[ii]-VK) - gKd*get_power(mKdinf[ii],4)*(V[ii]-VK) - gH*mHinf[ii]*(V[ii]-VH) - gleak*(V[ii]-Vleak) -gsyn*msyninf[ii]*(V[ii]-Vsyn) + Iapp;

      // 3. Compute each balance at each V
      //Computation of dVdot/dX*dX,inf/dV for each variable X

      BAL_mNa[ii] = gNa * (-3*get_power(mNainf[ii],2)*hNainf[ii]*(V[ii]-VNa) * dXinf(V[ii],25.5,-5.29));
      BAL_hNa[ii] = gNa * (-get_power(mNainf[ii],3)*(V[ii]-VNa) * dXinf(V[ii],48.9,5.18));
      BAL_mCaT[ii] = gCaT * (-3*get_power(mCaTinf[ii],2)*hCaTinf[ii]*(V[ii]-VCa) * dXinf(V[ii],27.1,-7.2));
      BAL_hCaT[ii] = gCaT * (-get_power(mCaTinf[ii],3)*(V[ii]-VCa) * dXinf(V[ii],32.1,5.5));
      BAL_mCaS[ii] = gCaS * (-3*get_power(mCaSinf[ii],2)*hCaSinf[ii]*(V[ii]-VCa) * dXinf(V[ii],33,-8.1));
      BAL_hCaS[ii] = gCaS * (-get_power(mCaSinf[ii],3)*(V[ii]-VCa) * dXinf(V[ii],60,6.2));
      BAL_mA[ii] = gA * (-3*get_power(mAinf[ii],2)*hAinf[ii]*(V[ii]-VK) * dXinf(V[ii],27.2,-8.7));
      BAL_hA[ii] = gA * (-get_power(mAinf[ii],3)*(V[ii]-VK) * dXinf(V[ii],56.9,4.9));
      BAL_mKd[ii] = gKd * (-4*get_power(mKdinf[ii],3)*(V[ii]-VK) * dXinf(V[ii],12.3,-11.8));
      BAL_mH[ii] = gH * (-dXinf(V[ii],70.0,6.0) * (V[ii]-VH));
      BAL_mKCa[ii] = gKCa*(-4*get_power(mKCainf[ii],3)*(V[ii]-VK) * dmKCainf(V[ii],Ca));
      //may need to make the next a minus sign!
      BAL_msyn[ii] = gsyn*(-dXinf(V[ii],V_half,V_slope)*(V[ii]-Vsyn)); //treat V1=V2 for the two neurons
      newBAL_msyn[ii] = new_gsyn*(-dXinf(V[ii],V_half,V_slope)*(V[ii]-Vsyn)); //treat V1=V2 for the two neurons
    }

  /*Algorithm that computes the values of a and b
    Choose the fast, slow and ultraslow timescales
    taumNa - Fastest regenerative
    taumH - Fastest restorative 
    tauhCaS - Slowest  (before it was tauhCaS)*/
  double tau_fast[2]; 
  double tau_slow[2];
  double tau_ultra[2];
  for (int ii=0;ii<2;ii++)
    {
      tau_fast[ii] = taumNa[ii];
      tau_slow[ii] = taumKd[ii];
      tau_ultra[ii] = taumH[ii];
      //printf("%d %f %f %f %f \n",ii,tau_fast[ii],tau_slow[ii],tau_ultra[ii],tauhCaS[ii]);
    }
  compare_tau(amNa,bmNa,taumNa,tau_fast,tau_slow,tau_ultra);
  compare_tau(ahNa,bhNa,tauhNa,tau_fast,tau_slow,tau_ultra);
  compare_tau(amCaT,bmCaT,taumCaT,tau_fast,tau_slow,tau_ultra);
  compare_tau(ahCaT,bhCaT,tauhCaT,tau_fast,tau_slow,tau_ultra);
  compare_tau(amCaS,bmCaS,taumCaS,tau_fast,tau_slow,tau_ultra);
  compare_tau(ahCaS,bhCaS,tauhCaS,tau_fast,tau_slow,tau_ultra);
  compare_tau(amA,bmA,taumA,tau_fast,tau_slow,tau_ultra);
  compare_tau(ahA,bhA,tauhA,tau_fast,tau_slow,tau_ultra);
  compare_tau(amKd,bmKd,taumKd,tau_fast,tau_slow,tau_ultra);
  compare_tau(amH,bmH,taumH,tau_fast,tau_slow,tau_ultra);
  compare_tau(amKCa,bmKCa,taumKCa,tau_fast,tau_slow,tau_ultra);
  compare_tau(amsyn,bmsyn,taumsyn,tau_fast,tau_slow,tau_ultra);

  double BALf[2], BALs[2], BALu[2];
  
  for (int ii=0;ii<2;ii++)
    {
      BALf[ii] = amNa[ii]*BAL_mNa[ii] + ahNa[ii]*BAL_hNa[ii] + amCaT[ii]*BAL_mCaT[ii] + ahCaT[ii]*BAL_hCaT[ii] + amCaS[ii]*BAL_mCaS[ii] + ahCaS[ii]*BAL_hCaS[ii] + amA[ii]*BAL_mA[ii] + ahA[ii]*BAL_hA[ii] + amKd[ii]*BAL_mKd[ii] + amH[ii]*BAL_mH[ii] + amKCa[ii]*BAL_mKCa[ii] + amsyn[ii]*BAL_msyn[ii];
      
      BALs[ii] = (bmNa[ii]-amNa[ii])*BAL_mNa[ii] + (bhNa[ii]-ahNa[ii])*BAL_hNa[ii] + (bmCaT[ii]-amCaT[ii])*BAL_mCaT[ii] + (bhCaT[ii]-ahCaT[ii])*BAL_hCaT[ii] + (bmCaS[ii]-amCaS[ii])*BAL_mCaS[ii] + (bhCaS[ii]-ahCaS[ii])*BAL_hCaS[ii] + (bmA[ii]-amA[ii])*BAL_mA[ii] + (bhA[ii]-ahA[ii])*BAL_hA[ii] + (bmKd[ii]-amKd[ii])*BAL_mKd[ii] + (bmH[ii]-amH[ii])*BAL_mH[ii] + (bmKCa[ii]-amKCa[ii])*BAL_mKCa[ii] + (bmsyn[ii]-amsyn[ii])*BAL_msyn[ii];

      BALu[ii] = (1-bmNa[ii])*BAL_mNa[ii] + (1-bhNa[ii])*BAL_hNa[ii] + (1-bmCaT[ii])*BAL_mCaT[ii] + (1-bhCaT[ii])*BAL_hCaT[ii] + (1-bmCaS[ii])*BAL_mCaS[ii] + (1-bhCaS[ii])*BAL_hCaS[ii] + (1-bmA[ii])*BAL_mA[ii] + (1-bhA[ii])*BAL_hA[ii] + (1-bmKd[ii])*BAL_mKd[ii] + (1-bmH[ii])*BAL_mH[ii] + (1-bmKCa[ii])*BAL_mKCa[ii] + (1-bmsyn[ii])*BAL_msyn[ii];
    }

  //for debugging
  /*  
  for (int ii=0;ii<2;ii++)
    {      
      printf("%f %f \n",amNa[ii],bmNa[ii]);
      printf("%f %f \n",ahNa[ii],bhNa[ii]);
      printf("%f %f \n",amCaT[ii],bmCaT[ii]);
      printf("%f %f \n",ahCaT[ii],bhCaT[ii]);
      printf("%f %f \n",amCaS[ii],bmCaS[ii]);
      printf("%f %f \n",ahCaS[ii],bhCaS[ii]);
      printf("%f %f \n",amA[ii],bmA[ii]);
      printf("%f %f \n",ahA[ii],bhA[ii]);
      printf("%f %f \n",amKd[ii],bmKd[ii]);
      printf("%f %f \n",amKCa[ii],bmKCa[ii]);
      printf("%f %f \n",amsyn[ii],bmsyn[ii]);
      
      printf("%f \n",BAL_mNa[ii]);
      printf("%f \n",BAL_hNa[ii]);
      printf("%f \n",BAL_mCaT[ii]);
      printf("%f \n",BAL_hCaT[ii]);
      printf("%f \n",BAL_mCaS[ii]);
      printf("%f \n",BAL_hCaS[ii]);
      printf("%f \n",BAL_mA[ii]);
      printf("%f \n",BAL_hA[ii]);
      printf("%f \n",BAL_mKd[ii]);
      printf("%f \n",BAL_mKCa[ii]);
      printf("%f \n",BAL_msyn[ii]);
      printf("%f \n",newBAL_msyn[ii]);  
    }
  */
    
  int X1 = 0; int X2 = 1;
  
  double A[4][4];
  double b[4];
  double c[4];
  //A*b = c with b = [Iapp gH gA gKCa]'
  A[0][0] = 0; 
  A[1][0] = 0; 
  A[2][0] = 0; 
  A[3][0] = 1; 
  //
  A[0][1] = (bmH[X1]-amH[X1]) * (-dXinf(V[X1],70.0,6.0) * (V[X1]-VH));
  A[1][1] = (bmH[X2]-amH[X2]) * (-dXinf(V[X2],70.0,6.0) * (V[X2]-VH));
  A[2][1] = (1-bmH[X1]) * (-dXinf(V[X1],70.0,6.0) * (V[X1]-VH));
  A[3][1] = -mHinf[X1]*(VTC-VH);
  //
  A[0][2] = (bmA[X1]-amA[X1]) * (-3*get_power(mAinf[X1],2)*hAinf[X1]*(V[X1]-VK) * dXinf(V[X1],27.2,-8.7)) + (bhA[X1]-ahA[X1]) * (-get_power(mAinf[X1],3)*(V[X1]-VK) * dXinf(V[X1],56.9,4.9));
  A[1][2] = (bmA[X2]-amA[X2]) * (-3*get_power(mAinf[X2],2)*hAinf[X2]*(V[X2]-VK) * dXinf(V[X2],27.2,-8.7)) + (bhA[X2]-ahA[X2]) * (-get_power(mAinf[X2],3)*(V[X2]-VK) * dXinf(V[X2],56.9,4.9));
  A[2][2] = (1-bmA[X1]) * (-3*get_power(mAinf[X1],2)*hAinf[X1]*(V[X1]-VK) * dXinf(V[X1],27.2,-8.7)) + (1-bhA[X1]) *  (-get_power(mAinf[X1],3)*(V[X1]-VK) * dXinf(V[X1],56.9,4.9));
  A[3][2] = -get_power(mAinf[X1],3)*hAinf[X1]*(VTC-VK);
  //
  A[0][3] = (bmKCa[X1]-amKCa[X1]) * (-4*get_power(mKCainf[X1],3)*(V[X1]-VK) * dmKCainf(V[X1],Ca));
  A[1][3] = (bmKCa[X2]-amKCa[X2]) * (-4*get_power(mKCainf[X2],3)*(V[X2]-VK) * dmKCainf(V[X2],Ca));
  A[2][3] = (1-bmKCa[X1]) * (-4*get_power(mKCainf[X1],3)*(V[X1]-VK) * dmKCainf(V[X1],Ca));
  A[3][3] = -get_power(mKCainf[X1],4) * (VTC-VK);

  /*
  printf("%f %f %f \n", A[0][1],A[0][2],A[0][3]);
  printf("%f %f %f \n", A[1][1],A[1][2],A[1][3]);
  printf("%f %f %f \n", A[2][1],A[2][2],A[2][3]);
  printf("%f %f %f \n", A[3][1],A[3][2],A[3][3]);
  */
  
  c[0] = BALs[X1] - ((bmNa[X1]-amNa[X1])*BAL_mNa[X1]  + (bhNa[X1]-ahNa[X1])*BAL_hNa[X1] + (bmCaT[X1]-amCaT[X1])*BAL_mCaT[X1] + (bhCaT[X1]-ahCaT[X1])*BAL_hCaT[X1] + (bmCaS[X1]-amCaS[X1])*BAL_mCaS[X1] + (bhCaS[X1]-ahCaS[X1])*BAL_hCaS[X1] + (bmKd[X1]-amKd[X1])*BAL_mKd[X1] +  (bmsyn[X1]-amsyn[X1])*newBAL_msyn[X1]);
  
  c[1] = BALs[X2] - ((bmNa[X2]-amNa[X2])*BAL_mNa[X2]  + (bhNa[X2]-ahNa[X2])*BAL_hNa[X2]  +  (bmCaT[X2]-amCaT[X2])*BAL_mCaT[X2] + (bhCaT[X2]-ahCaT[X2])*BAL_hCaT[X2] + (bmCaS[X2]-amCaS[X2])*BAL_mCaS[X2] + (bhCaS[X2]-ahCaS[X2])*BAL_hCaS[X2] + (bmKd[X2]-amKd[X2])*BAL_mKd[X2] + (bmsyn[X2]-amsyn[X2])*newBAL_msyn[X2]);
  
  c[2] = BALu[X1]  - ((1-bmNa[X1])*BAL_mNa[X1]  + (1-bhNa[X1])*BAL_hNa[X1] + (1-bmCaT[X1])*BAL_mCaT[X1] + (1-bhCaT[X1])*BAL_hCaT[X1] + (1-bmCaS[X1])*BAL_mCaS[X1] + (1-bhCaS[X1])*BAL_hCaS[X1] + (1-bmKd[X1])*BAL_mKd[X1] + (1-bmsyn[X1])*newBAL_msyn[X1]);
  
  c[3] = Istatic[X1] - (-gNa*get_power(mNainf[X1],3)*hNainf[X1]*(V[X1]-VNa) -gCaT*get_power(mCaTinf[X1],3)*hCaTinf[X1]*(V[X1]-VCa) -gCaS*get_power(mCaSinf[X1],3)*hCaSinf[X1]*(V[X1]-VCa) -gKd*get_power(mKdinf[X1],4)*(V[X1]-VK) -new_gsyn*msyninf[X1]*(V[X1]-Vsyn) -gleak*(V[X1]-Vleak));

  //printf("%f %f %f %f \n", c[0],c[1],c[2],c[3]);
  
  double Ainv[4][4];
  get_inverse(A,Ainv);

  /*
  printf("%f %f %f %f \n", Ainv[0][0],Ainv[0][1],Ainv[0][2],Ainv[0][3]);
  printf("%f %f %f %f \n", Ainv[1][0],Ainv[1][1],Ainv[1][2],Ainv[1][3]);
  printf("%f %f %f %f \n", Ainv[2][0],Ainv[2][1],Ainv[2][2],Ainv[2][3]);
  printf("%f %f %f %f \n", Ainv[3][0],Ainv[3][1],Ainv[3][2],Ainv[3][3]);
  */
  
  double gvec[4];
  int a1, a2;
  for (a1=0;a1<4; a1++)
    {
      gvec[a1] = 0;
    }
   for (a1=0;a1<4; a1++)
    {
      for (a2=0;a2<4; a2++)
	{
	  gvec[a1] += Ainv[a1][a2]*c[a2];
	}
    }
   printf("Iapp = %g \n",gvec[0]);
   printf("gH = %g \n",gvec[1]);
   printf("gA = %g \n",gvec[2]);
   printf("gKCa = %g \n",gvec[3]);
   Iapp = gvec[0];
   gH = gvec[1];
   gA = gvec[2];
   gKCa = gvec[3];
   printf("%f %f %f %f %f %f %f %f %f %f %f \n",Iapp,0.0,gNa,gCaT,gCaS,gA,gKCa,gKd,gH,gleak,new_gsyn);
   return 0;
}
