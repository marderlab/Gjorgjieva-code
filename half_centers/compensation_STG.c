//25 Sept 2015: this is actually not complete, many variables not defined. See compensation_STG3.c
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

// Parameters of the model: global
//p values
int p_Na = 3;
int p_A = 3;
int p_K = 4;
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
  double X = (5*Ca0)/(63*exp((5*V0)/63 + 283/126)*(Ca0 + 3)*(1/exp((5*V0)/63 + 283/126) + 1)*(1/exp((5*V0)/63 + 283/126) + 1));
  return X;
}

void compare_tau(double a[300001], double b[300001], double tau[300001], double tauf[300001], double taus[300001], double tauu[300001])
{
  for (int ii=0;ii<300001; ii++)
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
  double gleak = atof(argv[9]); //uS/nF    
  //synaptic
  double gsyn =  atof(argv[10]);
  double new_gsyn = atof(argv[11]);
  double corr_coef = 0;//atof(argv[12]); //correlation coefficient of the input
  double tau_c = 0;//atof(argv[13]); //correlation time const
  int num_run = 0;//atof(argv[14]);
  
  double Iapp = mean;
  //reversal potentials
  double VNa = 50; //mV or 30?
  double VA = -80; //mV
  double VK = -80; //mV
  double VCa = 80; //mV 
  double Vleak = -50; //mV
  double C = 0.628; //nF 
  double tauCa = 200; //ms
  //synaptic
  double Vsyn = -80; //mV
  double tausyn = 20; //ms

  char sp_file1 [999];
  sprintf (sp_file1, "Istatic.dat");
  FILE *spfile1;
  spfile1 = fopen(sp_file1,"w");
 
  //Initiate V
  //Compute the different time constants - all vectors, at each value of V
  double V[300001];

  double Ca = 10.0;
  double Istatic[300001];
  
  int nn = 0; //for determining Vosc
  double Vosc = 0;
  double Vzero[3]; Vzero[0]=0; Vzero[1]=0; Vzero[2]=0;
  for (int ii=0;ii<300001;ii++)
    {
      V[ii] = -100+ii*0.001;
       //Compute the static current (I/V curve)
      Istatic[ii] = -gNa*get_power(Xinf(V[ii],25.5,-5.29),3)*Xinf(V[ii],48.9,5.18)*(V[ii]-VNa) - gCaT*get_power(Xinf(V[ii],27.1,-7.2),3)*Xinf(V[ii],32.1,5.5)*(V[ii]-VCa) - gCaS*get_power(Xinf(V[ii],33,-8.1),3)*Xinf(V[ii],60,6.2)*(V[ii]-VCa) - gA*get_power(Xinf(V[ii],27.2,-8.7),3)*Xinf(V[ii],56.9,4.9)*(V[ii]-VK) -gKCa*get_power(m_KCainf(V[ii],Ca),4)*(V[ii]-VK) - gKd*get_power(Xinf(V[ii],12.3,-11.8),4)*(V[ii]-VK) - gleak*(V[ii]-Vleak) -gsyn*Xinf(V[ii],40,-2)*(V[ii]-Vsyn) + Iapp;

      /* Look for the zeros of the I/V curve. If there are 3 zeros, V_osc is the
	 more depolarized (the third). If there is only one zero, it corresponds
	 to V_osc */
      if (ii>0)
	{
	  if ((Istatic[ii-1]>0 && Istatic[ii]<=0) || (Istatic[ii-1]<0 && Istatic[ii]>=0))
	    {
	      printf("%d %f %f %f \n",ii-1,V[ii-1],Istatic[ii-1],Istatic[ii]);
	      Vzero[nn] = V[ii-1];
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
  double mKdinf1 = Xinf(VTC,12.3,-11.8);  double mKdinf2 = Xinf(Vosc,12.3,-11.8);
  double mAinf1 = Xinf(VTC,27.2,-8.7); double mAinf2 = Xinf(Vosc,27.2,-8.7);
  double hAinf1 = Xinf(VTC,56.9,4.9); double hAinf2 = Xinf(Vosc,56.9,4.9); 
  double mKCainf1 = m_KCainf(VTC,Ca); double mKCainf2 = m_KCainf(Vosc,Ca); 
  double mNaInf1 = Xinf(VTC,25.5,-5.29); double mNaInf2 = Xinf(Vosc,25.5,-5.29); 
  double hNaInf1 = Xinf(VTC,48.9,5.18); double hNaInf2 = Xinf(Vosc,48.9,5.18);
  double mCaTinf1 = Xinf(VTC,27.1,-7.2); double mCaTinf2 = Xinf(Vosc,27.1,-7.2);
  double hCaTinf1 = Xinf(VTC,32.1,5.5); double hCaTinf2 = Xinf(Vosc,32.1,5.5);
  double mCaSinf1 = Xinf(VTC,33,-8.1); double mCaSinf2 = Xinf(Vosc,33,-8.1);
  double hCaSinf1 = Xinf(VTC,60,6.2); double hCaSinf2 = Xinf(Vosc,60,6.2);
  double msyninf1 = Xinf(VTC,40,-2);  double msyninf2 = Xinf(Vosc,40,-2); //treat V1=V2 for the two neurons
  //treat V1=V2 for the two neurons

      
  double A[4][4];
  double b[4];
  double c[4];
  //A*b = c with b = [Iapp gKd gA gKCa]'
  A[0][0] = 0; 
  A[1][0] = 0; 
  A[2][0] = 0; 
  A[3][0] = 1; 
  //
  A[0][1] = (bmKd[X1]-amKd[X1]) * (-4*get_power(mKdinf1,3)*(VTC-VK) * dXinf(VTC,12.3,-11.8));
  A[1][1] = (bmKd[X2]-amKd[X2]) * (-4*get_power(mKdinf2,3)*(Vosc-VK) * dXinf(Vosc,12.3,-11.8));
  A[2][1] = (1-bmKd[X1]) *  (-4*get_power(mKdinf1,3)*(VTC-VK) * dXinf(VTC,12.3,-11.8));
  A[3][1] = -get_power(mKdinf1,4) * (VTC-VK);
  //
  A[0][2] = (bmA[X1]-amA[X1]) * (-3*get_power(mAinf1,2)*hAinf1*(VTC-VK) * dXinf(VTC,27.2,-8.7)) + (bhA[X1]-ahA[X1]) * (-get_power(mAinf1,3)*(VTC-VK) * dXinf(VTC,56.9,4.9));
  A[1][2] = (bmA[X2]-amA[X2]) * (-3*get_power(mAinf2,2)*hAinf2*(Vosc-VK) * dXinf(Vosc,27.2,-8.7)) + (bhA[X2]-ahA[X2]) * (-get_power(mAinf2,3)*(Vosc-VK) * dXinf(Vosc,56.9,4.9));
  A[2][2] = (1-bmA[X1]) * (-3*get_power(mAinf1,2)*hAinf1*(VTC-VK) * dXinf(VTC,27.2,-8.7)) + (1-bhA[X1]) *  (-get_power(mAinf1,3)*(VTC-VK) * dXinf(VTC,56.9,4.9));
  A[3][2] = -get_power(mAinf1,3)*hAinf1*(VTC-VK);
  //
  A[0][3] = (bmKCa[X1]-amKCa[X1]) * (-4*get_power(mKCainf1,3)*(VTC-VK) * dmKCainf(VTC,Ca));
  A[1][3] = (bmKCa[X2]-amKCa[X2]) * (-4*get_power(mKCainf2,3)*(Vosc-VK) * dmKCainf(Vosc,Ca));
  A[2][3] = (1-bmKCa[X1]) * (-4*get_power(mKCainf1,3)*(VTC-VK) * dmKCainf(VTC,Ca));
  A[3][3] = -get_power(mKCainf1,4) * (VTC-VK);

   double BAL_mNa1 = gNa * (-3*get_power(mNainf1,2)*hNainf1*(VTC-VNa) * dXinf(V[ii],25.5,-5.29));
   double BAL_mNa2 = gNa * (-3*get_power(mNainf2,2)*hNainf2*(Vosc-VNa) * dXinf(V[ii],25.5,-5.29));
   double BAL_hNa1 = gNa * (-get_power(mNainf1,3)*(VTC-VNa) * dXinf(VTC,48.9,5.18));
   double BAL_hNa2 = gNa * (-get_power(mNainf2,3)*(Vosc-VNa) * dXinf(Vosc,48.9,5.18));
   double BAL_mCaT1 = gCaT * (-3*get_power(mCaTinf1,2)*hCaTinf1*(VTC-VCa) * dXinf(VTC,27.1,-7.2));
   double BAL_mCaT2 = gCaT * (-3*get_power(mCaTinf2,2)*hCaTinf2*(Vosc-VCa) * dXinf(Vosc,27.1,-7.2));
   double BAL_hCaT1 = gCaT * (-get_power(mCaTinf1,3)*(VTC-VCa) * dXinf(VTC,32.1,5.5));
   double BAL_hCaT2 = gCaT * (-get_power(mCaTinf2,3)*(Vosc-VCa) * dXinf(Vosc,32.1,5.5));
   double BAL_mCaS1 = gCaS * (-3*get_power(mCaSinf1,2)*hCaSinf1*(VTC-VCa) * dXinf(VTC,33,-8.1));
   double BAL_mCaS2 = gCaS * (-3*get_power(mCaSinf2,2)*hCaSinf2*(Vosc-VCa) * dXinf(Vosc,33,-8.1));
   double BAL_hCaS1 = gCaS * (-get_power(mCaSinf1,3)*(VTC-VCa) * dXinf(VTC,60,6.2));
   double BAL_hCaS2 = gCaS * (-get_power(mCaSinf2,3)*(Vosc-VCa) * dXinf(Vosc,60,6.2));
   double BAL_mA1 = gA * (-3*get_power(mAinf1,2)*hAinf1*(VTC-VK) * dXinf(VTC,27.2,-8.7));
   double BAL_mA2 = gA * (-3*get_power(mAinf2,2)*hAinf2*(Vosc-VK) * dXinf(Vosc,27.2,-8.7));
   double BAL_hA1 = gA * (-get_power(mAinf1,3)*(VTC-VK) * dXinf(Vosc,56.9,4.9));
   double BAL_hA2 = gA * (-get_power(mAinf2,3)*(VTC-VK) * dXinf(Vosc,56.9,4.9));
   double BAL_mKd1 = gKd * (-4*get_power(mKdinf1,3)*(VTC-VK) * dXinf(VTC,12.3,-11.8));
   double BAL_mKd2 = gKd * (-4*get_power(mKdinf2,3)*(Vosc-VK) * dXinf(Vosc,12.3,-11.8));
   double BAL_mKCa1 = gKCa*(-4*get_power(mKCainf1,3)*(VTC-VK) * dmKCainf(VTC,Ca));
   double BAL_mKCa2 = gKCa*(-4*get_power(mKCainf2,3)*(Vosc-VK) * dmKCainf(Vosc,Ca));
   //may need to make the next a minus sign!
   double BAL_msyn1 = gsyn*(-dXinf(VTC,40,-2)*(VTC-Vsyn)); //treat V1=V2 for the two neurons
   double BAL_msyn2 = gsyn*(-dXinf(Vosc,40,-2)*(Vosc-Vsyn)); //treat V1=V2 for the two neurons
   double newBAL_msyn1 = new_gsyn*(-dXinf(VTC,40,-2)*(VTC-Vsyn)); //treat V1=V2 for the two neurons
   double newBAL_msyn2 = new_gsyn*(-dXinf(Vosc,40,-2)*(Vosc-Vsyn)); //treat V1=V2 for the two neurons


   BALf[ii] = amNa[ii]*BAL_mNa[ii] + ahNa[ii]*BAL_hNa[ii] + amCaT[ii]*BAL_mCaT[ii] + ahCaT[ii]*BAL_hCaT[ii] + amCaS[ii]*BAL_mCaS[ii] + ahCaS[ii]*BAL_hCaS[ii] + amA[ii]*BAL_mA[ii] + ahA[ii]*BAL_hA[ii] + amKd[ii]*BAL_mKd[ii] + amKCa[ii]*BAL_mKCa[ii] + amsyn[ii]*BAL_msyn[ii];
      
      BALs[ii] = (bmNa[ii]-amNa[ii])*BAL_mNa[ii] + (bhNa[ii]-ahNa[ii])*BAL_hNa[ii] + (bmCaT[ii]-amCaT[ii])*BAL_mCaT[ii] + (bhCaT[ii]-ahCaT[ii])*BAL_hCaT[ii] + (bmCaS[ii]-amCaS[ii])*BAL_mCaS[ii] + (bhCaS[ii]-ahCaS[ii])*BAL_hCaS[ii] + (bmA[ii]-amA[ii])*BAL_mA[ii] + (bhA[ii]-ahA[ii])*BAL_hA[ii] + (bmKd[ii]-amKd[ii])*BAL_mKd[ii] + (bmKCa[ii]-amKCa[ii])*BAL_mKCa[ii] + (bmsyn[ii]-amsyn[ii])*BAL_msyn[ii];

      BALu[ii] = (1-bmNa[ii])*BAL_mNa[ii] + (1-bhNa[ii])*BAL_hNa[ii] + (1-bmCaT[ii])*BAL_mCaT[ii] + (1-bhCaT[ii])*BAL_hCaT[ii] + (1-bmCaS[ii])*BAL_mCaS[ii] + (1-bhCaS[ii])*BAL_hCaS[ii] + (1-bmA[ii])*BAL_mA[ii] + (1-bhA[ii])*BAL_hA[ii] + (1-bmKd[ii])*BAL_mKd[ii] + (1-bmKCa[ii])*BAL_mKCa[ii] + (1-bmsyn[ii])*BAL_msyn[ii];

      
  c[0] = BALs[X1] - (bmNa[X1]-amNa[X1])*BAL_mNa1  + (bhNa[X1]-ahNa[X1])*BAL_hNa1 + (bmCaT[X1]-amCaT[X1])*BAL_mCaT1 + (bhCaT[X1]-ahCaT[X1])*BAL_hCaT1+ (bmCaS[X1]-amCaS[X1])*BAL_mCaS1+ (bhCaS[X1]-ahCaS[X1])*BAL_hCaS1 +  (bmsyn[X1]-amsyn[X1])*newBAL_msyn1;
  
  c[1] = BALs[X2] - (bmNa[X2]-amNa[X2])*BAL_mNa[X2]  + (bhNa[X2]-ahNa[X2])*BAL_hNa[X2]  +  (bmCaT[X2]-amCaT[X2])*BAL_mCaT[X2] + (bhCaT[X2]-ahCaT[X2])*BAL_hCaT[X2] + (bmCaS[X2]-amCaS[X2])*BAL_mCaS[X2] + (bhCaS[X2]-ahCaS[X2])*BAL_hCaS[X2] + (bmsyn[X2]-amsyn[X2])*newBAL_msyn[X2];
  
  c[2] = BALu[X1]  - (1-bmNa[X1])*BAL_mNa[X1]  + (1-bhNa[X1])*BAL_hNa[X1] + (1-bmCaT[X1])*BAL_mCaT[X1] + (1-bhCaT[X1])*BAL_hCaT[X2] + (1-bmCaS[X1])*BAL_mCaS[X1] + (1-bhCaS[X1])*BAL_hCaS[X1] + (1-bmsyn[X1])*newBAL_msyn[X1];
  
  c[3] = Istatic[X1] - (-gNa*get_power(mNainf[X1],3)*hNainf[X1]*(V[X1]-VNa) -gCaT*get_power(mCaTinf[X1],3)*hCaTinf[X1]*(V[X1]-VCa) -gCaS*get_power(mCaSinf[X1],3)*hCaSinf[X1]*(V[X1]-VCa) -new_gsyn*msyninf[X1]*(V[X1]-Vsyn) -gleak*(V[X1]-Vleak));

  double Ainv[4][4];
  get_inverse(A,Ainv);
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
	  gvec[a1] += Ainv[a2][a1]*c[a2];
	}
    }
   printf("Iapp = %g \n",gvec[0]);
   printf("gKd = %g \n",gvec[1]);
   printf("gA = %g \n",gvec[2]);
   printf("gKCa = %g \n",gvec[3]);

   fclose(spfile1);
  return 0;
}
