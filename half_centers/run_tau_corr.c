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
  double tau[3]; tau[0] = 1.0; tau[1] = 10.0; tau[2] = 100.0;
  double corr[4]; corr[0] = 0.6; corr[1] = 0.7; corr[2] = 0.8; corr[3] = 0.9; //corr[4] = 0.4; corr[5] = 0.5; corr[6] = 1.0;
  
  char cmd[999];
  
  //ll = system ("g++ model_7conduct.c -o m2.out -lgsl");   
  
  int ll;
  for (int kk=0;kk<3;kk++)
    {
      for (int jj=0;jj<4;jj++)
	{
	  printf("%g %g \n",tau[kk],corr[jj]);
	  
	  //sprintf(cmd, "sbatch run_script.txt %f %f %f %f %f %f %f", mean, sigma, gNa[ii], ratio1[jj], ratio2[kk], ratio3[tt], tau_c);
	  sprintf(cmd, "./m2.out 1.1 1.0 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 0.01 0.01 %f %f 0",corr[jj],tau[kk]);
	  ll = system(cmd);

	    
	  sprintf(cmd, "./m2.out 0 1.0 227.052 2.2781 2.8469 30.4321 1121.124 75.649 1.3194 0.1631 227.052 2.2781 2.8469 30.4321 1121.124 75.649 1.3194 0.1631 0.01 0.01 %f %f 0",corr[jj],tau[kk]);
	  /*
	  sprintf(cmd, "./m2.out 0 1.0 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 0.2 0.2 %f %f 0",corr[jj],tau[kk]);
	  sprintf(cmd, "./m2.out 0 1.0 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 0.005 0.005 %f %f 0",corr[jj],tau[kk]);
	  sprintf(cmd, "./m2.out 5 1.0 227.052 2.7781 3.3469 45.0    121.124 75.649 1.3194 0.1631 227.052 2.7781 3.3469 45.0    121.124 75.649 1.3194 0.1631 0.2 0.2 %f %f 0",corr[jj],tau[kk]);
	  sprintf(cmd, "./m2.out 1.1 1.0 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 0.01 0.01 %f %f 0",corr[jj],tau[kk]);
	  sprintf(cmd, "./m3.out 0 1.0 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 227.052 2.7781 3.3469 30.4321 121.124 75.649 1.3194 0.1631 0.01 0.01 %f %f 0",corr[jj],tau[kk]);
	  */
	 
	  ll = system(cmd);
	} 
    }    
  
  return 0;	      
}
