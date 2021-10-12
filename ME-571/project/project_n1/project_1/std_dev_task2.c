/***************************                                                    
std_dev reads problem size N as input argument and computes a mean a\
nd standard deviation.                                                        

Written by: Brian Kyanjo                                                    
            Department of Mathematics                                           
            Boise State University                                              
            2/15/2021                                                                                                                                         **************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//Data array using random number                                               
double random_number()
{
  double random_value;
  random_value = (double)rand()/(double)RAND_MAX;
  return random_value;
}

void main(int argc, char* argv[]){

  long N = atol(argv[1]);

  //array of random numbers
  double* x =  (double*) malloc(N*sizeof(double));
  int i;
  
  for (i = 0; i < N; ++i)
  {
    x[i] = (double)random_number();
  }

  // Welford's algorithm                       
  double* xbar =  (double*) malloc(N*sizeof(double));
  double xbar_o = x[0]; xbar[0] = xbar_o;
  double* M  =  (double*) malloc(N*sizeof(double));
  int k;
  double sigma;

   for (k = 1; k<=N-1; ++k)
   {
      xbar[k] = xbar[k-1] + (x[k] - xbar[k-1])/k;
      M[k] = M[k-1] + (x[k] - xbar[k])*(x[k] - xbar[k-1]);
   }

  //compute the standard deviation                                              
  sigma = sqrt(M[N-1]/N);

  printf("sigma = %e\n",sigma);

}
