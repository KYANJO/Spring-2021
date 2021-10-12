/***************************
trapez_integration computes an approximation of a definite integral of a function
given by a set of discrete points using trapezoidal rule

I = sum(0.5*dx*(u(i+1)+u(i)))

Written by: Michal A. Kopera
            Department of Mathematics
            Boise State University
            1/25/2021

 **************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*
trapez comoutes the approximation to the definite integral I to function u given by a set of discrete points

Inputs:
u  - an array of values of function u at discrete points
N  - number of points in the u array
dx - distance between points (assunes equidistant distribution)

Outputs:
I  - the value of the integral (returned through function name)
*/
double trapez(double* u, long N, double dx)
{
  long i; // local index for traversing arrays

  //initialize the integral variable I
  double I = 0;

  // go over interior points and compute second-order finite difference
  for (i=0; i<N-1; ++i)
    {
      I =  I + 0.5*dx*(u[i+1] + u[i]);
    }

  return I;
}

/*
init_u initializes array u with some data
Here we use a 4.0/(1+x*x) function

Inputs:
N  - the number of points
dx - the distance between points
 
Outputs:
u  - the values of u at specified points
*/

void init_u(double* u, double dx, long N)
{
  long i;
  double x;
  
  //for each point compute the value of our function at i*dx
  for (i=0; i<N; ++i)
    {
      x = i*dx;
      u[i] = 4.0/(1+x*x);
    }
}

/* main function 
Creates the function which derivative is computed and measures the execution time of the finite difference computation. 

Assumes that it receives N as the first (after the executable name) command line argument
Inputs:
N - number og points used to approximate function u

*/

int main(int argc, char* argv[])
{

  // get the number of points from input parameters
  long N = atol(argv[1]);
  
  //allocate arrays
  double* u    = (double*) malloc(N*sizeof(double));
  
  // compute the interval size dx - we integrate from 0 to 1
  double dx = 1.0/(N-1);

  // store the value of exact solution
  double pi_exact = M_PI;

  //start timer for initialization
  clock_t start_timer_init = clock();
  
  //create initial condition
  init_u(u,dx,N);

  //stop timer for initialization
  clock_t stop_timer_init = clock();
    
  //start timer for computation
  clock_t start_timer_comput = clock();

  //call finite difference function
  double I = trapez(u,N,dx);

  //stop timer for computation
  clock_t stop_timer_comput = clock();

  //compute the elapsed time
  double elapsed_time_init = (double)(stop_timer_init - start_timer_init)/CLOCKS_PER_SEC;
  double elapsed_time_comput = (double)(stop_timer_comput - start_timer_comput)/CLOCKS_PER_SEC;

  //compute the error
  double pi_error = abs(pi_exact-I);
  
  //print the elapsed time
  printf("Problem size N = %ld, I = %e,  Elapsed time: initialization = %e s, computation = %e s\n",N,I,elapsed_time_init, elapsed_time_comput);

  //free allocated memory
  free(u);
}
