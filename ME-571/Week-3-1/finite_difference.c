/***************************
finite_difference computes an approximation to derivative of a function
using 2nd order finite difference method

dudx = (u(x+dx) - u(x-dx))/(2*dx)

Written by: Michal A. Kopera
            Department of Mathematics
            Boise State University
            1/12/2021

Based on Example 6.6 in Chopp, D.L. 'Introduction to High Performance Scientific Computing"
 **************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/*
diff comoutes the approximation to the derivative dudx to function u given by a set of discrete points

Inputs:
u  - an array of values of function u at discrete points
N  - number of points in the u array
dx - distance between points (assunes equidistant distribution)

Outputs:
dudx - contains the finite difference approximation 
*/
void diff(double* u, long N, double dx, double*dudx )
{
  long i; // local index for traversing arrays

  //at the left end-point, use one-sided difference
  dudx[0] = (u[1] - u[0])/dx; 

  // go over interior points and compute second-order finite difference
  for (i=1; i<N-1; ++i)
    {
      dudx[i] = (u[i+1] - u[i-1])/dx/2.0;
    }

  //at the right end-point compute one-sided difference
  dudx[N-1] = (u[N-1] - u[N-2])/dx; 
}

/*
init_u initializes array u with some data
Here we use a sin function

Inputs:
N  - the number of points
dx - the distance between points
 
Outputs:
u  - the values of u at specified points
dudx - the values of derivatives of u at specified points
*/

void init_u(double* u, double* dudx, double dx, long N)
{
  long i;

  //for each point compute the value of our function at i*dx
  for (i=0; i<N; ++i)
    {
      u[i] = sin(i*dx);
      dudx[i] = cos(i*dx);
    }
}

/*
compute_erro computes the L1 norm of a vector v-v_exact

Inputs:
v - solution vector
v_exact - exact solution vector
N  - the number of points in the vector
 
Outputs:
err  - error (returned through function name)
*/

double compute_error(double* v, double* v_exact, long N)
{
  long i;

  double err = 0.0;
  //for each point compute the value of our function at i*dx
  for (i=0; i<N; ++i)
    {
      err = err + fabs(v[i]-v_exact[i]);
    }

  return err;
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
  double* dudx = (double*) malloc(N*sizeof(double));
  double* dudx_exact = (double*) malloc(N*sizeof(double));
  
  // compute the interval size dx; M_PI holds value of pi
  double dx = 2.0*M_PI/(N-1);

  //start timer for initialization
  clock_t start_timer_init = clock();
  
  //create initial condition
  init_u(u,dudx_exact,dx,N);

  //stop timer for initialization
  clock_t stop_timer_init = clock();
    
  //start timer for computation
  clock_t start_timer_comput = clock();

  //call finite difference function
  diff(u,N,dx,dudx);

  //compute the error
  double error = compute_error(dudx,dudx_exact,N_loc);

  //average error per point
  error = error/N;

  //stop timer for computation
  clock_t stop_timer_comput = clock();

  //compute the elapsed time
  double elapsed_time_init = (double)(stop_timer_init - start_timer_init)/CLOCKS_PER_SEC;
  double elapsed_time_comput = (double)(stop_timer_comput - start_timer_comput)/CLOCKS_PER_SEC;
   
  //print the elapsed time
  printf("Problem size N = %ld, Error = %e, Elapsed time: initialization = %e s, computation = %e s\n",N,error,elapsed_time_init, elapsed_time_comput);

  //free allocated memory
  free(u);
  free(dudx);
  free(dudx_exact);
}
