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
#include <mpi.h> //Include MPI library

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

void init_u(double* u, double* dudx, double dx, long N_loc, int irank)
{
  long i,i_loc;

  //for each point compute the value of our function at i*dx
  for (i_loc=0; i_loc<N_loc; ++i_loc)
    {
      //relation between local and global index
      i = N_loc*irank + i_loc;

      u[i_loc] = sin(i*dx);
      dudx[i_loc] = cos(i*dx);
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

  //Initialize MPI - create all available processes
  MPI_Init(&argc, &argv);

  //Variables to store communicator size and rank number
  int nproc, irank; 
  double wtime, wtime_global;
  //Check communicator size
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  //Check rank number
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);

  // get the number of points from input parameters
  long N = atol(argv[1]);

  // compute how many points each rank has
  long N_loc = N/nproc;
  
  //allocate arrays
  double* u    = (double*) malloc(N_loc*sizeof(double));
  double* dudx = (double*) malloc(N_loc*sizeof(double));
  double* dudx_exact = (double*) malloc(N_loc*sizeof(double));

  //start time loop                                                                
  wtime = MPI_Wtime();
  
  // compute the interval size dx; M_PI holds value of pi
  double dx = 2.0*M_PI/(N-1);
    
  //create initial condition
  init_u(u,dudx_exact,dx,N_loc,irank);

  //call finite difference function
  diff(u,N_loc,dx,dudx);

  wtime = MPI_Wtime()-wtime;

  MPI_Reduce(&wtime, &wtime_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(irank == 0){
    printf("%d,%f\n",N,wtime_global);
  }


  //compute the error
  //  double error = compute_error(dudx,dudx_exact,N_loc);

  // sum up errors from all processes

  //average error per point
  //  error = error/N_loc;
 
  //printf("%d,%f\n",N,elapsed_time_comput);
  //free allocated memory
  free(u);
  free(dudx);
  free(dudx_exact);

  //Finalize MPI - close all ranks
  MPI_Finalize();

}
