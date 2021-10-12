/***************************
finite_difference computes an approximation to derivative of a function
using 2nd order finite difference method

dudx = (u(x+dx) - u(x-dx))/(2*dx)

Written by: Michal A. Kopera
            Department of Mathematics
            Boise State University
Modified by:Brian Kyanjo
            1/12/2021

Based on Example 6.6 in Chopp, D.L. 'Introduction to High Performance Scientific Computing"
 **************************/
#include <stdio.h>
#include <cuda_runtime.h>

/*
diff comoutes the approximation to the derivative dudx to function u given by a set of discrete points

Inputs:
u  - an array of values of function u at discrete points
N  - number of points in the u array
dx - distance between points (assunes equidistant distribution)

Outputs:
dudx - contains the finite difference approximation 
*/

void checkResult(float *hostRef, float *gpuRef, const int N)
{
  double epsilon = 1.0E-8;
  bool match = 1;

  for (int i = 0; i < N; i++)
    {
      if (abs(hostRef[i] - gpuRef[i]) > epsilon)
        {
	  match = 0;
	  printf("Arrays do not match!\n");
	  printf("host %5.2f gpu %5.2f at current %d\n", hostRef[i],
		 gpuRef[i], i);
	  break;
        }
    }

  if (match) printf("Arrays match.\n\n");

  return;
}

void diffOnHost(float *u, const int N, double dx, float *dudx )
{
  int i; // local index for traversing arrays

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

__global__ void diffOnGPU(float *u, int N, double dx, float *dudx )
{
  int i = blockIdx.x * blockDim.x + threadIdx.x; // local index for traversing arrays                                      

  //at the left end-point, use one-sided difference                               
  if(i==0){
  dudx[i] = (u[1] - u[0])/dx;
 }
 else if(i<N-1){
  // go over interior points and compute second-order finite difference           
  dudx[i] = (u[i+1] - u[i-1])/dx/2.0;
  }
  else if(i==N-1){
  //at the right end-point compute one-sided difference                            
  dudx[i] = (u[i] - u[i-1])/dx;
}
}

/*
init_u initializes array u with some data
Here we use a sin function

Inputs:
N  - the number of points
dx - the distance between points
 
Outputs:
u  - the values of u at specified points
*/

void init_u(float* u, double dx, int N)
{
  int i;

  //for each point compute the value of our function at i*dx
  for (i=0; i<N; ++i)
    {
      u[i] = sin(i*dx);
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
  printf("%s Starting...\n", argv[0]);

  // get the number of points from input parameters                   
 int nElem = atol(argv[1]);

  printf("Grid size %d\n", nElem);

  // malloc host memory
  size_t nBytes = nElem * sizeof(float);

  float *h_u, *hostRef, *gpuRef;
  h_u     = (float *)malloc(nBytes);
  hostRef = (float *)malloc(nBytes);
  gpuRef  = (float *)malloc(nBytes);

// compute the interval size dx; M_PI holds value of pi
  double dx = 2.0*M_PI/(nElem-1);

  // initialize data 
  init_u(h_u,dx, nElem);

  memset(hostRef, 0, nBytes);
  memset(gpuRef,  0, nBytes);

  // malloc device global memory
  float *d_u, *d_dudx;
  cudaMalloc((float**)&d_u, nBytes);
  cudaMalloc((float**)&d_dudx, nBytes);
 
  // transfer data from host to device
  cudaMemcpy(d_u, h_u, nBytes, cudaMemcpyHostToDevice);
  cudaMemcpy(d_dudx, gpuRef, nBytes, cudaMemcpyHostToDevice);

  int nBlocks = 2;
  int nThreads = 1024;

  //call finite difference function on the device
  diffOnGPU<<<nBlocks, nThreads>>>(d_u,nElem,dx,d_dudx);
  printf("Execution configure <<<%d, %d>>>\n", nBlocks, nThreads);

  // copy kernel result back to host side
 // cudaMemcpy(gpuRef, d_u, nBytes, cudaMemcpyDeviceToHost);
  cudaMemcpy(gpuRef, d_dudx, nBytes, cudaMemcpyDeviceToHost);

  //call finite difference function at the host side for results check                                 
  diffOnHost(h_u,nElem,dx,hostRef);

  // check device results
  checkResult(hostRef, gpuRef, nElem);

  // free device global memory
  cudaFree(d_u);
  cudaFree(d_dudx);
 
  // free host memory
  free(h_u);
  free(hostRef);
  free(gpuRef);

  cudaDeviceReset();
  return(0);

}
