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
#include "common.h"
#include <cuda_runtime.h>

#define THREADS_PER_BLOCK 1024

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
	    // break;
        }
    }

if (match) printf("Arrays match.\n\n");

    return;
}



void diffOnHost(float *u,const int N, float dx,float *dudx )
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


__global__ void diffOnGPU(float *u, const int N, float dx, float *dudx )
{
 int i = blockIdx.x * blockDim.x + threadIdx.x; 


  if(i==0)   //at the left end-point, use one-sided difference
  {
    dudx[i] = (u[i+1] - u[i])/dx;
    }
  else if(i==N-1)   //at the right end-point compute one-sided difference
  {
    dudx[i] = (u[i] - u[i-1])/dx;
  }
    else   // go over interior points and compute second-order finite difference
    {
      dudx[i] = (u[i+1] - u[i-1])/dx/2.0;
    }

}

__global__ void diffOnGPU_shared(float *u, const int N, float dx, float *dudx )
{

  float u_shared[THREADS_PER_BLOCK]; //declare shared array
 int i_glb = blockIdx.x * blockDim.x + threadIdx.x; 
 int i_loc = threadIdx.x;
 
 //load data into shared array
 if(i_glb<N){
   u_shared[i_loc] = u[i_glb];
 }
 
 __syncthreads(); //make sure copying to shared array is done
 
 //operate on shared array instead
  if(i_glb==0)   //at the global left end-point, use one-sided difference
  {
    dudx[i_glb] = (u[i_glb+1] - u[i_glb])/dx;
    }
  else if(i_glb==N-1)   //at the global right end-point compute one-sided difference
  {
    dudx[i_glb] = (u[i_glb] - u[i_glb-1])/dx;
  }
    else   // go over interior points and compute second-order finite difference
    {
      if(i_loc>0 && i_loc < THREADS_PER_BLOCK-1) { //if not the first and last thread in a block
	dudx[i_glb] = (u_shared[i_loc+1] - u_shared[i_loc-1])/dx/2.0;
      }else{
	dudx[i_glb] = (u[i_glb+1] - u[i_glb-1])/dx/2.0;
      }
    }

}

__global__ void diffOnGPU_shared1(float *u, const int N, float dx, float *dudx )
{

  float u_shared[THREADS_PER_BLOCK+2]; //declare shared array including ghost points
 int i_glb = blockIdx.x * blockDim.x + threadIdx.x; 
 int i_loc = threadIdx.x+1; //shift the local index to account for ghost points
 
 //load data into shared array
 if(i_glb<N){
   u_shared[i_loc] = u[i_glb];
 }

 if(threadIdx.x==0){
   u_shared[0] = u[i_glb-1]; //first thread loads the ghost from neighboring left block into ghost
 }

 if(threadIdx.x==THREADS_PER_BLOCK-1 && i_glb<N){
   u_shared[THREADS_PER_BLOCK+1] = u[i_glb+1]; //last thread loads the ghost from neighboring right block into ghost
 }
 
 __syncthreads(); //make sure copying to shared array is done
 
 //operate on shared array
 if(i_glb==0)   //at the global left end-point, use one-sided difference
   {
     dudx[i_glb] = (u[i_glb+1] - u[i_glb])/dx;
   }
 else if(i_glb==N-1)   //at the global right end-point compute one-sided difference
   {
     dudx[i_glb] = (u[i_glb] - u[i_glb-1])/dx;
   }
 else   // go over interior points and compute second-order finite difference
   {
     dudx[i_glb] = (u_shared[i_loc+1] - u_shared[i_loc-1])/dx/2.0;
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

void init_u(float *u, double dx, const int N)
{
 int i;

  //for each point compute the value of our function at i*dx
  for (i=0; i<N; ++i)
    {
      u[i] = sin(i*dx);
    }
}

__global__ void initOnGPU(float *u, double dx, const int N)
{

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i<N)
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

// set up device
    int dev = 0;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("Using Device %d: %s\n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));


  // get the number of points from input parameters
  int N = atoi(argv[1]);
  
  double istart, iElaps, istarth, iElapsh, iElap_copy, iElaps_shared; 
  
  istart = seconds();
  
  //allocate host memory
  size_t nBytes = N * sizeof(float);
  
  float *h_u,  *hostRef, *gpuRef;

  h_u     = (float *)malloc(nBytes);
  hostRef = (float *)malloc(nBytes);
  gpuRef  = (float *)malloc(nBytes);
  
  memset(hostRef, 0, nBytes);
  memset(gpuRef,  0, nBytes);
  
  
    // malloc device global memory
   float *d_u,*d_dudx;
   CHECK( cudaMalloc((float**)&d_u,nBytes));
   CHECK( cudaMalloc((float**)&d_dudx,nBytes));
  
  
  
  // compute the interval size dx; M_PI holds value of pi
  double dx = 2.0*M_PI/(N-1);



  //create initial condition on host
  init_u(h_u,dx,N);

  istarth = seconds();
  
    //call finite difference function on host
  diffOnHost(h_u,N,dx,hostRef);
 
 iElapsh = seconds() - istarth;
 
 
 // invoke kernel at host side
int iLen = THREADS_PER_BLOCK;
 dim3 block (iLen); //how many threads in a block
 dim3 grid  ((N + block.x -1) / block.x);


 
initOnGPU<<<grid,block>>>(d_u,dx,N);
CHECK(cudaDeviceSynchronize());
  /*
int nBlocks = 1;
 int nThreads = 1024; */


//time the global memory kernel
  istart = seconds();
  diffOnGPU<<<grid, block>>>(d_u,N,dx,d_dudx);
 CHECK(cudaDeviceSynchronize());
//diffOnGPU<<<nBlocks, nThreads>>>(d_u,N,dx,d_dudx);
 // copy kernel result back to host side

iElaps = seconds() - istart;

CHECK(cudaMemcpy(gpuRef, d_dudx, nBytes, cudaMemcpyDeviceToHost));
// check device results
 checkResult(hostRef, gpuRef, N);


 //time the shared memory kernel
 istart = seconds();
 //diffOnGPU_shared<<<grid, block>>>(d_u,N,dx,d_dudx);
 CHECK(cudaDeviceSynchronize());
 iElaps_shared = seconds() - istart;


 
 istart = seconds();
 
CHECK(cudaMemcpy(gpuRef, d_dudx, nBytes, cudaMemcpyDeviceToHost));
 
iElap_copy = seconds() - istart;
 
// check kernel error
    CHECK(cudaGetLastError()) ;
 
  // check device results
    checkResult(hostRef, gpuRef, N);


    printf("N = %d, threads = %d, blocks = %d,  timelapsed host = %f device_global = %f device_local = %f copy = %f\n ",N,block.x,grid.x,iElapsh,iElaps,iElaps_shared,iElap_copy);
  
  
  //free device global memory
 CHECK( cudaFree(d_u));
 CHECK( cudaFree(d_dudx));
  
  //free host memory
  free(h_u);
  free(hostRef);
  free(gpuRef);
  
  cudaDeviceReset();
 
 return(0);
}
