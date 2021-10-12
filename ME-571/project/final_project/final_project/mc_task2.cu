/*
Author: Brian KYANJO
Date:   May 1st, 2021
Class:  ME571

Description:
------------
Monte Carlo integration implementation using CUDA
*/

#include <cuda_runtime.h>
#include <curand.h>
#include "common.h"
#include <stdio.h>
#include <curand_kernel.h>

#define THREADS_PER_BLOCK 128
#define SEED 60

__global__ void mcOnGPU(float *f, const int N, const int nb, curandState *states)
{   
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
	int n = N/nb;
	
	double cc = 0;
	//initialse curand
	curand_init((SEED << 20) + ix, 0, 0, &states[ix]);
	//curandState state = states[ix];
	if (ix<n)
	{	
		for(int i=0; i<nb; i++)
		{
			double xran = curand_uniform_double (&states[ix]);
			cc += cos(-log(xran));
		}
	}
		atomicAdd(f,cc); 
		__syncthreads();				
}

__global__ void integralOnGPU(float *f,double *Int ,const int N)
{
	Int[0] = abs(f[0]/N);
}


int main(int argc, char **argv)
{	
	// problem size
	long long int N = atol(argv[1]);
	int T = THREADS_PER_BLOCK;
	
	//random number generator
	curandState *States;
	
	// malloc host memory
	double gpuRef;
	
	//start timing
	double iStart = seconds();

	// malloc device global memory
    float *d_f;
    double *d_Int;
    CHECK(cudaMalloc((void **)&d_f, sizeof(double)));
    CHECK(cudaMalloc((void **)&d_Int, sizeof(double)));
    
	//invoke the kernel
	int B = ((N + T -1)/T);
	if(B > 65535) B = 65535;

	int nb = ceil((N*1.0)/(B*T));
	
	//states allocate memory
	CHECK(cudaMalloc( (void **)&States, (B*T)*sizeof(curandState)));

	// start kernel time
	double iStartc = seconds();
	mcOnGPU<<<B,T>>>(d_f, N, nb, States);
	CHECK(cudaDeviceSynchronize());

	integralOnGPU<<<1,1>>>(d_f,d_Int ,N);
	CHECK(cudaDeviceSynchronize());

	//end kernel time
	double iElaps_c = seconds() - iStartc;
	 // check kernel error
	CHECK(cudaGetLastError());
	
	// copy kernel result back to host side
    CHECK(cudaMemcpy(&gpuRef, d_Int, sizeof(double), cudaMemcpyDeviceToHost));

	double iElaps_s = seconds() - iStart;
	
	//error achived
	double error = abs(0.5 - gpuRef);
	//printf("%lld,%f,%e,%f\n",N,gpuRef,error,iElaps_g);
	printf("%lld,%f,%e,%f,%f\n",N,gpuRef,error,iElaps_s,iElaps_c);
	
	//free device memory
	CHECK(cudaFree(States));
	CHECK(cudaFree(d_f));
    CHECK(cudaFree(d_Int));

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
