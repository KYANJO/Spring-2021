#include <cuda_runtime.h>
#include <curand.h>
#include "common.h"
#include <stdio.h>
#include <curand_kernel.h>

#define THREADS_PER_BLOCK 256

__global__ void reductionOnGPU(float *f, const int N, const int nb, curandState *states)
{   
	__shared__ double u[THREADS_PER_BLOCK];
	//extern __shared__ float u[];
	int i_glb = blockIdx.x * blockDim.x + threadIdx.x; 
	unsigned int i_loc = threadIdx.x;
	int ib = blockDim.x; 
	//unsigned int i_locb = blockIdx.x;
	unsigned int i;

	int n = N/nb;
	//double cc = 0;
	//initialse curand
	curand_init(1234, i_glb, 0, &states[i_glb]);
	//curandState state = states[ix];

	if (i_loc<n)
	{
		double xran = curand_uniform_double (&states[i_glb]);
		u[i_loc] = cos(-log(xran));
		__syncthreads();
		//u[i_loc] = cc;
	}
	//u[i_loc] = cc;
	//__syncthreads();
	
	for (i = 1; i<ib; i *=2)	
	{
		__syncthreads();
		if(i_loc % (2*i) == 0)
		{
			u[i_loc] += u[i + i_loc];
		}
	}
	__syncthreads();
	//Compute a global sum
	if(i_loc==0)
	{
		//f[blockIdx.x] = u[0];
		atomicAdd(f,u[0]); 
		//printf("%f\n", u[0]);
	}
	__syncthreads();
	//printf("%f\n", u[0]);
					
}

__global__ void integralOnGPU(float *f,double *Int ,const int N)
{
	*Int = abs(*f/N);
}


int main(int argc, char **argv)
{	
	// problem size
	int N = atoi(argv[1]);
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
	//int B = ((N + T -1)/T);
	//if(B > 65535) B = 65535;
	int B = 64;
	int nb = N/(B*T);
	
	//states allocate memory
	CHECK(cudaMalloc( (void **)&States, (B*T)*sizeof(curandState)));

	reductionOnGPU<<<B,T>>>(d_f, N, nb,States);
	CHECK(cudaDeviceSynchronize());

	integralOnGPU<<<1,1>>>(d_f,d_Int ,N);
	CHECK(cudaDeviceSynchronize());

	 // check kernel error
	 CHECK(cudaGetLastError());

	 double iElaps_g = seconds() - iStart;
	
	// copy kernel result back to host side
    CHECK(cudaMemcpy(&gpuRef, d_Int, sizeof(double), cudaMemcpyDeviceToHost));

	//error achived
	double error = 0.5 - gpuRef;
	printf("%d,%f,%e,%f\n",N,gpuRef,error,iElaps_g);
	
	//free device memory
	CHECK(cudaFree(States));
	CHECK(cudaFree(d_f));
    CHECK(cudaFree(d_Int));

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
