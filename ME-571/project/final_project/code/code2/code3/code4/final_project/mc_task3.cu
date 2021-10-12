#include <cuda_runtime.h>
#include <curand.h>
#include "common.h"
#include <stdio.h>
#include <curand_kernel.h>

#define THREADS_PER_BLOCK 32
#define SEED 60
//Generate data
__global__ void mcrandom(double *udata, const int N, const int nb, curandState *states)
{
	unsigned int i_glb = blockIdx.x * blockDim.x + threadIdx.x; 
	int n = N/nb;
	
	//initialse curand
	curand_init((SEED << 20) + i_glb, 0, 0, &states[i_glb]);
	if (i_glb<n)
	{
		for(int i=0; i<nb; i++)
		{
		double xran = curand_uniform_double (&states[i_glb]);
		udata[i_glb] += cos(-log(xran));
		}
		
	}
	
}

//reduction kernel
__global__ void reductionOnGPU(double *udata,float *f)
{   
	__shared__ double u[THREADS_PER_BLOCK];
	unsigned int i_glb = blockIdx.x * blockDim.x + threadIdx.x; 
	unsigned int i_loc = threadIdx.x;
	int ib = blockDim.x; 
	unsigned int i;

	//load memory
	u[i_loc] = udata[i_glb];
	__syncthreads();

	//reduction in shared memory
	for (i = 1; i<ib; i *=2)	
	{
		int index = 2*i*i_loc;
		//__syncthreads();
		if (index < blockDim.x) 
		{
			u[index] += u[i + index];
		}
		__syncthreads();
	}

	if(i_loc==0)
	{
		atomicAdd(f,u[0]);
	}
	
}

__global__ void integralOnGPU(float *f, double *Int ,const int N)
{
	//global mean
	*Int = *f/N;
}


int main(int argc, char **argv)
{	
	// problem size
	long int N = atol(argv[1]);
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
	double *d_udata;
   
	CHECK(cudaMalloc((void **)&d_Int, sizeof(double)));
	
	//invoke the kernel
	int B = ((N + T -1)/T);
	if(B > 65535) B = 65535;
	
	int nb = ceil((N*1.0)/(B*T));
	
	//states allocate memory
	CHECK(cudaMalloc( (void **)&States, (B*T)*sizeof(curandState)));
	CHECK(cudaMalloc((void **)&d_udata, (B*T)*sizeof(double)));
	CHECK(cudaMalloc((void **)&d_f, (B*T)*sizeof(float)));
	
	//double iStart = seconds();
	mcrandom<<<B,T>>>(d_udata, N, nb, States);
	CHECK(cudaDeviceSynchronize());
	
	reductionOnGPU<<<B,T>>>(d_udata,d_f);
	CHECK(cudaDeviceSynchronize());
	

	integralOnGPU<<<1,1>>>(d_f,d_Int ,N);
	CHECK(cudaDeviceSynchronize());

	double iElaps_g = seconds() - iStart;
	
	 // check kernel error
	CHECK(cudaGetLastError());

	 //double iElaps_g = seconds() - iStart;
	
	// copy kernel result back to host side
    CHECK(cudaMemcpy(&gpuRef, d_Int, sizeof(double), cudaMemcpyDeviceToHost));

	//error achived
	double error = abs(0.5 - gpuRef);
	printf("%lld,%f,%e,%f\n",N,gpuRef,error,iElaps_g);
	
	//free device memory
	CHECK(cudaFree(States));
	CHECK(cudaFree(d_f));
    CHECK(cudaFree(d_Int));
	CHECK(cudaFree(d_udata));

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
