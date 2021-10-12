#include <cuda_runtime.h>
#include <curand.h>
#include "common.h"
#include <stdio.h>
#include <curand_kernel.h>

#define THREADS_PER_BLOCK 64

__global__ void mcrandom(double *udata, const int N, const int nb, curandState *states)
{
	unsigned int i_glb = blockIdx.x * blockDim.x + threadIdx.x; 
	int n = N/nb;
	//initialse curand
	//double cc = 0;
	curand_init(1234, i_glb, 0, &states[i_glb]);
	if (i_glb<n)
	{
		for(int i=0; i<nb; i++)
		{
		double xran = curand_uniform_double (&states[i_glb]);
		udata[i_glb] += cos(-log(xran));
		}
		
	}

	//udata[i_glb] = cc;
	
}

__global__ void reductionOnGPU(float *f,const int N, const int nb, curandState *states)
{   
	__shared__ double u[THREADS_PER_BLOCK];
	//extern __shared__ double u[];
	unsigned int i_glb = blockIdx.x * blockDim.x + threadIdx.x; 
	unsigned int i_loc = threadIdx.x;
	int ib = blockDim.x; 
	

	int n = N/nb;
	//initialse curand
	double cc = 0;
	curand_init(1234, i_glb, 0, &states[i_glb]);
	if (i_glb<n)
	{
		for(int i=0; i<nb; i++)
		{
		double xran = curand_uniform_double (&states[i_glb]);
		cc += cos(-log(xran));
		}
		
	}

	//load memory
	u[i_loc] = cc;
	
	__syncthreads();
	
	//reduction in shared memory
	unsigned int i;
	for (i = 1; i<ib; i *=2)	
	{
		if(i_loc % (2*i) == 0)
		{
			//__syncthreads();
			u[i_loc] += u[i + i_loc];
		}
		__syncthreads();
	}
	//printf("%f\n",u[i_glb]);
/*

	for (i = 1; i<ib; i *=2)	
	{
		int index = 2*i*i_loc;
		if (index < blockDim.x) 
		{
			u[index] += u[i + index];
		}
		__syncthreads();
	}
	*/
	//append top the a global memory
	if(i_loc==0)
	{
		//f[blockIdx.x] = u[0];
		__syncthreads();

		atomicAdd(f,u[0]/N);
		//__syncthreads();
	}
	//printf("%f,%d\n",*f,N);
	
}

__global__ void integralOnGPU(float *f, double *Int ,const int N)
{
	printf("%f,%d\n",*f,N);
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
	float gpuRef;
	
	//start timing
	double iStart = seconds();

	// malloc device global memory
	//float *d_f;
	float *d_Int;
	//double *d_udata;
	//float *d_c;
   
	CHECK(cudaMalloc((void **)&d_Int, sizeof(float)));
	
	//invoke the kernel
	int B = ((N + T -1)/T);
	if(B > 65535) B = 65535;
	
	int nb = ceil((N*1.0)/(B*T));
	
	//states allocate memory
	CHECK(cudaMalloc( (void **)&States, (B*T)*sizeof(curandState)));
	//CHECK(cudaMalloc((void **)&d_udata, (B*T)*sizeof(double)));
	//CHECK(cudaMalloc((void **)&d_f, N*sizeof(float)));
    
	//mcrandom<<<B,T>>>(d_udata, N, nb, States);
	//CHECK(cudaDeviceSynchronize());
	reductionOnGPU<<<B,T>>>(d_Int, N, nb, States);
	CHECK(cudaDeviceSynchronize());

	//integralOnGPU<<<1,1>>>(d_f,d_Int ,N);
	//CHECK(cudaDeviceSynchronize());

	 // check kernel error
	CHECK(cudaGetLastError());

	 double iElaps_g = seconds() - iStart;
	
	// copy kernel result back to host side
    CHECK(cudaMemcpy(&gpuRef, d_Int, sizeof(float), cudaMemcpyDeviceToHost));

	//error achived
	double error = 0.5 - gpuRef;
	//printf("%f\n",gpuRef);
	printf("%lld,%f,%e,%f\n",N,gpuRef,error,iElaps_g);
	
	//free device memory
	CHECK(cudaFree(States));
	//CHECK(cudaFree(d_f));
    CHECK(cudaFree(d_Int));
	//CHECK(cudaFree(d_udata));
	//CHECK(cudaFree(d_c));

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
