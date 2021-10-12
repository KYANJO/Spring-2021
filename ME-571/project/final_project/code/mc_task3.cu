#include <cuda_runtime.h>
#include <curand.h>
#include "common.h"
#include <stdio.h>
#include <curand_kernel.h>

#define THREADS_PER_BLOCK 256

__global__ void mcrandom(double *udata, const int N, const int nb, curandState *states)
{
	unsigned int i_glb = blockIdx.x * blockDim.x + threadIdx.x; 
	int n = N/nb;
	//initialse curand
	double cc = 0;
	curand_init(1234, i_glb, 0, &states[i_glb]);
	if (i_glb<n)
	{
		double xran = curand_uniform_double (&states[i_glb]);
		cc += cos(-log(xran));
	}
	udata[i_glb] = cc;
}

__global__ void reductionOnGPU(double *udata, float *f, float *c)
{   
	__shared__ double u[THREADS_PER_BLOCK];
	//extern __shared__ double u[];
	unsigned int i_glb = blockIdx.x * blockDim.x + threadIdx.x; 
	unsigned int i_loc = threadIdx.x;
	int ib = blockDim.x; 
	unsigned int i;

	//load memory
	u[i_loc] = udata[i_glb];
	__syncthreads();

	/*
	//reduction in shared memory
	for (i = 1; i<ib; i *=2)	
	{
		if(i_loc % (2*i) == 0)
		{
			//__syncthreads();
			u[i_loc] += u[i + i_loc];
		}
		__syncthreads();
	}
*/

	for (i = 1; i<ib; i *=2)	
	{
		int index = 2*i*i_loc;
		if (index < blockDim.x) 
		{
			u[index] += u[i + index];
		}
		__syncthreads();
	}
	
	//append top the a global memory
	if(i_loc==0)
	{
		c[blockIdx.x] = u[0];
		//printf("%f\n",*f);
	}
	//printf("%f\n",*c);
	//atomicAdd(f,c);

	
	/*
	for (i = 1; i<ib; i *=2)	
	{
		int index = 2*i*i_loc;
		__syncthreads();
		if (index < blockDim.x) 
		{
			c[index] += c[i + index];
		}
		__syncthreads();
	}
		
	//if(i_loc==0)
	//{
		f[blockIdx.x] = c[0];
		//printf("%f\n",*f);
	//}
    */
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
	double *d_udata;
	float *d_c;
    CHECK(cudaMalloc((void **)&d_f, sizeof(double)));
	CHECK(cudaMalloc((void **)&d_Int, sizeof(double)));
	CHECK(cudaMalloc((void **)&d_udata, T));
	CHECK(cudaMalloc((void **)&d_c, T));
	
	//invoke the kernel
	//int B = ((N + T -1)/T);
	//if(B > 65535) B = 65535;
	int B = 64;
	int nb = N/(B*T);
	
	//states allocate memory
	CHECK(cudaMalloc( (void **)&States, (B*T)*sizeof(curandState)));
    
	mcrandom<<<B,T>>>(d_udata, N, nb, States);
	CHECK(cudaDeviceSynchronize());
	reductionOnGPU<<<B,T>>>(d_udata,d_f,d_c);
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
	CHECK(cudaFree(d_udata));
	CHECK(cudaFree(d_c));

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
