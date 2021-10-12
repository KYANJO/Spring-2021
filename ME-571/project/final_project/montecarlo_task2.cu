#include <cuda_runtime.h>
#include <curand.h>
#include "common.h"
#include <stdio.h>
#include <curand_kernel.h>

#define THREADS_PER_BLOCK 256

__global__ void mcOnGPU(float *f, const int N, const int nb, curandState *states)
{   
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
	int n = N/nb;
	float cc = 0.0f;
	//initialse curand
	curand_init(1234, ix, 0, &states[ix]);
	//curandState state = states[ix];
	if (ix<n)
	{
		float xran = curand_uniform (&states[ix]);
		cc += cos(-log(xran));
		
		//__syncthreads();

		atomicAdd(f,cc); 
		__syncthreads();
	}
					
}

__global__ void integralOnGPU(float *f,float *Int ,const int N)
{
	Int[0] = abs(f[0]/N);
}


int main(int argc, char **argv)
{	
	// problem size
	int N = atoi(argv[1]);
	int T = THREADS_PER_BLOCK;
	
	//random number generator
	curandState *States;
	
	// malloc host memory
    float gpuRef;
	
	// malloc device global memory
    float *d_f;
    float *d_Int;
    CHECK(cudaMalloc((void **)&d_f, sizeof(float)));
    CHECK(cudaMalloc((void **)&d_Int, sizeof(float)));
    
	//invoke the kernel
	int B = ((N + T -1)/T);
	if(B > 65535) B = 65535;
	int nb = N/(B*T);
	
	//states allocate memory
	CHECK(cudaMalloc( (void **)&States, (B*T)*sizeof(curandState)));

	mcOnGPU<<<B,T>>>(d_f, N, nb, States);
	CHECK(cudaDeviceSynchronize());

	integralOnGPU<<<1,1>>>(d_f,d_Int ,N);
	CHECK(cudaDeviceSynchronize());

	// copy kernel result back to host side
    CHECK(cudaMemcpy(&gpuRef, d_Int, sizeof(float), cudaMemcpyDeviceToHost));

	//error achived
	float error = 0.5 - gpuRef;
	printf("%d,%f,%e\n",N,gpuRef,error);
	
	//free device memory
	CHECK(cudaFree(States));
	CHECK(cudaFree(d_f));
    CHECK(cudaFree(d_Int));

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
