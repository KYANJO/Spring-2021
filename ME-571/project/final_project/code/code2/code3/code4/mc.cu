/******************************************************************************
 * Implements the Monte-Carlo method of performing numerical integral for complex
 * integrals on the GPU
 * 
 * Author: Yao GAhounzo
 * Date : 04/21/2021
 */


 #include <curand.h>
 #include <curand_kernel.h>
 #include "common.h"
 #include <cuda_runtime.h>
 #include <stdio.h>
 #include <math.h>
 
 #define THREADS_PER_BLOCK 64

 // performs the montecarlo integration on the GPU
 
 __global__ void monteCarloOnGPU(curandState *state, float *fx, float *Integ, const int N, const int BT){	
     
    __shared__ float partialSum[THREADS_PER_BLOCK];
    //__shared__ float sumThread;
    float sumThread = 0;
    double x;
    
    unsigned int id = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int id_loc = threadIdx.x;

    int Nx = N/BT;
     
    if(id < Nx){	

		// Initialize CURAND
		curand_init(id, id, 0, &state[id]);

		for(int i = 0; i < BT; i++){

			x  = curand_uniform(&state[id]);
			sumThread += cos(-log(x));
			 
		}

    }      
    
    //__syncthreads();
    partialSum[id_loc] = sumThread;
    __syncthreads();
    
    // Reduction in shared memory
    for (unsigned int i = 1; i < blockDim.x; i *= 2){
        if (id_loc % (2*i) == 0) {
            partialSum[id_loc] += partialSum[id_loc + i];
        }
        __syncthreads();
    }
    
    //__syncthreads();
    if (id_loc == 0) fx[blockIdx.x] = partialSum[0];


 }
 
 
 __global__ void Integral(float *fx, float *Integ, const int N){

     int id = blockIdx.x;
    __syncthreads();
 
    atomicAdd(Integ, fx[id]);
     
   // __syncthreads();
 
    *Integ = *Integ / N; 

 }
 
 int main(int argc, char **argv){
     
     printf("\n");
     printf("%s Starting...\n", argv[0]);
 
     // set up device
     int dev = 0;
     cudaDeviceProp deviceProp;
     CHECK(cudaGetDeviceProperties(&deviceProp, dev));
     printf("Using Device %d: %s\n\n", dev, deviceProp.name);
     CHECK(cudaSetDevice(dev));
 
 
     int N = atoi(argv[1]);
     
     double iStart, iElaps, iStart_comp, iElaps_comp;
     iStart = seconds();
     // host memory
     //int size = sizeof(float);
     float gpuRef,  error;
     float exact = 0.5;
 
     // malloc device memory
     float *Integ_d, *fx_d;
     CHECK(cudaMalloc((float **)&Integ_d, sizeof(float)));
     //CHECK(cudaMalloc((float **)&fx_d, sizeof(float)));
 
     // invoke kernel
     //int nThreads = 64;
     int nThreads = THREADS_PER_BLOCK;
     int nBlocks = ((N + nThreads - 1)/nThreads);
     //int nBlocks = N/nThreads;
 
     if(nBlocks > 65535){
         nBlocks = 65535;
     }
     
     int BT = N/(nThreads*nBlocks);
     if(BT == 0) BT = 1;
     int Nx = N/BT;
 
     curandState *devStates;
     //float *random_d;
     CHECK(cudaMalloc((void **)&devStates, Nx*sizeof(curandState)));
     //CHECK(cudaMalloc((void **)&random_d, Nx*sizeof(float)));
     CHECK(cudaMalloc((float **)&fx_d, N*sizeof(float)));
 
     // On GPU
     iStart_comp = seconds();
     monteCarloOnGPU<<<nBlocks, nThreads>>>(devStates, fx_d, Integ_d, N, BT);
     CHECK(cudaDeviceSynchronize());
     iElaps_comp = seconds() - iStart_comp;
     
     iStart_comp = seconds();
     Integral<<<2, 1024>>>(fx_d, Integ_d, N);
     CHECK(cudaDeviceSynchronize());
     iElaps_comp = seconds() - iStart_comp;
     
     // check kernel error
     CHECK(cudaGetLastError());
 
     // copy kernel result back to host side
     CHECK(cudaMemcpy(&gpuRef, Integ_d, sizeof(float), cudaMemcpyDeviceToHost));
 
     // compute error
     error = abs(exact - gpuRef);
 
     iElaps = seconds() - iStart;
 
     printf("N = %ld, Integral = %f, error = %e, elapsed_time = %f s, time_comp = %f s\n\n",N,gpuRef,error,iElaps,iElaps_comp);
 
     // free device global memory
     CHECK(cudaFree(Integ_d));
     CHECK(cudaFree(fx_d));
 
     // reset device
     CHECK(cudaDeviceReset());
     
     return 0;
 }
 