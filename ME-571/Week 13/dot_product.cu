#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>

#define THREADS_PER_BLOCK 256

void checkResult(int hostRef, int gpuRef)
{
  //double epsilon = 1.0E-7;

    if (hostRef!=gpuRef)//abs(hostRef - gpuRef)/abs(hostRef) > epsilon)
      {
	printf("Numbers do not match!\n");
	printf("host %d gpu %d\n", hostRef, gpuRef);
      }else{
      printf("Results match.\n\n");
    }

    return;
}


void initialData(int *ip, int size)
{
    // generate different seed for random number
    time_t t;
    srand((unsigned) time(&t));

    for (int i = 0; i < size; i++)
    {
      ip[i] = (int)( rand() & 0xFF ) / 10.0f;
    }

    return;
}

void dotProductOnHost(int *A, int *B, int *C, const int N)
{
    for (int idx = 0; idx < N; idx++)
    {
        *C += A[idx] * B[idx];
    }
}

__global__ void dotProductOnGPU_globalAtomic(int *A, int *B, int *C, const int N)
{
  int temp = 0;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i<N) {
      temp = A[i] * B[i]; //each thread computes one entry of temp (no data race)
      atomicAdd(C,temp); //add the result to global sum
    }
    
}


__global__ void dotProductOnGPU(int *A, int *B, int *C, const int N)
{
    __shared__ int temp[THREADS_PER_BLOCK]; //define a shared array with entries for all threads in a block
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i<N) {
      temp[threadIdx.x] = A[i] * B[i]; //each thread computes one entry of temp (no data race)
    }

    __syncthreads();

    if(threadIdx.x ==0) { //thread 0 sums up all the entries (can be done much better with some fancy reduction)
      int sum = 0;
      for (int j=0; j<THREADS_PER_BLOCK; j++)
	sum += temp[j];
      atomicAdd(C,sum); //add the resulting (block local) sum to global C
    }

    
}

/*__global__ void dotProductOnGPU_sharedAtomic(int *A, int *B, int *C, const int N)
{
    __shared__ int temp[THREADS_PER_BLOCK]; //define a shared array with entries for all threads in a block
    __shared__ int c_block;
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i<N) {
      temp[threadIdx.x] = A[i] * B[i]; //each thread computes one entry of temp (no data race)
    }

    atomicAdd_block(c_block,temp[threadIdx.x]);
    __syncthreads();

    
    if(threadIdx.x ==0) { //thread 0 sums up all the entries (can be done much better with some fancy reduction)
      atomicAdd(C,c_block); //add the resulting (block local) sum to global C
    }

    
    }*/


int main(int argc, char **argv)
{
  printf("%s Starting...\n", argv[0]);
  
  // set up device
  int dev = 0;
  cudaDeviceProp deviceProp;
  CHECK(cudaGetDeviceProperties(&deviceProp, dev));
  printf("Using Device %d: %s\n", dev, deviceProp.name);
  CHECK(cudaSetDevice(dev));
  
  // set up data size of vectors
  int nElem = 1 << 25;
  printf("Vector size %d\n", nElem);
  
  // malloc host memory
  size_t nBytes = nElem * sizeof(int);
  
  int *h_A, *h_B;
  h_A     = (int *)malloc(nBytes);
  h_B     = (int *)malloc(nBytes);
  int host_result;
  int gpu_result;
  
  double iStart, iElaps;
  
  // initialize data at host side
  iStart = seconds();
  initialData(h_A, nElem);
  initialData(h_B, nElem);
  iElaps = seconds() - iStart;
  printf("initialData Time elapsed %f sec\n", iElaps);
  host_result = 0;
  gpu_result = 0;

  // add vector at host side for result checks
  iStart = seconds();
  dotProductOnHost(h_A, h_B, &host_result, nElem);
  iElaps = seconds() - iStart;
  printf("dotProductOnHost Time elapsed %f sec\n", iElaps);

  // malloc device global memory
  int *d_A, *d_B, *d_C;
  CHECK(cudaMalloc((int **)&d_A, nBytes));
  CHECK(cudaMalloc((int **)&d_B, nBytes));
  CHECK(cudaMalloc((int **)&d_C, sizeof(int)));
  
  // transfer data from host to device
  CHECK(cudaMemcpy(d_A, h_A, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(d_B, h_B, nBytes, cudaMemcpyHostToDevice));
  CHECK(cudaMemset(d_C, 0, sizeof(int))); //set d_C to zero on the device

  // invoke kernel at host side
  int nThreads = THREADS_PER_BLOCK;
  int nBlocks =  ((nElem + nThreads - 1) / nThreads);
  
  iStart = seconds();
  //  dotProductOnGPU_globalAtomic<<<nBlocks, nThreads>>>(d_A, d_B, d_C, nElem);
  dotProductOnGPU<<<nBlocks, nThreads>>>(d_A, d_B, d_C, nElem);
  CHECK(cudaDeviceSynchronize());
  iElaps = seconds() - iStart;
  printf("dotProductOnGPU <<<  %d, %d  >>>  Time elapsed %f sec\n", nBlocks,
	 nThreads, iElaps);

  
  // check kernel error
  CHECK(cudaGetLastError()) ;
  
  // copy kernel result back to host side
  CHECK(cudaMemcpy(&gpu_result, d_C, sizeof(int), cudaMemcpyDeviceToHost));
  
  // check device results
  checkResult(host_result, gpu_result);
  
  // free device global memory
  CHECK(cudaFree(d_A));
  CHECK(cudaFree(d_B));
  CHECK(cudaFree(d_C));
  
  // free host memory
  free(h_A);
  free(h_B);
  
  return(0);
}
