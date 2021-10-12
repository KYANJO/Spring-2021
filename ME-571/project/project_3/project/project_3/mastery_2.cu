#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#define THREADS_PER_BLOCK 256
/*
 * This example demonstrates a matrix-vector multiplication on the CPU.
 */

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
 
 void euclideanOnHost(int *b, int *normb, const int N)
 {
     int cc = 0.0f;

        for (int i = 0; i < N; i++)
        {   
            cc += b[i]*b[i];
        }
        *normb = pow(cc,0.5);       
}

__global__ void euclideanOnGPU(int *b, int *normb, const int N)
{
    
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    //int *c;
    int temp = 0.0f;
    
    __syncthreads();
    if(i<N)
    {
        temp = b[i]*b[i];
        atomicAdd(normb,temp); //add the result to the global sum
    }
    __syncthreads();
        *normb = pow(*normb,0.5);
    
}

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
 
 int main(int argc, char **argv)
 {
 
     // set up data size of matrix
     int N = 1 << 5;
    
     int nx  = N;
     size_t nBytes = nx * sizeof(int);
     printf("Vector size %d\n", N);
 
     // malloc host memory
    int *h_b, hostRef, gpuRef;
    h_b = (int *)malloc(nBytes);

     // initialize data at host side
    initialData(h_b,nx);

    hostRef = 0;
    gpuRef = 0;

    double iStart = seconds();
    euclideanOnHost(h_b, &hostRef, N);
    double iElaps_h = seconds() - iStart;
    printf("euclideanOnHost elapsed %f sec\n", iElaps_h);

    // malloc device global memory
    int *d_b, *d_norm;
    CHECK(cudaMalloc((void **)&d_b, nBytes));
    CHECK(cudaMalloc((void **)&d_norm, sizeof(int)));

    // transfer data from host to device
    CHECK(cudaMemcpy(d_b, h_b, nBytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemset(d_norm, 0, sizeof(int)));

    // invoke kernel at host side
    int nThreads = THREADS_PER_BLOCK;
    int nBlocks =  ((N+ nThreads - 1) / nThreads);

    iStart = seconds();
    euclideanOnGPU<<<nBlocks, nThreads>>>(d_b,d_norm, N);
    CHECK(cudaDeviceSynchronize());
    double iElaps_g = seconds() - iStart;
    printf("euclideanOnGPU <<<(%d,%d), (%d,%d)>>> Time elapsed %f sec\n", nBlocks,
    nThreads, iElaps_g);

    //printf("%d,%d,%f,%f\n",N,nb, iElaps_h,iElaps_g);       

    // check kernel error
    CHECK(cudaGetLastError());

    // copy kernel result back to host side
    CHECK(cudaMemcpy(&gpuRef, d_norm, sizeof(int), cudaMemcpyDeviceToHost));

    // check device results
    checkResult(hostRef, gpuRef);

    // free device global memory
    CHECK(cudaFree(d_b));
    

    // free host memory
    free(h_b);

    // reset device
    //CHECK(cudaDeviceReset());

    return (0);
}
