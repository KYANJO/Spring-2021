#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#define THREADS_PER_BLOCK 256
/*
 * This example demonstrates a matrix-vector multiplication on the CPU.
 */

 void initialData(float *ip, int size)
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
 
//Euclidean on Host
 float euclidean(float *b, const int N)
 {
     float cc = 0.0f;

        for (int i = 0; i < N; i++)
        {   
            cc += b[i]*b[i];
        }
         
        float normb = pow(cc,0.5);   
        return normb;
}

void normOnHost(float *b, float *normb, float *V, const int N)
{
      for(int i = 0; i<N; i++)
      {
        V[i] = b[i]/(*normb);
      }
}


//Normalisation kernel
__global__ void normOnGPU(float *b, float normb, float *V, const int N)
{  
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  
      if(i<N)
      {
        V[i] = b[i]/normb;
      }     
}

void checkResult(float*hostRef, float*gpuRef, const int N)
{
    double epsilon = 1.0E-6;
    bool match = 1;
    for (int i = 0; i < N; i++)
    {
      if (abs(hostRef[i] - gpuRef[i])/abs(hostRef[i]) > epsilon)
        {
            match = 0;
            printf("host %f gpu %f, err = %e\n", hostRef[i], gpuRef[i], abs(hostRef[i]-gpuRef[i])/abs(hostRef[i]));
            break;
        }
    }
    if (match)
        printf("Arrays match.\n\n");
    else
        printf("Arrays do not match.\n\n");
}
 
 int main(int argc, char **argv)
 {
 
     // set up data size of matrix
     int N = 1 << 5;
    
     int nx  = N;
     size_t nBytes = nx * sizeof(int);
     printf("Vector size %d\n", N);
 
     // malloc host memory
    float *h_b;
    float *hostRef, *gpuRef;
    h_b = (float *)malloc(nBytes);
    hostRef = (float *)malloc(nBytes);
    gpuRef = (float *)malloc(nBytes);

     // initialize data at host side
    initialData(h_b,nx);

    float normb = euclidean(h_b, N);

    memset(hostRef, 0, nBytes);
    memset(gpuRef, 0, nBytes);

    double iStart = seconds();
    normOnHost(h_b, &normb, hostRef, N);
    double iElaps_h = seconds() - iStart;
    printf("normOnHost elapsed %f sec\n", iElaps_h);

    // malloc device global memory
    float *d_b;
    float *d_norm;
    float d_normb = normb;
    CHECK(cudaMalloc((void **)&d_b, nBytes));
    CHECK(cudaMalloc((void **)&d_norm, nBytes));
    
    // transfer data from host to device
    CHECK(cudaMemcpy(d_b, h_b, nBytes, cudaMemcpyHostToDevice));
    
    // invoke kernel at host side
    int nThreads = THREADS_PER_BLOCK;
    int nBlocks =  ((N+ nThreads - 1) / nThreads);

    iStart = seconds();
    normOnGPU<<<nBlocks, nThreads>>>(d_b, d_normb, d_norm, N);
    CHECK(cudaDeviceSynchronize());
    double iElaps_g = seconds() - iStart;
    printf("normOnGPU <<<(%d,%d), (%d,%d)>>> Time elapsed %f sec\n", nBlocks,
    nThreads, iElaps_g);

    //printf("%d,%d,%f,%f\n",N,nb, iElaps_h,iElaps_g);       

    // check kernel error
    CHECK(cudaGetLastError());

    // copy kernel result back to host side
    CHECK(cudaMemcpy(gpuRef, d_norm, nBytes, cudaMemcpyDeviceToHost));

    // check device results
    checkResult(hostRef, gpuRef, nx);

    // free device global memory
    CHECK(cudaFree(d_b));
    CHECK(cudaFree(d_norm));
    

    // free host memory
    free(h_b);
    free(hostRef);
    free(gpuRef);

    // reset device
    //CHECK(cudaDeviceReset());

    return (0);
}
