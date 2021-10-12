#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#include "kernel.cu"

__global__ void test_init(float *u, float *u_new, int N)
{
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int j = tjreadIdx.y + blockIdx.y * blockDIm.y;

  unsigned int idx = j * N + i;

  if(i<N && j<N) {
    u[idx] = 0.0;
    u_new[idx] = 1.0;
  }
  
}

int main(int argc, char **argv)
{
    printf("%s Starting...\n", argv[0]);

    // set up device
    int dev = 0;
    cudaDeviceProp deviceProp;
    CHECK(cudaGetDeviceProperties(&deviceProp, dev));
    printf("Using Device %d: %s\n", dev, deviceProp.name);
    CHECK(cudaSetDevice(dev));

    float error

    // set up problem size
    int N = 6
    int nxy = N * N;
    int nBytes = nxy * sizeof(float);
    printf("Problem size: nx %d ny %d\n", N, N);
    
    // malloc device global memory
    float *d_u, *d_u_new, *d_error;
    CHECK(cudaMalloc((void **)&d_u     , nBytes));
    CHECK(cudaMalloc((void **)&d_u_new , nBytes));
    CHECK(cudaMalloc((void **)&d_error , sizeof(float)));
    
    //set-up blocks and threads
    int dimx = 3;
    int dimy = 3;
    dim3 block(dimx, dimy);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);

    //initialize u
    
    test_init<<<grid,block>>>(d_u, d_u_new, N)
    CHECK(cudaDeviceSynchronize());
    
    // compute error
    computeError<<<grid, block>>>(d_error, d_u, d_u_new, N);
    CHECK(cudaDeviceSynchronize());
    
    CHECK(cudaMemcpy(error, d_error, sizeof(float), cudaMemcpyDeviceToHost));
    printf("Error after initialization: %e, expected value: 1.0\n",error);
    
    // update

    updateSolution<<<grid,block>>>(d_u, d_u_new, N);
    CHECK(cudaDeviceSynchronize());
    
    //compute error again
    computeError<<<grid, block>>>(d_error, d_u, d_u_new, N);
    CHECK(cudaDeviceSynchronize());
    
    CHECK(cudaMemcpy(error, d_error, sizeof(float), cudaMemcpyDeviceToHost));
    printf("Error after update: %e, expected value: 0.0\n",error);
    
    // check kernel error
    CHECK(cudaGetLastError());

    // free device global memory
    CHECK(cudaFree(d_u));
    CHECK(cudaFree(d_u_new));
    CHECK(cudaFree(d_error));

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
