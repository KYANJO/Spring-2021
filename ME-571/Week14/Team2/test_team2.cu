#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#include "kernel.cu"

__global__ void test_kernels(float *u, float *u_new, float *f, int N)
{
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int j = tjreadIdx.y + blockIdx.y * blockDIm.y;

  unsigned int idx = j * N + i;

  if(i<N && j<N) {
    float h = 1/(N-1);
    float x = i*h;
    float y = j*h;
    float ue = 0;
    float fe = -8*M_PI*M_PI*sin(2*pi*x)*sin(2*pi*y);
    float ue_new = -0.25*h*h*fe;
    printf("Thread (%d, %d), Block (%d,%d), u = %f, u_new = %f, f = %f, expected values: %f, %f, %f\n",threadIdx.x,threadIdx.u,blockIdx.x,blockIdx.y,u[idx],u_new[idx],f[idx],ue,ue_new,fe);
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

    // set up problem size
    int N = 6
    int nxy = N * N;
    int nBytes = nxy * sizeof(float);
    printf("Problem size: nx %d ny %d\n", N, N);
    
    // malloc device global memory
    float *d_u, *d_f, *d_u_new;
    CHECK(cudaMalloc((void **)&d_u     , nBytes));
    CHECK(cudaMalloc((void **)&d_u_new , nBytes));
    CHECK(cudaMalloc((void **)&d_f     , nBytes));


    //set-up blocks and threads
    int dimx = 3;
    int dimy = 3;
    dim3 block(dimx, dimy);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
    
    //initialize arrays
    initData<<<grid,block>>>(d_u, d_f, N)
    CHECK(cudaDeviceSynchronize());
    //take one Poisson update 
      poissonKernel<<<grid, block>>>(d_u, d_u_new, d_f, N);
    CHECK(cudaDeviceSynchronize());

    // check kernel error
    CHECK(cudaGetLastError());

    test_kernels<<<grid, block>>>(d_u, d_u_new, d_f, int N)
    // free device global memory
    CHECK(cudaFree(d_u));
    CHECK(cudaFree(d_f));
    CHECK(cudaFree(d_u_new));

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
