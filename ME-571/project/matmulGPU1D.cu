#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>

/*
 * This example demonstrates a matrix-matrix multiplication on the CPU.
 */

void initialData(float *ip, const float ival, int size)
{
    for (int i = 0; i < size; i++)
    {
        ip[i] = (float)(rand() & 0xFF) / 100.0f;
    }

    return;
}

void matmulOnHost(float *A, float *B, float *C, const int N)
{
  int   id, ida, idb;
  float cc;
  
    for (int iy = 0; iy < N; iy++)
    {
        for (int ix = 0; ix < N; ix++)
        {
    	    cc = 0;
    	    for (int k = 0; k < N; k++){
        	    ida = iy*N + k;
        	    idb = k *N + ix;
                cc += A[ida]*B[idb];
	        }
    	    id = iy*N+ix;
    	    C[id] = cc;
        }

    }

    return;
}

void printMatrix(float *C, const int nx, const int ny)
{
    float *ic = C;

    for (int iy = 0; iy < ny; iy++)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            printf("%f ", ic[ix]);

        }

        ic += nx;
        printf("\n");
    }

    return;
}

void checkResult(float *hostRef, float *gpuRef, const int N)
{
    double epsilon = 1.0E-8;
    bool match = 1;

    for (int i = 0; i < N; i++)
    {
        if (abs(hostRef[i] - gpuRef[i]) > epsilon)
        {
            match = 0;
            printf("host %f gpu %f\n", hostRef[i], gpuRef[i]);
            break;
        }
    }

    if (match)
        printf("Arrays match.\n\n");
    else
        printf("Arrays do not match.\n\n");
}

__global__ void matmulOnGPU1D(float *A, float *B, float *C, const int N)
{
  int   id;
  int ida, idb;
  float cc;
  
  unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
  
    if (ix < N)

        for (int iy = 0; iy < N; iy++){
            

            cc = 0;
            for (int k = 0; k < N; k++){
                ida = iy*N + k;
                idb = k *N + ix;
                cc += A[ida]*B[idb];
            }

            id = iy * N + ix;
            C[id] = cc;
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

    // set up data size of matrix
    int N = 1 << 11;
   
    int nxy = N * N;
    int nBytes = nxy * sizeof(float);
    printf("Matrix size: N %d\n",N);

    // malloc host memory
    float *h_A, *h_B, *hostRef, *gpuRef;
    h_A = (float *)malloc(nBytes);
    h_B = (float *)malloc(nBytes);
    hostRef = (float *)malloc(nBytes);
    gpuRef = (float *)malloc(nBytes);

    // initialize data at host side
    double iStart = seconds();
    initialData(h_A,  2.0f, nxy);
    initialData(h_B,  0.5f, nxy);
    double iElaps = seconds() - iStart;
    printf("Matrix initialization elapsed %f sec\n", iElaps);

    //memset(C, 0, nBytes);

    memset(hostRef, 0, nBytes);
    memset(gpuRef, 0, nBytes);

    iStart = seconds();
    matmulOnHost(h_A, h_B, hostRef, N);
    iElaps = seconds() - iStart;
    printf("matmulOnHost elapsed %f sec\n", iElaps);

    float *d_A, *d_B, *d_C;
    CHECK(cudaMalloc((void **)&d_A, nBytes));
    CHECK(cudaMalloc((void **)&d_B, nBytes));
    CHECK(cudaMalloc((void **)&d_C, nBytes));

    CHECK(cudaMemcpy(d_A, h_A, nBytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_B, h_B, nBytes, cudaMemcpyHostToDevice));

    int dimx = 1;
    int dimy = 1;
    dim3 block(dimx, dimy);
    dim3 grid((N + block.x - 1) / block.x, (N + block.y - 1) / block.y);

    iStart = seconds();
    matmulOnGPU1D<<<grid, block>>>(d_A, d_B, d_C, N);
    CHECK(cudaDeviceSynchronize());
    iElaps = seconds() - iStart;
    printf("matmulOnGPU1D <<<  (%d,%d), (%d,%d)  >>> elapsed %f sec\n",
           grid.x, grid.y, block.x, block.y, iElaps);


    // adjust block size
    block.x = 16;
    block.y = 1;
    grid.x  = (N + block.x - 1) / block.x;
    grid.y  = (N + block.y - 1) / block.y;

    iStart = seconds();
    matmulOnGPU1D<<<grid, block>>>(d_A, d_B, d_C, N);
    CHECK(cudaDeviceSynchronize());
    iElaps = seconds() - iStart;
    printf("matmulOnGPU1D <<<  (%d,%d), (%d,%d)  >>> elapsed %f sec\n",
           grid.x, grid.y, block.x, block.y, iElaps);

       
    CHECK(cudaMemcpy(gpuRef, d_C, nBytes, cudaMemcpyDeviceToHost));
    checkResult(hostRef, gpuRef, nxy);

    CHECK(cudaFree(d_A));
    CHECK(cudaFree(d_B));
    CHECK(cudaFree(d_C));

    free(h_A);
    free(h_B);
    free(hostRef);
    free(gpuRef);

    CHECK(cudaDeviceReset());

    return (0);
}
