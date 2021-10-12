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
// grid 2D block 2D
__global__ void matmulOnGPU(float *A, float *B, float *C, const int N)
{
    unsigned int ix = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
    unsigned int idx;
    unsigned int ida, idb;

    float cc = 0.0f;

    if (ix < N && iy < N)
    {
        for (int k = 0; k < N; k++)
        {
            ida = iy*N + k;
            idb = k *N + ix;
            cc += A[ida]*B[idb];
        }
    }
    idx = iy * N + ix;
    C[idx] = cc;
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

int main(int argc, char* argv[])
{

    // set up data size of matrix
    int N = 1 << 11;
   
    int nxy = N * N;
    int nBytes = nxy * sizeof(float);
    //printf("Matrix size: nx %d ny %d\n", N, N);

    // malloc host memory
    float *h_A, *h_B, *hostRef, *gpuRef;
    h_A = (float *)malloc(nBytes);
    h_B = (float *)malloc(nBytes);
    hostRef = (float *)malloc(nBytes);
    gpuRef = (float *)malloc(nBytes);


    // initialize data at host side
    initialData(h_A,2.0f,nxy);
    initialData(h_B,0.5f, nxy);

    memset(hostRef, 0, nBytes);
    memset(gpuRef, 0, nBytes);

    double iStart = seconds();
    matmulOnHost(h_A, h_B, hostRef, N);
    double iElaps_h = seconds() - iStart;
    //printf("matmul elapsed %f sec\n", iElaps);

    // malloc device global memory
    float *d_A, *d_B, *d_C;
    CHECK(cudaMalloc((void **)&d_A, nBytes));
    CHECK(cudaMalloc((void **)&d_B, nBytes));
    CHECK(cudaMalloc((void **)&d_C, nBytes));

    // transfer data from host to device
    CHECK(cudaMemcpy(d_A, h_A, nBytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_B, h_B, nBytes, cudaMemcpyHostToDevice));
    //CHECK(cudaMemcpy(d_C, gpuRef, nBytes, cudaMemcpyHostToDevice));

    // invoke kernel at host side
    int nb = atol(argv[1]); //block size

    int dimx = nb;
    int dimy = nb;
    dim3 block(dimx, dimy);
    dim3 grid((N + block.x - 1) / block.x, (N + block.y - 1) / block.y);

    iStart = seconds();
    matmulOnGPU<<<grid, block>>>(d_A, d_B, d_C, N);
    CHECK(cudaDeviceSynchronize());
    double iElaps_g = seconds() - iStart;
    //printf("matmulOnGPU <<<(%d,%d), (%d,%d)>>> elapsed %f sec\n", grid.x,
          // grid.y,
           //block.x, block.y, iElaps);
    printf("%d,%d,%f,%f\n",N,nb, iElaps_h,iElaps_g);
    // check kernel error
    CHECK(cudaGetLastError());

    // copy kernel result back to host side
    CHECK(cudaMemcpy(gpuRef, d_C, nBytes, cudaMemcpyDeviceToHost));

    // check device results
    //checkResult(hostRef, gpuRef, nxy);

    // free device global memory
    CHECK(cudaFree(d_A));
    CHECK(cudaFree(d_B));
    CHECK(cudaFree(d_C));

    // free host memory
    free(h_A);
    free(h_B);
    free(hostRef);
    free(gpuRef);

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
