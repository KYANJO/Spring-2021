#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>


/*
 * This example demonstrates a matrix-vector multiplication on the CPU.
 */

 void initialData(float *ip, const float ival, int size)
 {
     for (int i = 0; i < size; i++)
     {
         ip[i] = (float)(rand() & 0xFF) / 100.0f;
     }
 
     return;
 }
 
 void matvecmulOnHost(float *A, float *x, float *b, const int N)
 {
     for (int i = 0; i < N; i++)
     {
         b[i] = 0;
     }
     
     for (int ix = 0; ix < N; ix++)
     {
        for (int iy = 0; iy < N; iy++)
        {   
            int id = iy*N+ix;
            b[ix] += A[id]*x[iy];
        }
     }        
}

__global__ void matvecmulOnGPU(float *A, float *x, float *b, const int N)
{
    unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
    float cc = 0.0f;

    if(id>N) return;
        for (int k = 0; k < N; k++)
        {   
            int idx = k*N + id;
            cc += A[idx]*x[k];
        }
        b[id] = cc;
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
 
 int main(int argc, char **argv)
 {
 
     // set up data size of matrix
     int N = 1 << 5;
    
     int nxa = N * N;
     int nx  = N;
     int nBytes = nxa * sizeof(float);
     int nBytx = nx * sizeof(float);
     printf("Matrix size: nx %d ny %d\n", N, N);
 
     // malloc host memory
    float *h_A, *h_x, *hostRef, *gpuRef;
    h_A = (float *)malloc(nBytes);
    h_x = (float *)malloc(nBytx);
    hostRef = (float *)malloc(nBytx);
    gpuRef = (float *)malloc(nBytx);

 
     // initialize data at host side
    initialData(h_A,2.0f,nxa);
    initialData(h_x,0.5f, nx);

    memset(hostRef, 0, nBytx);
    memset(gpuRef, 0, nBytx);

    double iStart = seconds();
    matvecmulOnHost(h_A, h_x, hostRef, N);
    double iElaps_h = seconds() - iStart;
    printf("matvecmul elapsed %f sec\n", iElaps_h);

    // malloc device global memory
    float *d_A, *d_x, *d_b;
    CHECK(cudaMalloc((void **)&d_A, nBytes));
    CHECK(cudaMalloc((void **)&d_x, nBytx));
    CHECK(cudaMalloc((void **)&d_b, nBytx));

    // transfer data from host to device
    CHECK(cudaMemcpy(d_A, h_A, nBytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_x, h_x, nBytx, cudaMemcpyHostToDevice));
   
    // invoke kernel at host side
    int nb = atol(argv[1]); //block size
    int dimx = nb;
    dim3 block(dimx,1);
    dim3 grid((N + block.x - 1) / block.x,1);

    iStart = seconds();
    matvecmulOnGPU<<<grid, block>>>(d_A, d_x, d_b, N);
    CHECK(cudaDeviceSynchronize());
    double iElaps_g = seconds() - iStart;
    printf("matvecmulOnGPU <<<(%d,%d), (%d,%d)>>> elapsed %f sec\n", grid.x,
           grid.y,
           block.x, block.y, iElaps_g);
    
    //printf("%d,%d,%f,%f\n",N,nb, iElaps_h,iElaps_g);       

    // check kernel error
    CHECK(cudaGetLastError());

    // copy kernel result back to host side
    CHECK(cudaMemcpy(gpuRef, d_b, nBytx, cudaMemcpyDeviceToHost));

    // check device results
    checkResult(hostRef, gpuRef, nx);

    // free device global memory
    CHECK(cudaFree(d_A));
    CHECK(cudaFree(d_x));
    CHECK(cudaFree(d_b));

    // free host memory
    free(h_A);
    free(h_x);
    free(hostRef);
    free(gpuRef);

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
