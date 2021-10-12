#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#define THREADS_PER_BLOCK 256
/*
 * This example demonstrates h_Amatrix-matrix multiplication on the CPU.
 */

void initialData(float *ip, const float ival, int size)
{
    for (int i = 0; i < size; i++)
    {
        ip[i] = (float)(rand() & 0xFF) / 100.0f;
    }

    return;
}

 void powerIterationOnHost(float *xint, float *A, float *x, float *V, const int N)
 {
    int n = 100;
    x[0] = *xint;
    
    for(int j=0; j<n; j++)
    {
        //matrix vect product
        int nx  = N;
        int nByt = nx * sizeof(float);
        float *b;
        b = (float *)malloc(nByt);

        for (int ix = 0; ix < N; ix++)
        { 
            float temp = 0.0f;
            for (int iy = 0; iy < N; iy++)
            {   
                int id = iy*N+ix;
                temp += A[id]*x[iy];
            }
            b[ix] = temp;
        }
        
        //euclidean norm
        float cc = 0.0f;

       for (int i = 0; i < N; i++)
       {   
           cc += b[i]*b[i];
       }
        
       float normb = pow(cc,0.5);   
       printf("norm = %f\n", normb);
        //normalize
        for(int i = 0; i<N; i++)
        {
            V[i+1] = b[i]/(normb);
        }
        //printf("V = %f", V[50]);
        //check for convergence
        float epsilon = 1e-6;
        if (fabs(V[j+1] - V[j]) < epsilon)
        {
            break;
        }

    }
    
 }

__global__ void powerIterationOnGPU(float *xint, float *A, float *x, float *V, const int N)
 {
    int n = 100;
    x[0] = *xint;
    int nx  = N;
    int nByt = nx * sizeof(float);
    float *normb;
    float *b;
    b = (float *)malloc(nByt);
    normb = (float *)malloc(sizeof(float));

    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    for(int j=0; j<n; j++)
    //if(i<n)
    {
        //matrix vect product
        float cc = 0.0f;

        if(i<N)
        {
            for (int k = 0; k < N; k++)
            {   
                int idx = k*N + i;
                cc += A[idx]*x[k];
            }
            b[i] = cc;
        }
        //euclidean norm
        int temp = 0.0f;
    
        __syncthreads();
        if(i<N)
        {
            temp = b[i]*b[i];
            atomicAdd(normb,temp); //add the result to the global sum
        }
        __syncthreads();
            *normb = pow(*normb,0.5);
            //printf("norm = %f\n", *normb);
        //normalize
        if(i<N)
        {
            V[i+1] = b[i]/(*normb);
        }   
        
        //printf("V = %f\n", V[50]);

        //check for convergence
        float epsilon = 1e-6;
        if (fabs(V[j+1] - V[j]) < epsilon)
        {
            break;
        }

    }
    
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

int main(int argc, char **argv)
{

    // set up dath_Asize of matrix
    int N = 1 << 5;
   
    int nxy = N * N;
    int nx  = N;
    int nByt = nx * sizeof(float);
    int nBytes = nxy * sizeof(float);
    //printf("Matrix size: nx %d ny %d\n", N, N);

    // malloc host memory
    float *h_A, *h_x, *h_xint, *hostRef, *gpuRef;
    h_A= (float *)malloc(nBytes);
    h_x = (float *)malloc(nByt);
    h_xint = (float *)malloc(nByt);
    hostRef = (float *)malloc(nByt);
    gpuRef = (float *)malloc(nByt);

    double iStart = seconds();
    // initialize dath_Aat host side
    initialData(h_A,  2.0f, nxy);
    initialData(h_x,  0.5f, nx);
    initialData(h_xint,  0.5f, nx);

    memset(hostRef, 0, nByt);
    memset(gpuRef, 0, nByt);
 
    powerIterationOnHost(h_xint, h_A, h_x, hostRef, N);
    double iElaps_h = seconds() - iStart;
    //Elapsed time
    printf("N = %d\t Elapsed time = %f\n",N,iElaps_h);

    // malloc device global memory
    float *d_A;
    float *d_x;
    float *d_xint;
    float *d_V;
    CHECK(cudaMalloc((void **)&d_A, nBytes));
    CHECK(cudaMalloc((void **)&d_x, nByt));
    CHECK(cudaMalloc((void **)&d_xint, nByt));
    CHECK(cudaMalloc((void **)&d_V, nByt));

    // transfer data from host to device
    CHECK(cudaMemcpy(d_A, h_A, nBytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_x, h_x, nByt, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_xint, h_xint, nByt, cudaMemcpyHostToDevice));

    // invoke kernel at host side
    int nb = atol(argv[1]); //block size
    int dimx = nb;
    dim3 block(dimx,1);
    dim3 grid((N + block.x - 1) / block.x,1);

    iStart = seconds();
    powerIterationOnGPU<<<grid, block>>>(d_xint, d_A, d_x, d_V, N);
    CHECK(cudaDeviceSynchronize());
    double iElaps_g = seconds() - iStart;
    printf("powerIterationOnGPU <<<(%d,%d), (%d,%d)>>> elapsed %f sec\n", grid.x,
           grid.y,
           block.x, block.y, iElaps_g);

     // check kernel error
     CHECK(cudaGetLastError());

     // copy kernel result back to host side
    CHECK(cudaMemcpy(gpuRef, d_V, nByt, cudaMemcpyDeviceToHost));

    // check device results
    checkResult(hostRef, gpuRef, nx);

    // free device global memory
    CHECK(cudaFree(d_A));
    CHECK(cudaFree(d_x));
    CHECK(cudaFree(d_xint));
    CHECK(cudaFree(d_V));

    free(h_A);
    free(h_x);
    free(h_xint);
    free(hostRef);
    free(gpuRef);

     // reset device
     CHECK(cudaDeviceReset());

    return (0);
}
