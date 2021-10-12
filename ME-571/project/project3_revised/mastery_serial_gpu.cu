#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#define THREADS_PER_BLOCK 1024

 void initialData(float *ip, const float ival, int size)
 {
     for (int i = 0; i < size; i++)
     {
         ip[i] = (float)(rand() & 0xFF) / 100.0f;
     }
 
     return;
 }
 

//matrix vector multiplication on Host
void matvecmulOnHost(float *A, float *x, float *b, const int N)
{
    for (int i = 0; i < N; i++) b[i] = 0;
    
    for (int ix = 0; ix < N; ix++)
    {
       for (int iy = 0; iy < N; iy++)
       {   
           int id = iy*N+ix;
           b[ix] += A[id]*x[iy];
       }
    }        
}

//Euclidean norm on Host
void euclideanOnHost(float *b, float *normb, const int N)
 {
     float cc = 0.0f;

        for (int i = 0; i < N; i++)  cc += b[i]*b[i];
    
        *normb = pow(cc,0.5);       
}

//Normalisation on Host
void normOnHost(float *b, float *normb, float *V, const int N)
{
      for(int i = 0; i<N; i++) V[i] = b[i]/(*normb);
    
}


//matrix vector multiplication on kernel
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

//Euclidean norm on host
__global__ void euclideanOnGPU(float *b, float *cc, float *normb, const int N)
{
    float temp = 0;
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    cc[i] = b[i]*b[i];

    for(int i = 0; i< N; i++) temp += cc[i];
    
    *normb  = pow(temp,0.5);
    
}

//Normalisation kernel
__global__ void normOnGPU(float *b, float *normb, float *V, const int N)
{  
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  
      if(i<N) V[i] = b[i]/(*normb);
        
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
            printf("host %f gpu %f\n", hostRef[i], gpuRef[i]);
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
    int N = 1 << 11;
    
    int nxa = N * N;
    int nx  = N;
    int nBytes = nxa * sizeof(float);
    int nBytx = nx * sizeof(float);
    //printf("Matrix size: nx %d ny %d\n", N, N);

    double iStart, iElaps_g, iend_h;
    iStart = seconds();

    // malloc host memory
   float *h_A, *h_x, *h_b, *hostRef, *gpuRef;
   h_A = (float *)malloc(nBytes);
   h_x = (float *)malloc(nBytx);
   h_b = (float *)malloc(nBytx);
   //h_normb = (float *)malloc(sizeof(float));
   hostRef = (float *)malloc(nBytx);
   gpuRef = (float *)malloc(nBytx);

   memset(hostRef, 0, nBytx);
   memset(gpuRef, 0, nBytx);
   memset(h_b, 0, nBytx);
   float h_normb = 0;

   // initialize data at host side
   initialData(h_A,2.0f,nxa);
   initialData(h_x,0.5f, nx);
   initialData(hostRef,2.0f, nx);

   //On host
   // number of simulations
   int n = 100;

   for (int i = 0; i < n; i++)
   {
    //Matrix vector multiplication on Host
    matvecmulOnHost(h_A, hostRef, h_b, N);

    //Euclidean norm on host
    euclideanOnHost(h_b, &h_normb, N);

    //Normalisation on Host
    normOnHost(h_b, &h_normb, hostRef, N);
   }
   iend_h = seconds() - iStart; 

   //start timing gpu
   iStart = seconds();
   // malloc device global memory
   float *d_A, *d_xo, *d_x, *d_b, *d_normb, *d_c;
   CHECK(cudaMalloc((void **)&d_A, nBytes));
   CHECK(cudaMalloc((void **)&d_xo, nBytx));
   CHECK(cudaMalloc((void **)&d_x, nBytx));
   CHECK(cudaMalloc((void **)&d_b, nBytx));
   CHECK(cudaMalloc((void **)&d_c, nBytx));
   CHECK(cudaMalloc((void **)&d_normb, sizeof(float)));

   // transfer data from host to device
   CHECK(cudaMemcpy(d_A, h_A, nBytes, cudaMemcpyHostToDevice));
   CHECK(cudaMemcpy(d_x, h_x, nBytx, cudaMemcpyHostToDevice));
   CHECK(cudaMemcpy(d_xo, d_x, nBytx, cudaMemcpyDeviceToDevice));

   CHECK(cudaMemset(d_normb, 0, sizeof(float)));

    // invoke kernel at host side
    int nthreads = THREADS_PER_BLOCK;
    int nb = atol(argv[1]); //block size
    int dimx = nb;
    int dimy = nb;
    dim3 block(dimx,dimy);
    dim3 grid((N + block.x - 1) / block.x,(N + block.y - 1) / block.y);
    dim3 nblocks((N + nthreads - 1) / nthreads);

    for (int i = 0; i < n; i++)
   {
       //Matrix vector multiplication on GPU
    matvecmulOnGPU<<<grid, block>>>(d_A, d_x, d_b, N);
    CHECK(cudaDeviceSynchronize());

    //Euclidean norm on GPU
    euclideanOnGPU<<<nblocks, nthreads>>>(d_b, d_c, d_normb, N);
    CHECK(cudaDeviceSynchronize());

    //Normalisation on Host
    normOnGPU<<<nblocks, nthreads>>>(d_b, d_normb, d_x, N);
    CHECK(cudaDeviceSynchronize());
   }

   // check kernel error
   CHECK(cudaGetLastError());

    // copy kernel result back to host side
    CHECK(cudaMemcpy(gpuRef, d_x, nBytx, cudaMemcpyDeviceToHost));
    iElaps_g = seconds() - iStart;
    //Elapsed time
    printf("N = %d\t serial_time = %f\t kernel_time = %f\n",N,iend_h,iElaps_g);

    // check device results
    checkResult(hostRef, gpuRef, nx);

    // free device global memory
    CHECK(cudaFree(d_A));
    CHECK(cudaFree(d_x));
    CHECK(cudaFree(d_b));

   //free host memory
   free(h_A);
   free(h_b);
   free(hostRef);
   free(gpuRef);

   CHECK(cudaDeviceReset());

    return (0);

}


