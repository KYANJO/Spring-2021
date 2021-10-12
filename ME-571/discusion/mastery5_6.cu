#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>

/*
 * This example demonstrates a matrix-matrix multiplication on the CPU.
 */

 #define THREADS_PER_BLOCK 1024

 void initialData(float *ip, const float ival, const int size)
{
    for (int i = 0; i < size; i++)
    {
        ip[i] = (float)(rand() & 0xFF) / 100.0f;
    }

    return;
}

void matvecOnHost(float *A, float *x, float *b, const int N)
{
    
    for (int iy = 0; iy < N; iy++)
    {
        float M = 0;
        for (int ix = 0; ix < N; ix++)
        {
            
            M = M + A[iy*N + ix]*x[ix];
            
        }

        b[iy] = M;
    }

    return;
}

void EucldNormOnHost(float *b, float *res, const int N)
{
    float sum = 0.0;

    for (int k = 0; k < N; k++)
    {      
        sum = sum + b[k]*b[k];
    }

    *res = sqrt(sum);

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

void printVector(float *C, const int nx)
{
    float *ic = C;

    
    for (int ix = 0; ix < nx; ix++)
    {
        printf("%f ", ic[ix]);

    }
    printf("\n");
    
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
            printf("host %f gpu %f\n", hostRef[i], gpuRef[i]);
            break;
        }
    }

    if (match)
        printf("Arrays match.\n\n");
    else
        printf("Arrays do not match.\n\n");
}

void checkResultNorm(float hostRes, float gpuRes)
{
    //double epsilon = 1.0E-6;
    //(abs(hostRes - gpuRes)/abs(hostRes) > epsilon){//

    if (hostRes!=gpuRes){
        printf("Numbers do not match!\n");
        printf("host %d gpu %d\n", hostRes, gpuRes);
    }else{
      printf("Results match.\n\n");
    }

    return; 
}

void VecNormalizationOnHost(float *b, float *c, float *normb, const int N){

    for (int i = 0; i< N; i++){
        c[i] = b[i] / *normb;
    }
    
    return;
}

__global__ void matvecOnGPU(float *A, float *x, float *b, const int N)
{
    unsigned int iy = threadIdx.y + blockIdx.y * blockDim.y;
    
    if (iy< N){

        float bi = 0;
        
        for (int ix = 0; ix < N; ix++){

            bi= bi+ A[iy*N + ix]*x[ix];
        }
        b[iy] = bi;
    }
    //__syncthreads();     
}


__global__ void EucldNormOnGPU(float *b, float *res, float *c, const int N)
{

    float temp = 0;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    /*
    c[i] = b[i]*b[i];
    //_syncthreads();
    for(int k = 0; k< N; k++){

        temp += c[k];
    }
    *res  = sqrt(temp);
    */
    
    __syncthreads();

    if(i<N) {

      temp = b[i]*b[i]; //each thread computes one entry of temp (no data race)
      atomicAdd(c,temp); //add the result to global sum

    }
    
    __syncthreads();
    //if(i == 0)
      *res = sqrtf(*c);
      
}

__global__ void VecNormalizationOnGPU(float *b, float *c, float *normb, const int N){

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i<N){
        c[i] = b[i] / *normb;
    }
    
}

// Serial power iteration

void powerIterationOnHost(float *A, float *x0, float *xn, const int N){

    float resNorm;
    float *bn  = (float *)malloc(N);
    int  n;

    for (int i = 0; i < N; i++) {
        xn[i] = x0[i];
    }
    
    for (n = 0; n < 100; n++){

        // Matrix vector product
        
        matvecOnHost(A, xn, bn, N);

        // Euclidean norm
        
        EucldNormOnHost(bn, &resNorm, N);

        // Normalization

        VecNormalizationOnHost(bn, xn, &resNorm, N);

    return;
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
    int nnBytes = N*sizeof(float);
    printf("Matrix size: N %d N %d\n",N,N);

    // malloc host memory
    float *h_A, *h_x, *hostRef_xn, *gpuRef_xn, *h_b, h_norm;
    
    h_A = (float *)malloc(nBytes);
    h_x = (float *)malloc(nnBytes);
    h_b = (float *)malloc(nnBytes);
    hostRef_xn = (float *)malloc(nnBytes);
    gpuRef_xn = (float *)malloc(nnBytes);
    h_norm = 0;

    // initialize data at host side
    double iStart = seconds();
    initialData(h_A,  2.0f, nxy);
    initialData(h_x,  0.5f, N);
    
    double iElaps = seconds() - iStart;
    printf("Matrix initialization elapsed %f sec\n\n", iElaps);

    memset(hostRef_xn, 0, nnBytes);
    memset(gpuRef_xn, 0, nnBytes);
    memset(h_b, 0, nnBytes);


    float *d_A, *d_x, *d_b, *d_norm, *d_xn, *d_c;
    CHECK(cudaMalloc((void **)&d_A, nBytes));
    CHECK(cudaMalloc((void **)&d_x, nnBytes));
    CHECK(cudaMalloc((void **)&d_b, nnBytes));
    CHECK(cudaMalloc((void **)&d_xn, nnBytes));
    CHECK(cudaMalloc((void **)&d_norm, sizeof(float)));
    CHECK(cudaMalloc((void **)&d_c, nnBytes));

    CHECK(cudaMemcpy(d_A, h_A, nBytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_x, h_x, nnBytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemset(d_norm, 0, sizeof(float)));

    int dimx = 32;
    int dimy = 32;
    dim3 block(dimx, dimy);
    dim3 grid((N + block.x - 1) / block.x, (N+block.y - 1)/block.y);
    int threads = 1024;
    dim3 blocks  ((N + threads - 1) / threads);

    // invoke kernel at host side
    //int nThreads = THREADS_PER_BLOCK;
    //int nBlocks =  ((N + nThreads - 1) / nThreads);


    printf("******* Power iterative *******\n\n");
    // On Host
    /*iStart = seconds();
    powerIterationOnHost(h_A, h_x, hostRef_xn, N);
    iElaps = seconds() - iStart;
    printf("powerIterationOnHost elapsed %f sec\n\n", iElaps);
    */
    
    for (int i = 0; i < N; i++){

        hostRef_xn[i] = h_x[i];
    }
   
    for (int n = 0; n < 100; n++){

        matvecOnHost(h_A, hostRef_xn, h_b, N);
        
        EucldNormOnHost(h_b, &h_norm, N);
        
        VecNormalizationOnHost( h_b, hostRef_xn, &h_norm, N);
        
    }
    
    // On Device
    
    CHECK(cudaMemcpy(d_xn, d_x, N*sizeof(float), cudaMemcpyDeviceToDevice));
    
    for (int n = 0; n < 100; n++){

        matvecOnGPU<<<grid, block>>>(d_A, d_x, d_b, N);
        cudaDeviceSynchronize();
        EucldNormOnGPU<<<blocks, threads>>>(d_b, d_norm,d_c, N);
        cudaDeviceSynchronize();
        VecNormalizationOnGPU<<<blocks, threads>>>(d_b, d_x, d_norm, N);
        cudaDeviceSynchronize();
    }

    CHECK(cudaMemcpy(gpuRef_xn, d_x, nnBytes, cudaMemcpyDeviceToHost));
    checkResult(hostRef_xn, gpuRef_xn,N);
    
    CHECK(cudaFree(d_A));
    CHECK(cudaFree(d_x));
    CHECK(cudaFree(d_b));


    free(h_A);
    free(h_x);
    free(hostRef_xn);
    free(gpuRef_xn);
    

    CHECK(cudaDeviceReset());

    return (0);

}