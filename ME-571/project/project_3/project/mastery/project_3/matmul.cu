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

void matmul(float *A, float *B, float *C, const int N)
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

int main(int argc, char **argv)
{

    // set up data size of matrix
    int N = 1 << 5;
   
    int nxy = N * N;
    int nBytes = nxy * sizeof(float);
    //printf("Matrix size: nx %d ny %d\n", N, N);

    // malloc host memory
    float *A, *B, *C;
    A = (float *)malloc(nBytes);
    B = (float *)malloc(nBytes);
    C = (float *)malloc(nBytes);

    double iStart = seconds();
    // initialize data at host side
    initialData(A,  2.0f, nxy);
    initialData(B,  0.5f, nxy);

    memset(C, 0, nBytes);
 
    matmul(A, B, C, N);
    double iElaps = seconds() - iStart;
    //Elapsed time
    printf("N = %d\tElapsed time = %f\n",N,iElaps);

    free(A);
    free(B);
    free(C);

    return (0);
}
