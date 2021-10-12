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
            for (int iy = 0; iy < N; iy++)
            {   
                int id = iy*N+ix;
                b[ix] += A[id]*x[iy];
            }
        }

        //euclidean norm
        float cc = 0.0f;

       for (int i = 0; i < N; i++)
       {   
           cc += b[i]*b[i];
       }
        
       float normb = pow(cc,0.5);   
      
        //normalize
        for(int i = 0; i<N; i++)
        {
            V[i] = b[i]/(normb);
        }

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
    int nx  = N;
    int nByt = nx * sizeof(float);
    int nBytes = nxy * sizeof(float);
    //printf("Matrix size: nx %d ny %d\n", N, N);

    // malloc host memory
    float *A, *x, *xint, *V;
    A = (float *)malloc(nBytes);
    x = (float *)malloc(nByt);
    xint = (float *)malloc(nByt);
    V = (float *)malloc(nByt);

    double iStart = seconds();
    // initialize data at host side
    initialData(A,  2.0f, nxy);
    initialData(x,  0.5f, nx);
    initialData(xint,  0.5f, nx);

    memset(V, 0, nByt);
 
    powerIterationOnHost(xint, A, x, V, N);
    double iElaps = seconds() - iStart;
    //Elapsed time
    printf("N = %d\t Elapsed time = %f\n",N,iElaps);

    free(A);
    free(x);
    free(xint);
    free(V);

    return (0);
}
