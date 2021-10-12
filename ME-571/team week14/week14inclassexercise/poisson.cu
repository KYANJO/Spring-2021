#include "common.h"
#include <cuda_runtime.h>
#include <stdio.h>
#include <math.h>

#include "kernels.cu"

#define ID_2D(i,j,nx) ((j)*(nx+2)+(i))


__global__ void test_kernels(float *u, float *u_new, float *f, int N)
{
  unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;

  unsigned int idx = j * N + i;

  if(i<N && j<N) {
    float h = 1/(N-1);
    float x = i*h;
    float y = j*h;
    float ue = 0;
    float fe = -8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y);
    float ue_new = -0.25*h*h*fe;
    printf("Thread (%d, %d), Block (%d,%d), u = %f, u_new = %f, f = %f, expected values: %f, %f, %f\n",threadIdx.x,threadIdx.y,blockIdx.x,blockIdx.y,u[idx],u_new[idx],f[idx],ue,ue_new,fe);
  }
  
}

double L2_error(int N, double L, double* data, double* reference){

    double error=0.0; 
      int id;
      
      for(int j=1;j<N+1;j++){
        for(int i=1;i<N+1;i++){
            id = ID_2D(i,j,N);
            error += pow((data[id]-reference[id]),2);
        }
      }
  
      return sqrt(error)/N;
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
    int N = 6;
    int i, j;
    // define domain size (assume both x and y directions are the same)
    double L = 1.0;
    //compute the grid spacing
    double h = L/(N-1);
    int nxy = N * N;
    int nx = N;
    int ny = N;
    int nBytes = nxy * sizeof(float);
    printf("Problem size: nx %d ny %d\n", N, N);
    
    //allocate on host
    float *h_u, *h_u_exact, *h_error, *x, *y;
    h_u = (float *)malloc(nBytes);
    h_u_exact = (float *)malloc(nBytes);
    h_error = (float *)malloc(sizeof(float));
    x       = (float *)malloc((N)*sizeof(float));
    y       = (float *)malloc((N)*sizeof(float));

    // malloc device global memory
    float *d_u, *d_f, *d_u_new, *d_error;
    CHECK(cudaMalloc((void **)&d_u     , nBytes));
    CHECK(cudaMalloc((void **)&d_f     , nBytes));
    CHECK(cudaMalloc((void **)&d_u_new , nBytes));
    CHECK(cudaMalloc((void **)&d_error     ,sizeof(float)));

    int id;
    //initialize
    for(i=0;i<N;i++){
        x[i] = i*h;
        y[i] = i*h;//y coordinates are the same as x, but need not be
      }
      
    //initialize array u
    for(j=0;j<N;j++){
        for(i=0;i<N;i++){
            
            id = ID_2D(i,j,N);
            h_u_exact[id] = sin(2*M_PI*x[i])*sin(2*M_PI*y[j]);
            // h_f[id] = -8*M_PI*M_PI*sin(2*M_PI*x[i])*sin(2*M_PI*y[j]);
        }
    }

    //set-up blocks and threads
    int dimx = 32;
    int dimy = 32;
    dim3 block(dimx, dimy);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
    
    //initialize arrays
    initData<<<grid,block>>>(d_u, d_f, N);
    CHECK(cudaDeviceSynchronize());

    //initialize error
    h_error[0] = 100.0;
    double epsilon = 1e-8;

    while(h_error[0] > epsilon){
        //take one Poisson update 
        poissonKernel<<<grid, block>>>(d_u, d_u_new, d_f, N);
        CHECK(cudaDeviceSynchronize());
        computeError<<<grid, block>>>(d_error, d_u, d_u_new, N);
        CHECK(cudaDeviceSynchronize());
        updateSolution<<<grid, block>>>(d_u, d_u_new, N);
        CHECK(cudaDeviceSynchronize());
        CHECK(cudaMemcpy(h_error, d_error, sizeof(float), cudaMemcpyDeviceToHost));
    }

    CHECK(cudaMemcpy(h_u, d_u, nBytes, cudaMemcpyDeviceToHost));

    L2_error(N, L, h_u, h_u_exact);
    // check kernel error
    CHECK(cudaGetLastError());
    test_kernels<<<grid, block>>>(d_u, d_u_new, d_f, N);
    // free device global memory
    CHECK(cudaFree(d_u));
    CHECK(cudaFree(d_f));
    CHECK(cudaFree(d_u_new));
    CHECK(cudaFree(d_error));

    free(h_u);
    free(h_u_exact);
    free(h_error);
    free(x);
    free(y);

    // reset device
    CHECK(cudaDeviceReset());

    return (0);
}
