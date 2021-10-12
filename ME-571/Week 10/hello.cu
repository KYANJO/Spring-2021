//Just your regular Hello World file
// to be compiled with nvcc rather than gcc

#include <stdio.h>

__global__ void helloFromGPU(void) {
  printf("Hello World from GPU, block %d, thread %d\n",blockIdx.x,threadIdx.x);
}


int main(void) {
  printf("Hello World from CPU!\n");

  helloFromGPU<<<1, 10>>>();
  cudaDeviceReset();

  return 0;
}
