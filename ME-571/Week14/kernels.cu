__global__ void initData(float *u, float *f, int N)
{

    printf("Hello from initData: Thread (%d, %d), Block (%d,%d)\n",threadIdx.x,threadIdx.u,blockIdx.x,blockIdx.y);
}

__global__ void poissonKernel(float *u, float *u_new, float *f, int N)
{

    printf("Hello from poissonKernel: Thread (%d, %d), Block (%d,%d)\n",threadIdx.x,threadIdx.u,blockIdx.x,blockIdx.y);
}

__global__ void computeError(float error, float *u, float *u_new, int N)
{

    printf("Hello from computeError: Thread (%d, %d), Block (%d,%d)\n",threadIdx.x,threadIdx.u,blockIdx.x,blockIdx.y);
}


__global__ void updateSolution(float *u, float *u_new, int N)
{

    printf("Hello from updateSolution: Thread (%d, %d), Block (%d,%d)\n",threadIdx.x,threadIdx.u,blockIdx.x,blockIdx.y);
}
