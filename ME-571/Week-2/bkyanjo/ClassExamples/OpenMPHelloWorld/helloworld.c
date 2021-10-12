#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

//To set number of threads beforehand, use env var
//export OMP_NUM_THREADS=4
int main(int argc, char **argv) {
	#pragma omp parallel
	{
	printf("This will print a lot\n");
	int rank = omp_get_thread_num();
	int n = omp_get_max_threads();
        printf("Core #%d out of %d says hello world\n", rank, n);
	}
	printf("This should print only once\n");
}
