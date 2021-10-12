#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char **argv) {
	int n, rank;
	MPI_Status status;
	int error = MPI_Init(&argc,&argv);
	printf("This will print a lot\n");
	MPI_Comm_size(MPI_COMM_WORLD, &n);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        printf("Core #%d out of %d says hello world\n", rank, n);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		printf("This should only print once\n");
		MPI_Finalize();
	}
	return 0;
}
