#include <stdio.h>
#include <mpi.h> //Include MPI library

int main(int argc, char *argv[])
{

  //Initialize MPI - create all available processes
  MPI_Init(&argc, &argv);

  //Variables to store communicator size and rank number
  int nproc, irank; 

  //Check communicator size
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  //Check rank number
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);

  //Print hello message from each rank
  printf("Hello, World! from rank %d of %d\n",irank,nproc);

  //Finalize MPI - close all ranks
  MPI_Finalize();
  return 0;

}
