/***************************
std_dev reads problem size N as input argument and computes a mean and standard deviation.

Written by: Michal A. Kopera
            Department of Mathematics
            Boise State University
            2/1/2021

 **************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


 //Data array using random number
 double random_number()
 {
   double random_value;
   random_value = (double)rand()/(double)RAND_MAX;
   return random_value;
 }

void main(int argc, char* argv[]){

  MPI_Init(&argc, &argv);

  int i, nproc, irank, root_rank = 0;

  long N = atol(argv[1]);

  double* x_loc;
  
  //check communicator size
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  //Check rank number
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);

  //start timer
  double start_timer = MPI_Wtime();
  
  //start timer
  double timer_bcast_start = MPI_Wtime();

  MPI_Bcast(&N,1,MPI_INT,root_rank,MPI_COMM_WORLD);

  double timer_bcast = MPI_Wtime() - timer_bcast_start;
 
  double timer_comput_start = MPI_Wtime();

  //Compute the number of points per rank - N_loc
  int N_loc;
  N_loc = N/nproc;

  x_loc = (double*) malloc(N_loc*sizeof(double));

  double timer_comput = MPI_Wtime() - timer_comput_start;

  //creatig an array of random numbers
  
  for(i=0; i<N_loc; i++)
    {
      x_loc[i] = (double)random_number();
    }
  
  timer_comput_start = MPI_Wtime();

  srand(irank+1);
  double mu = 0,mu1 = 0;
  for(i=0;i<N_loc;i++){
    mu += x_loc[i];
  }

  timer_comput = timer_comput + (MPI_Wtime() - timer_comput_start);

  //compute the global sum
  double timer_allreduce_start = MPI_Wtime();

  MPI_Allreduce(&mu,&mu1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  double timer_allreduce = MPI_Wtime() - timer_allreduce_start; 

  //compute the mean time
  timer_comput_start = MPI_Wtime();

  mu1 = mu1/N;

  //compute the sum of deviations
  double sigma=0,sigma1 = 0;
  for(i=0;i<N_loc;i++){
    sigma += (x_loc[i] - mu1)*(x_loc[i] - mu1);
  }

  timer_comput = timer_comput + (MPI_Wtime() - timer_comput_start);

  // COMPUTE GLOBAL SUM OF SIGMA_LOC ------------------------
  double timer_reduce_start = MPI_Wtime();
  MPI_Reduce(&sigma,&sigma1,1,MPI_DOUBLE,MPI_SUM,root_rank,MPI_COMM_WORLD);
  
  double timer_reduce = MPI_Wtime() - timer_reduce_start;

  // COMPUTE THE STANDARD DEVIATIONS ------------------------
  
  sigma1 = sqrt(sigma1/N);
  
//stop timer for computation
  double stop_timer = MPI_Wtime();

  double elapsed_timer = stop_timer - start_timer;
  double timer_comm = timer_bcast + timer_allreduce + timer_reduce;

  int ntimers = 6;
  double timing_array[6] = {elapsed_timer, timer_comput, timer_comm, timer_bcast, timer_allreduce, timer_reduce};
  double timing_array_global[6];

  MPI_Reduce(timing_array, timing_array_global, ntimers, MPI_DOUBLE,MPI_MAX,root_rank,MPI_COMM_WORLD);
 
  if(irank == root_rank){
    printf("%d, %e",nproc,sigma1);
    for(i=0; i<ntimers; i++)
      {
	printf(",%f",timing_array[i]);
      }
    printf("\n");
  }

  MPI_Finalize();
}
