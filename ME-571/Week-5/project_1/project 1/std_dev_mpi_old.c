/***************************
std_dev reads in the data file provided as input argument and computes a mean and standard deviation.

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

void main(int argc, char* argv[]){

  MPI_Init(&argc, &argv);

  int N, i, nproc, irank, root_rank = 0;

  double* x;
  
  // 1) --------- READ THE DATA ------------
  // make sure only rank 0 does this
  
  //check communicator size
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  //Check rank number
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);

  if(irank == root_rank){
    // printf("reading from file: %s \n",argv[1]);

  FILE *file_id;
  file_id = fopen(argv[1], "r");
  
  //read the first line which contains the number of data entries
  fscanf(file_id, "%d", &N);
  //printf("N = %d\n",N);

  x = (double*) malloc(N*sizeof(double));

  //read N data elements
  int a;
  for (i=0; i<N; i++){
    fscanf(file_id,"%d",&a);
    x[i] = (double)a;
  }

  //colse the file
  fclose(file_id);
  }

  //start timer
  double start_timer = MPI_Wtime();
  
  // 2) --------- LET ALL RANKS KNOW HOW MANY POINTS THERE ARE -----------
  //start timer
  double timer_bcast_start = MPI_Wtime();

  MPI_Bcast(&N,1,MPI_INT,root_rank,MPI_COMM_WORLD);

  double timer_bcast = MPI_Wtime() - timer_bcast_start;
 
  // 3) --------- COMPUTE THE NUMBER OF POINTS PER RANK - N_LOC ----------
  double timer_comput_start = MPI_Wtime();

  int N_loc;
  N_loc = N/nproc;

  double* x_loc = (double*) malloc(N_loc*sizeof(double));
  
  double timer_comput = MPI_Wtime() - timer_comput_start;

  // 4) --------- SEND EACH PROCESS A CHUNK OF DATA ----------------------
  double timer_scatter_start = MPI_Wtime();
  MPI_Scatter(x,N_loc,MPI_DOUBLE,x_loc,N_loc,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);
  double timer_scatter = MPI_Wtime() - timer_scatter_start;


  // 5) --------- COMPUTE THE SUM OF DATA POINTS -------------------------
  timer_comput_start = MPI_Wtime();
  // remember to use appropriate size of the array
  double mu = 0,mu1 = 0;
  for(i=0;i<N_loc;i++){
    mu += x_loc[i];
  }
  timer_comput = timer_comput + (MPI_Wtime() - timer_comput_start);

  // 6) --------- COMPUTE THE GLOBAL SUM ---------------------------------
  double timer_allreduce_start = MPI_Wtime();

  MPI_Allreduce(&mu,&mu1,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  double timer_allreduce = MPI_Wtime() - timer_allreduce_start; 

  // 7) --------- COMPUTE THE MEAN VALUE ----------------------------------
  timer_comput_start = MPI_Wtime();

  mu1 = mu1/N;

  // 8) --------- COMPUTE THE SUM OF DEVIATIONS --------------------------
  //remember to use an appropriate loop size and variable name
  double sigma=0,sigma1 = 0;
  for(i=0;i<N_loc;i++){
    sigma += (x_loc[i] - mu1)*(x_loc[i] - mu1);
  }
  timer_comput = timer_comput + (MPI_Wtime() - timer_comput_start);

  // 8) --------- COMPUTE GLOBAL SUM OF SIGMA_LOC ------------------------
  double timer_reduce_start = MPI_Wtime();
  MPI_Reduce(&sigma,&sigma1,1,MPI_DOUBLE,MPI_SUM,root_rank,MPI_COMM_WORLD);

  double timer_reduce = MPI_Wtime() - timer_reduce_start;

  // 9) --------- COMPUTE THE STANDARD DEVIATIONS ------------------------
  //only on rank 0
  sigma1 = sqrt(sigma1/N);

  //stop timer for computation
  double stop_timer = MPI_Wtime();

  double elapsed_timer = stop_timer - start_timer;
  double timer_comm = timer_bcast + timer_scatter + timer_allreduce + timer_reduce;

  int ntimers = 7;
  double timing_array[7] = {elapsed_timer, timer_comput, timer_comm, timer_bcast, timer_scatter, timer_allreduce, timer_reduce};
  double timing_array_global[7];

  MPI_Reduce(timing_array, timing_array_global, ntimers, MPI_DOUBLE,MPI_MAX,root_rank,MPI_COMM_WORLD);

  if(irank == root_rank){
    printf("%d",nproc);
    for(i=0; i<ntimers; i++)
      {
	printf(",%f",timing_array[i]);
      }
    printf("\n");
  } 
  MPI_Finalize();
}
