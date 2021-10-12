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

  int N, i;

  // 1) --------- READ THE DATA ------------
  // make sure only rank 0 does this
 
  FILE *file_id;
  file_id = fopen(argv[1], "r");

  //read the first line which contains the number of data entries
  fscanf(file_id, "%d", &N);
  printf("N = %d\n",N);

  //allocate array for data
  double* x    = (double*) malloc(N*sizeof(double));

  //read N data elements
  int a;
  for (i=0; i<N; i++){
    fscanf(file_id,"%d",&a);
    x[i] = (double)a;
  }

  //colse the file
  fclose(file_id);

  //start timer
  double start_timer = MPI_Wtime();

  
  // 2) --------- LET ALL RANKS KNOW HOW MANY POINTS THERE ARE ------------

  // 3) --------- COMPUTE THE NUMBER OF POINTS PER RANK - N_LOC -----------

  // 4) --------- SEND EACH PROCESS A CHUNK OF DATA -----------------------

  // 5) --------- COMPUTE THE SUM OF DATA POINTS --------------------------
  // remember to use appropriate size of the array
  double mu = 0;
  for(i=0;i<N;i++){
    mu += x[i];
  }

  // 6) --------- COMPUTE THE GLOBAL SUM ----------------------------------

  // 7) --------- COMPUTE THE MEAN VALUE ----------------------------------
  mu = mu/N;

  // 8) --------- COMPUTE THE SUM OF DEVIATIONS ---------------------------
  //remember to use an appropriate loop size and variable name
  double sigma=0;
  for(i=0;i<N;i++){
    sigma += (x[i] - mu)*(x[i] - mu);
  }

  
  
  // 8) --------- COMPUTE GLOBAL SUM OF SIGMA_LOC -------------------------

  // 9) --------- COMPUTE THE STANDARD DEVIATIONS -------------------------
  //only on rank 0
  sigma = sqrt(sigma/N);

  //stop timer for computation
  double stop_timer = MPI_Wtime();

  double elapsed_time = stop_timer - start_timer;
  printf("Data size = %ld, mean = %f, std_dev = %f, elapsed_time = %f s \n",N,mu,sigma,elapsed_time);
  
}
