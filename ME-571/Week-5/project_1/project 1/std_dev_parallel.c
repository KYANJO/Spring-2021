/***************************                                                    
std_dev reads problem size N as input argument and computes a mean a\
nd standard deviation.                                                        

Written by: Brian Kyanjo                                                    
            Department of Mathematics                                           
            Boise State University                                              
            2/15/2021                          
********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>  //include parallel library

//Data array using random number                                               
double random_number()
{
  double random_value;
  random_value = (double)rand()/(double)RAND_MAX;
  return random_value;
}

void main(int argc, char* argv[]){
  
  MPI_Init(&argc, &argv);

  int N_loc, nproc, irank, root_rank = 0;
  
  //check communicator size                                                   
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  //Check rank number                                                          
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);

  long N = atol(argv[1]);

  //start timer                                                                
  double start_timer = MPI_Wtime();
  //compute the number of points per rank - N_loc
  double timer_comput_start = MPI_Wtime();

  N_loc = N/nproc;

  //array of random numbers
  double* x =  (double*) malloc(N_loc*sizeof(double));
  int i;
  srand(irank+1);
  for (i = 0; i < N_loc; i++)
  {
    x[i] = (double)random_number();
  }

  double timer_comput = MPI_Wtime() - timer_comput_start;
  timer_comput_start = MPI_Wtime();

  // Welford's algorithm                       
  double* xbar =  (double*) malloc(N_loc*sizeof(double));
  double xb_0 = x[0]; xbar[0] = xb_0;
  double* M =  (double*) malloc(N_loc*sizeof(double));
  double M_glob;
  int k,Ni;
  double sigma,xbar_g,sigma1;

   for (k = 1; k<N_loc; k++)
   {
      xbar[k] = xbar[k-1] + (x[k] - xbar[k-1])/k;
      M[k] = M[k-1] + (x[k] - xbar[k])*(x[k] - xbar[k-1]);
   }

   double* xbar_glob =(double*) malloc(nproc*sizeof(double));
   double xbar_k=xbar[N_loc-1];

   timer_comput = timer_comput + (MPI_Wtime() - timer_comput_start);
   
   //compute the glbal mean
   double timer_gather_start = MPI_Wtime();

   MPI_Gather(&xbar_k,1,MPI_DOUBLE,xbar_glob,1,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);

   double timer_gather = MPI_Wtime() - timer_gather_start;

   //compute the global square sum
   double timer_reduce_start = MPI_Wtime();

   MPI_Reduce(&M[N_loc-1],&M_glob,1,MPI_DOUBLE,MPI_SUM,root_rank,MPI_COMM_WORLD);

   double timer_reduce = MPI_Wtime() - timer_reduce_start;

   timer_comput_start = MPI_Wtime();   
   double xbar_o;
   
     xbar_o = xbar_glob[0];
     for (i = 1; i<nproc; i++)
     {
       Ni = (i+1)*N_loc;
       xbar_g = xbar_o + (-xbar_o + xbar_glob[i])*(i*N_loc)/Ni;
       M_glob += (pow((xbar_o - xbar_glob[i]),2))*(i*pow(N_loc,2)/Ni);
       xbar_o = xbar_g;
    }
     
   timer_comput = timer_comput + (MPI_Wtime() - timer_comput_start);

  //compute the standard deviation 
   sigma = sqrt(M_glob/N);

   //stop timer for computation                                                
   double stop_timer = MPI_Wtime(); 

   double elapsed_timer = stop_timer - start_timer;
   double timer_comm = timer_gather + timer_reduce;

   int ntimers = 5;
   double timing_array[5] = {elapsed_timer, timer_comput, timer_comm, timer_gather, timer_reduce};
   double timing_array_global[5];

   MPI_Reduce(timing_array, timing_array_global, ntimers, MPI_DOUBLE,MPI_MAX,root_rank,MPI_COMM_WORLD);

   if(irank == root_rank)
   {
     printf("nproc = %d, mu = %e, sigma = %e",nproc, xbar_o,sigma);
     for(i=0; i<ntimers; i++)
       {
	 printf(",%f",timing_array[i]);
       }
     printf("\n");
   }
  
    MPI_Finalize();
}
