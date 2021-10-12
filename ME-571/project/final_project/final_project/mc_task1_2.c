/*
Author: Brian KYANJO
Date:   May 1st, 2021
Class:  ME571

Description:
------------
Monte Carlo integration implementation using MPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
                                       
 //Data array using random number                                               
 double random_number()
 {
   double random_value;
   random_value = (double)rand()/(double) RAND_MAX;
   return random_value;
 }

void main(int argc, char* argv[]){

    MPI_Init(&argc, &argv);

    double lam = 1.0;

    int i, nproc, irank, root_rank = 0;

    long N = atol(argv[1]);

    double* x_loc;

    //check communicator size                                                     
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    //Check rank number                                                           
    MPI_Comm_rank(MPI_COMM_WORLD, &irank);

    //Compute the number of points per rank - N_loc                               
    int N_loc;
    N_loc = N/nproc;

    x_loc = (double*) malloc(N_loc*sizeof(double));

    //creatig an array of random numbers                                          

    for(i=0; i<N_loc; i++)
        {
        x_loc[i] = (double)random_number();
        }

    srand(time(NULL));
    double u = 0,u1 = 0;
    for(i=0;i<N_loc;i++){
        u += cos(-log(x_loc[i]));
    }
     
    MPI_Reduce(&u,&u1,1,MPI_DOUBLE,MPI_SUM,root_rank,MPI_COMM_WORLD);

    u1 = u1/N;

    //Exact value
    double uexact = lam/(1+pow(lam,2));
    double error  =  uexact - u1;

    if(irank == root_rank){
        printf("nproc = %d\t uexact = %f\t umc = %f\t error = %e\n",nproc,uexact,u1,error);
    }

    MPI_Finalize();

}