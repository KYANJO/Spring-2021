#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "finite_difference_mpi.c"

//function headers (definitions at the end of the file)
double func(double x);
double func_first_der(double x); double func_second_der(double x);

int main(int argc, char * argv[]){

  double dx,dudx_error_v1,dudx_error_v2,dudx_error_v3;
  double *u,*dudx;
  double *dudx_ref;
  double *x;
  int npts, j, npts_loc;

  int irank, nproc;

  int root_rank = 0;
  
  MPI_Init(&argc,&argv); //initialize MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);

  npts = getArgs_mpi(argc,argv,irank,MPI_COMM_WORLD); //get command-line argument

  //set dx based on number of points
  dx = 1.0/(npts-1); 
  
  //decide how many points is owned by local rank
  npts_loc = npts/nproc;

  //set-up coordinates and reference solutions
  x = malloc(sizeof(double)*npts_loc);
  dudx_ref = malloc(sizeof(double)*npts_loc);
  u = malloc(sizeof(double)*npts_loc);
  dudx = malloc(sizeof(double)*npts_loc);

  // initialize arrays and reference solution
  for(j=0;j<npts_loc;j++){
    x[j] = (irank*npts_loc+j)*dx; // create x coordinates 
    u[j] = func(x[j]); //create function values
    dudx_ref[j] = func_first_der(x[j]); //create reference for first derivative
  }  

#ifdef DEBUG
  printf("rank=%d, npts =%d, npts_loc=%d\n",irank,npts,npts_loc);
  for(j=0;j<npts_loc;j++)
    printf("rank=%d, x[%d] = %f, dudx_ref = %f\n",irank,j,x[j],dudx_ref[j]);
#endif  

  //Begin calculation
  
  //  SO_FirstDeriv_1D_MPI_v1(npts_loc,dx,u,dudx,irank,nproc); //calc. 1st deriv.

  
  double time_v1_local, time_v2_local, time_v3_local;
  for(j=0;j<2;j++){ //repeat each experiment two times and take the second result for timing
    time_v1_local = MPI_Wtime();
    SO_FirstDeriv_1D_MPI_v1(npts_loc,dx,u,dudx,irank,nproc); //calc. 1st deriv. 
    time_v1_local = MPI_Wtime() - time_v1_local;
    // Computation of the L2 error 
    dudx_error_v1 = L2_error_MPI(npts_loc,dx,dudx,dudx_ref); //check error
  }
  
  for(j=0;j<2;j++){ //repeat each experiment two times and take the second result for timing
    time_v2_local = MPI_Wtime();
    SO_FirstDeriv_1D_MPI_v2(npts_loc,dx,u,dudx,irank,nproc); //calc. 1st deriv. 
    time_v2_local = MPI_Wtime() - time_v2_local;
    // Computation of the L2 error 
    dudx_error_v2 = L2_error_MPI(npts_loc,dx,dudx,dudx_ref); //check error
  }
  
  for(j=0;j<2;j++){ //repeat each experiment two times and take the second result for timing
    time_v3_local = MPI_Wtime();
    SO_FirstDeriv_1D_MPI_v3(npts_loc,dx,u,dudx,irank,nproc); //calc. 1st deriv. 
    time_v3_local = MPI_Wtime() - time_v3_local;
    // Computation of the L2 error 
    dudx_error_v3 = L2_error_MPI(npts_loc,dx,dudx,dudx_ref); //check error
  }

  //Wrap up timing and print results
  int ntimers = 3;
  double timing_array_local[3] = {time_v1_local,time_v2_local,time_v3_local};
  double timing_array_global[3];

  double error_array[3] = {dudx_error_v1,dudx_error_v2,dudx_error_v3};
  
  MPI_Reduce(timing_array_local, timing_array_global, ntimers, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);
  
  if(irank==root_rank){
    printf("%d\t%d",nproc,npts);
    for(j=0; j<ntimers; j++){
      printf(", \t%f\t%e",timing_array_global[j],error_array[j]); 
    }
    printf("\n");
  }
  
  //Deallocation
  free(u);
  free(dudx);
  free(x);
  free(dudx_ref);

  MPI_Finalize();

  return 0;
} 


double func(double x){
  return(x*x*x*x);
}
double func_first_der(double x){ 
  return(4*x*x*x);
}
double func_second_der(double x){
  return(12*x*x);
}
