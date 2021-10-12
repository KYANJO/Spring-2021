#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "finite_difference_mpi.c"

//function headers (definitions at the end of the file)
double func(double x);
double func_first_der(double x); double func_second_der(double x);

int main(int argc, char * argv[]){

  double dx,ux_error;
  double *u,*u_x;
  double *u_x_ref;
  double *x;
  int npts, j, npts_loc;

  int irank, nproc;

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
  u_x_ref = malloc(sizeof(double)*npts_loc);
  u = malloc(sizeof(double)*npts_loc);
  u_x = malloc(sizeof(double)*npts_loc);

  // initialize arrays and reference solution
  for(j=0;j<npts_loc;j++){
    x[j] = (irank*npts_loc+j)*dx; // create x coordinates 
    u[j] = func(x[j]); //create function values
    u_x_ref[j] = func_first_der(x[j]); //create reference for first derivative
  }  

#ifdef DEBUG
  printf("rank=%d, npts =%d, npts_loc=%d\n",irank,npts,npts_loc);
  for(j=0;j<npts_loc;j++)
    printf("rank=%d, x[%d] = %f, u_x = %f\n",irank,j,x[j],u_x_ref[j]);
#endif  

  //Begin calculation
  
  SO_FirstDeriv_1D_MPI (npts_loc,dx,u,u_x,irank,nproc); //calc. 1st deriv. 
  
  // Computation of the L2 error 
  ux_error = L2_error_MPI(npts_loc,dx,u_x,u_x_ref); //check for ghosts
  if(irank==0)
    printf("%d\t%e\n",npts,ux_error);
  
  //Deallocation
  free(u);
  free(u_x);
  free(x);
  free(u_x_ref);

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
