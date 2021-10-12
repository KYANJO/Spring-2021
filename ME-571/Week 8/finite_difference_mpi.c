#include<mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void SO_FirstDeriv_1D_MPI_v1(int npts_loc, double dx, double *u,
			  double *dudx, int irank, int nproc){ 
  double two_invdx = 1.0/(2.0*dx);
  double ghost_left, ghost_right;
  MPI_Status status;

  //Interior points
  for(int i=1;i<npts_loc-1;i++) //middle points (not boundary or ghosts)
    dudx[i] = (u[i+1]-u[i-1])*two_invdx;

  if(irank>0){
    //receive ghost from left
    MPI_Recv(&ghost_left,1,MPI_DOUBLE,irank-1,101, MPI_COMM_WORLD, &status);

    //compute left rank-boundary
    dudx[0] = (u[1]-ghost_left)*two_invdx;
  }

  if(irank==0){
    //compute left point using one-sided formula
    dudx[0] = (-3.0*u[0] + 4.0*u[1] - u[2])*two_invdx; 
  }

  if(irank<nproc-1) {
    //send ghost to right
    ghost_right = u[npts_loc-1];
    MPI_Send(&ghost_right,1,MPI_DOUBLE,irank+1,101,MPI_COMM_WORLD);

    //receive ghost from right
    MPI_Recv(&ghost_right,1,MPI_DOUBLE,irank+1,102,MPI_COMM_WORLD, &status); 
    
    //compute right rank-boundary
    dudx[npts_loc-1] = (ghost_right-u[npts_loc-2])*two_invdx;
  }

  if(irank==nproc-1) {
    dudx[npts_loc-1] = (3.0*u[npts_loc-1] - 4.0*u[npts_loc-2] + u[npts_loc-3])*two_invdx; 
  }

  if(irank>0){
    //send ghost to left
    ghost_left = u[0]; 
    MPI_Send(&ghost_left,1,MPI_DOUBLE,irank-1,102,MPI_COMM_WORLD);

  }
  
   return; 
}

void SO_FirstDeriv_1D_MPI_v2(int npts_loc, double dx, double *u,
			  double *dudx, int irank, int nproc){ 
  double two_invdx = 1.0/(2.0*dx);
  double ghost_left, ghost_right;
  MPI_Status status;

  //Interior points
  for(int i=1;i<npts_loc-1;i++) //middle points (not boundary or ghosts)
    dudx[i] = (u[i+1]-u[i-1])*two_invdx;

  if(irank>0){
    //receive ghost from left
    MPI_Recv(&ghost_left,1,MPI_DOUBLE,irank-1,101, MPI_COMM_WORLD, &status);
  }

  if(irank<nproc-1) {
    //send ghost to right
    ghost_right = u[npts_loc-1];
    MPI_Send(&ghost_right,1,MPI_DOUBLE,irank+1,101,MPI_COMM_WORLD);
  }

  if(irank==0){
    //compute left point using one-sided formula
    dudx[0] = (-3.0*u[0] + 4.0*u[1] - u[2])*two_invdx; 
  }else{
    //compute left rank-boundary
    dudx[0] = (u[1]-ghost_left)*two_invdx;
  }

  if(irank<nproc-1) {
    //receive ghost from right
    MPI_Recv(&ghost_right,1,MPI_DOUBLE,irank+1,102,MPI_COMM_WORLD, &status); 
  }
  if(irank>0){
    //send ghost to left
    ghost_left = u[0]; 
    MPI_Send(&ghost_left,1,MPI_DOUBLE,irank-1,102,MPI_COMM_WORLD);
  }

  if(irank==nproc-1) {
    dudx[npts_loc-1] = (3.0*u[npts_loc-1] - 4.0*u[npts_loc-2] + u[npts_loc-3])*two_invdx; 
  }else{
    //compute right rank-boundary
    dudx[npts_loc-1] = (ghost_right-u[npts_loc-2])*two_invdx;
  }
    
   return; 
}

void SO_FirstDeriv_1D_MPI_v3(int npts_loc, double dx, double *u,
			  double *dudx, int irank, int nproc){ 
  double two_invdx = 1.0/(2.0*dx);
  double ghost_left, ghost_right;
  MPI_Status status;

  //Interior points
  for(int i=1;i<npts_loc-1;i++) //middle points (not boundary or ghosts)
    dudx[i] = (u[i+1]-u[i-1])*two_invdx;

  //Communicate to the right
  ghost_right = u[npts_loc-1];
  if(irank>0 && irank<nproc-1){
    MPI_Sendrecv(&ghost_right,1,MPI_DOUBLE,irank+1,101,&ghost_left,1,MPI_DOUBLE,irank-1,101,MPI_COMM_WORLD,&status);
  }
  else if(irank==0 && nproc>1){
    MPI_Send(&ghost_right,1,MPI_DOUBLE,irank+1,101,MPI_COMM_WORLD);
  }
  else if(irank==nproc-1 && nproc>1){
    MPI_Recv(&ghost_left,1,MPI_DOUBLE,irank-1,101, MPI_COMM_WORLD, &status);
  }

  //Compute left end-points
  if(irank==0){
    //compute left domain boundary point using one-sided formula
    dudx[0] = (-3.0*u[0] + 4.0*u[1] - u[2])*two_invdx; 
  }else{
    //compute left rank-boundary
    dudx[0] = (u[1]-ghost_left)*two_invdx;
  }

  //Communicate to the left
  ghost_left = u[0];
  if(irank>0 && irank<nproc-1){
    MPI_Sendrecv(&ghost_left,1,MPI_DOUBLE,irank-1,102,&ghost_right,1,MPI_DOUBLE,irank+1,102,MPI_COMM_WORLD,&status);
  }
  else if(irank==0 && nproc>1){
    MPI_Recv(&ghost_right,1,MPI_DOUBLE,irank+1,102,MPI_COMM_WORLD, &status);
  }
  else if(irank==nproc-1 && nproc>1){
    MPI_Send(&ghost_left,1,MPI_DOUBLE,irank-1,102, MPI_COMM_WORLD);
  }

  //Compute right end-points
  if(irank==nproc-1) {
    //compute at the domain boundary using one-sided formula
    dudx[npts_loc-1] = (3.0*u[npts_loc-1] - 4.0*u[npts_loc-2] + u[npts_loc-3])*two_invdx; 
  }else{
    //compute right rank-boundary
    dudx[npts_loc-1] = (ghost_right-u[npts_loc-2])*two_invdx;
  }
    
   return; 
}

void SO_FirstDeriv_1D_MPI(int npts_loc, double dx, double *u,
			  double *u_x, int irank, int nproc){ 
  double two_invdx = 1.0/(2.0*dx);
  double ghost;
  MPI_Status status;

  //Interior points
  for(int i=1;i<npts_loc-1;i++) //middle points (not boundary or ghosts)
    u_x[i] = (u[i+1]-u[i-1])*two_invdx;

  // Physical boundaries
  //left
  if(irank == 0)
    u_x[0] = (-3.0*u[0] + 4.0*u[1] - u[2])*two_invdx; 

  //right
  if(irank == (nproc-1))
    u_x[npts_loc-1] = (3.0*u[npts_loc-1] - 4.0*u[npts_loc-2] + u[npts_loc-3])*two_invdx; 

  // Inter-process boundary
  if(irank == 0){ //first rank communicates only to the right
    ghost = u[npts_loc-1]; //right ghost point
    MPI_Send(&ghost,1,MPI_DOUBLE,1,101,MPI_COMM_WORLD); //send ghost to rank 1
    MPI_Recv(&ghost,1,MPI_DOUBLE,1,102,MPI_COMM_WORLD, &status); //receive ghost from rank 1
    u_x[npts_loc-1] = (ghost - u[npts_loc-2])*two_invdx; //compute derivative at rank boundary 
  }
  else if(irank == (nproc-1)){ //Last rank communicates only to the left
    //receive ghost from left
    MPI_Recv(&ghost,1,MPI_DOUBLE,irank-1,101,MPI_COMM_WORLD, &status);
    //compute left rank-boundary
    u_x[0] = (u[1]-ghost)*two_invdx;

    //send ghost to left
    ghost = u[0];
    MPI_Send(&ghost,1,MPI_DOUBLE,irank-1,102,MPI_COMM_WORLD);

  } else{
    //receive ghost from left
    MPI_Recv(&ghost,1,MPI_DOUBLE,irank-1,101, MPI_COMM_WORLD, &status);
    //compute left rank-boundary
    u_x[0] = (u[1]-ghost)*two_invdx;

    //send ghost to left
    ghost = u[0]; 
    MPI_Send(&ghost,1,MPI_DOUBLE,irank-1,102,MPI_COMM_WORLD);

    //send ghost to right
    ghost = u[npts_loc-1];
    MPI_Send(&ghost,1,MPI_DOUBLE,irank+1,101,MPI_COMM_WORLD);

    //receive ghost from right
    MPI_Recv(&ghost,1,MPI_DOUBLE,irank+1,102,MPI_COMM_WORLD, &status); 
    //compute right rank-boundary
    u_x[npts_loc-1] = (ghost-u[npts_loc-2])*two_invdx;
  }
  return; 
}

void SO_FirstDeriv_1D_MPI_SR(int npts_loc, double dx, double *u,
			  double *u_x, int irank, int nproc){ 
  double two_invdx = 1.0/(2.0*dx);
  double ghost;
  MPI_Status status;

  //Interior points
  for(int i=1;i<npts_loc-1;i++) //middle points (not boundary or ghosts)
    u_x[i] = (u[i+1]-u[i-1])*two_invdx;

  // Physical boundaries
  //left
  if(irank == 0)
    u_x[0] = (-3.0*u[0] + 4.0*u[1] - u[2])*two_invdx; 

  //right
  if(irank == (nproc-1))
    u_x[npts_loc-1] = (3.0*u[npts_loc-1] - 4.0*u[npts_loc-2] + u[npts_loc-3])*two_invdx; 

  //Inter-process boundary
  if(irank == 0){//first rank communicates only to the right
    ghost = u[npts_loc-1]; //right ghost point

    //send ghost to right neighbor and receive ghost from it
    MPI_Sendrecv_replace(&ghost,1,MPI_DOUBLE,1,101,1,102,MPI_COMM_WORLD,&status);

    //compute derivative at right rank boundary 
    u_x[npts_loc-1] = (ghost - u[npts_loc-2])*two_invdx; 
  }
  else if(irank == (nproc-1)){ //Last rank communicates only to the left

    ghost = u[0]; //left ghost point
    //send ghost to left and receive ghost from it
    MPI_Sendrecv_replace(&ghost,1,MPI_DOUBLE,irank-1,102,irank-1,101,MPI_COMM_WORLD,&status);

    //compute derivative at left rank-boundary
    u_x[0] = (u[1]-ghost)*two_invdx;

  } else{

    //left ghost
    ghost = u[0]; 
    //send ghost to left rank receive ghost from left rank
    MPI_Sendrecv_replace(&ghost,1,MPI_DOUBLE,irank-1,102,irank-1,101,MPI_COMM_WORLD, &status);
    //compute left rank-boundary
    u_x[0] = (u[1]-ghost)*two_invdx;


    //right ghost 
    ghost = u[npts_loc-1];
    //send ghost to right and revceive ghost from right
    MPI_Sendrecv_replace(&ghost,1,MPI_DOUBLE,irank+1,101,irank+1,102,MPI_COMM_WORLD,&status);
    //compute right rank-boundary
    u_x[npts_loc-1] = (ghost-u[npts_loc-2])*two_invdx;
  }
  return; 
}


void SO_FirstDeriv_1D_MPI_SRv1(int npts_loc, double dx, double *u,
			  double *u_x, int irank, int nproc){ 
  double two_invdx = 1.0/(2.0*dx);
  double ghost;
  MPI_Status status;

  //Interior points
  for(int i=1;i<npts_loc-1;i++) //middle points (not boundary or ghosts)
    u_x[i] = (u[i+1]-u[i-1])*two_invdx;

  // Physical boundaries
  //left
  if(irank == 0)
    u_x[0] = (-3.0*u[0] + 4.0*u[1] - u[2])*two_invdx; 

  //right
  if(irank == (nproc-1))
    u_x[npts_loc-1] = (3.0*u[npts_loc-1] - 4.0*u[npts_loc-2] + u[npts_loc-3])*two_invdx; 

  //left side communication
  if(irank!=0){
    ghost = u[0]; 
    //send ghost to left rank receive ghost from left rank
    MPI_Sendrecv_replace(&ghost,1,MPI_DOUBLE,irank-1,102,irank-1,101,MPI_COMM_WORLD, &status);
    //compute left rank-boundary
    u_x[0] = (u[1]-ghost)*two_invdx;
  }

  //right side communication
  if(irank!=(nproc-1)){
    ghost = u[npts_loc-1];
    //send ghost to right and revceive ghost from right
    MPI_Sendrecv_replace(&ghost,1,MPI_DOUBLE,irank+1,101,irank+1,102,MPI_COMM_WORLD,&status);
    //compute right rank-boundary
    u_x[npts_loc-1] = (ghost-u[npts_loc-2])*two_invdx;
  }
  return; 
}

double L2_error_MPI(int npts_loc, double dx, double* data, double* reference){
    double error=0.0; 
    double error_total = 0.0;

    for(int j=0;j<npts_loc;j++){
      error += dx*pow((data[j]-reference[j]),2);
    }
    MPI_Allreduce(&error,&error_total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return sqrt(error_total);
}

int getArgs_mpi(int argc, char *argv[], int irank, MPI_Comm comm)
{
  int npts;

  if(irank==0){
    if ( argc != 2 ) /* argc should be 2 for correct execution */
      {
	//If not correct number of arguments, assume n=1000
	printf("Incorrect number of arguments. Usage: ./derivation N \n Assuming N=1024.\n");
	npts=1024;
      }
    else
      {
	//Use input argument
	npts = atoi(argv[1]);	
      }
  }
  MPI_Bcast(&npts,1,MPI_INT,0,comm);
  return npts;
}
