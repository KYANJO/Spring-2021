/*
The code solves a 2D Poisson equation in a square domain using finite difference method

Written by: Michal A. Kopera
            Department of Mathematics
	    Boise State University
	    03/06/2021
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define ID_2D(i,j,nx) ((j)*(nx+2)+(i))

// *** Function interfaces ***

// getArgs_mpi reads input parameters Nfrom command line
void getArgs(int *N, int argc, char *argv[], int irank, MPI_Comm comm);

//functions f and g
double* f(int N, double* u, double* v);

double* g(int N, double* u, double* v);

// write results to file
void write_results(double *u, int N, int N_loc, int nproc, int irank, double h, double L);

// *** Main function ***
int main(int argc, char * argv[]){

  int N, i, j, id, id_ghost;
  int id_lft, id_rgt, id_top, id_bot;
  
  int irank,nproc;
  
  //initialize MPI
  MPI_Init(&argc,&argv); 
  MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  
   //compute row and column index in the rank grid
  int q = (int)sqrt(nproc); //assume nproc is a square number
  int rank_row = irank/q;
  int rank_col = irank%q;

  double eps = 0.02, a =0.75, b = 0.01;

  //read command line arguments
  getArgs(&N, argc, argv, irank, MPI_COMM_WORLD);

  // define domain size (assume both x and y directions are the same)
  double L = 20.0;
 
  //tfinal
  int tfinal = 20;
  
  //compute the grid spacing
  double h = 2*L/(N-1);
  
   //compute local problem size
  double N_loc = N/q;
  
  //allocate arrays
  double *u       = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *v       = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *f       = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *g       = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *u_new   = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *v_new   = malloc((N_loc+2)*(N_loc+2)*sizeof(double));
  double *x       = malloc((N_loc+2)*sizeof(double));
  double *y       = malloc((N_loc+2)*sizeof(double));
  
  //allocate buffers
  double *send_buffer_rgt = malloc(N_loc*sizeof(double));
  double *send_buffer_lft = malloc(N_loc*sizeof(double));
  double *send_buffer_bot = malloc(N_loc*sizeof(double));
  double *send_buffer_top = malloc(N_loc*sizeof(double));

  double *recv_buffer_rgt = malloc(N_loc*sizeof(double));
  double *recv_buffer_lft = malloc(N_loc*sizeof(double));
  double *recv_buffer_bot = malloc(N_loc*sizeof(double));
  double *recv_buffer_top = malloc(N_loc*sizeof(double));
  
  double *send_buffer_rgt_v = malloc(N_loc*sizeof(double));
  double *send_buffer_lft_v = malloc(N_loc*sizeof(double));
  double *send_buffer_bot_v = malloc(N_loc*sizeof(double));
  double *send_buffer_top_v = malloc(N_loc*sizeof(double));

  double *recv_buffer_rgt_v = malloc(N_loc*sizeof(double));
  double *recv_buffer_lft_v = malloc(N_loc*sizeof(double));
  double *recv_buffer_bot_v = malloc(N_loc*sizeof(double));
  double *recv_buffer_top_v = malloc(N_loc*sizeof(double));

  //initialize x and y
  for(i=1;i<N_loc+1;i++){
    x[i] = (i-(N_loc+1)/2.0)*h + rank_col*N_loc*h;
    y[i] = (i-(N_loc+1)/2.0)*h + rank_col*N_loc*h;//y coordinates are the same as x, but need not be
  }
  
  //initialize array u
  for(j=1;j<N_loc+1;j++){
      for(i=1;i<N_loc+1;i++){
	id = ID_2D(i,j,N_loc);
	
	if(y[j]>=0){
	  u[id] = 1;
	}
	else if(y[j]<0){
	  u[id] = 0;
	}
	
	if(x[i]>=0){
	  v[id] = 1;
        }
	else if(x[i]<0){
          v[id] = 0;
        }
      }
  }

  for(int j=1;j<N_loc+1;j++){
    for(int i=1;i<N_loc+1;i++){
      id = ID_2D(i,j,N_loc);
      f[id] = (1.0/eps)*(u[id])*(1.0 - u[id])*(u[id] - ((v[id] + b)/a));
      g[id] = u[id] - v[id];
    }
  }

  //initialize BC in ghost cells
  if(rank_row==0){
   for(i=1;i<N_loc+1;i++){
    id = ID_2D(i,0,N_loc); //bottom ghost
    int idn =ID_2D(i,1,N_loc);
    u[id] = u[idn];
    }
    }
    
    if(rank_row==q-1){
    for(i=1;i<N_loc+1;i++){
    id = ID_2D(i,N_loc+1,N_loc); //top ghost
    int idn =ID_2D(i,N_loc,N_loc);
    u[id] = u[idn];
    }
    }

    if(rank_col==0){
    for(i=1;i<N_loc+1;i++){
    id = ID_2D(0,i,N_loc); //left ghost
    int idn = ID_2D(1,i,N_loc);
    u[id] = u[idn];
    }
    }
    
    if(rank_row==q-1){
    for(i=1;i<N_loc+1;i++){
    id = ID_2D(N_loc+1,i,N_loc); //right ghost
    int idn =ID_2D(N_loc,i,N_loc);
    u[id] = u[idn];
    }
    }

  double dt = 0.0001;
  int ntime = (int)tfinal/dt;
  for(int time=1;time<=ntime;time++){
    
    //go over interior points 2:N-1
    for(j=2;j<N_loc;j++){
      for(i=2;i<N_loc;i++){
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top

	u_new[id] = u[id] +  (dt/pow(h,2))*(f[id] + u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);	
	v_new[id] = v[id] +  (dt/pow(h,2))*g[id];
      }
    }
    
    //enforce Dirichlet BC by copying values from ghosts
    if(rank_row ==0){
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,1,N_loc); //bottom boundary
      id_ghost = ID_2D(i,0,N_loc); //bottom ghost
      u_new[id] = u[id_ghost];
      v_new[id] = v[id_ghost];  
    } 
    }
    
    
      if(rank_row == q-1){
      for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,N_loc,N_loc); //top boundary
      id_ghost = ID_2D(i,N_loc+1,N_loc); //top ghost
      u_new[id] = u[id_ghost];
      v_new[id] = v[id_ghost];
      }
      }

      if(rank_col == 0){
      for(i=1;i<N_loc+1;i++){
      id = ID_2D(1,i,N_loc); //left boundary
      id_ghost = ID_2D(0,i,N_loc); //left ghost
      u_new[id] = u[id_ghost];    
      v_new[id] = v[id_ghost];
      }
      }

      if(rank_col==q-1){
      for(i=1;i<N_loc+1;i++){
      id = ID_2D(N_loc,i,N_loc); //right boundary
      id_ghost = ID_2D(N_loc+1,i,N_loc); //right ghost
      u_new[id] = u[id_ghost];
      v_new[id] = v[id_ghost];
      }
      }
    
      // *** COMMUNICATE

    // *** pack data to send buffers from ghosts
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(0,i,N_loc);
      send_buffer_lft[i] = u[id];
      send_buffer_lft_v[i] = v[id];
    }

    for(i=1;i<N_loc+1;i++){
      id = ID_2D(N_loc+1,i,N_loc);
      send_buffer_rgt[i] = u[id];
      send_buffer_rgt_v[i] = v[id];
    }

    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,0,N_loc);
      send_buffer_bot[i] = u[id];
      send_buffer_bot_v[i] = v[id];
    }

    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,N_loc+1,N_loc);
      send_buffer_top[i] = u[id];
      send_buffer_top_v[i] = v[id];
    }
    
    
// IMPLEMENTING SENDING AND RECIEVING //

    int rank_lft = irank-1;
    int rank_rgt = irank+1;
    int rank_top = irank+q;
    int rank_bot = irank-q;
    
     //Send data to the left
    if(rank_col==0){
      MPI_Recv(recv_buffer_rgt,N_loc,MPI_DOUBLE,rank_rgt,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(recv_buffer_rgt_v,N_loc,MPI_DOUBLE,rank_rgt,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    else if(rank_col==q-1){
      MPI_Send(send_buffer_lft,N_loc,MPI_DOUBLE,rank_lft,101,MPI_COMM_WORLD);
       MPI_Send(send_buffer_lft_v,N_loc,MPI_DOUBLE,rank_lft,101,MPI_COMM_WORLD);
    }
    else{
      MPI_Sendrecv(send_buffer_lft,N_loc,MPI_DOUBLE,rank_lft,101,recv_buffer_rgt,N_loc,MPI_DOUBLE,rank_rgt,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Sendrecv(send_buffer_lft_v,N_loc,MPI_DOUBLE,rank_lft,101,recv_buffer_rgt_v,N_loc,MPI_DOUBLE,rank_rgt,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    //Send data to the right
    if(rank_col==q-1){
      MPI_Recv(recv_buffer_lft,N_loc,MPI_DOUBLE,rank_lft,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(recv_buffer_lft_v,N_loc,MPI_DOUBLE,rank_lft,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    else if(rank_col==0){
      MPI_Send(send_buffer_rgt,N_loc,MPI_DOUBLE,rank_rgt,102,MPI_COMM_WORLD);
       MPI_Send(send_buffer_rgt_v,N_loc,MPI_DOUBLE,rank_rgt,102,MPI_COMM_WORLD);
    }
    else{
      MPI_Sendrecv(send_buffer_rgt,N_loc,MPI_DOUBLE,rank_rgt,102,recv_buffer_lft,N_loc,MPI_DOUBLE,rank_lft,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Sendrecv(send_buffer_rgt_v,N_loc,MPI_DOUBLE,rank_rgt,102,recv_buffer_lft_v,N_loc,MPI_DOUBLE,rank_lft,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    //Send data to the bottom
    if(rank_row==0){
      MPI_Recv(recv_buffer_top,N_loc,MPI_DOUBLE,rank_top,103,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
       MPI_Recv(recv_buffer_top_v,N_loc,MPI_DOUBLE,rank_top,103,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    else if(rank_row==q-1){
      MPI_Send(send_buffer_bot,N_loc,MPI_DOUBLE,rank_bot,103,MPI_COMM_WORLD);
      MPI_Send(send_buffer_bot_v,N_loc,MPI_DOUBLE,rank_bot,103,MPI_COMM_WORLD);
    }
    else{
      MPI_Sendrecv(send_buffer_bot,N_loc,MPI_DOUBLE,rank_bot,103,recv_buffer_top,N_loc,MPI_DOUBLE,rank_top,103,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Sendrecv(send_buffer_bot_v,N_loc,MPI_DOUBLE,rank_bot,103,recv_buffer_top_v,N_loc,MPI_DOUBLE,rank_top,103,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }

    //Send data to the top
    if(rank_row==q-1){
      MPI_Recv(recv_buffer_bot,N_loc,MPI_DOUBLE,rank_bot,104,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(recv_buffer_bot_v,N_loc,MPI_DOUBLE,rank_bot,104,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    else if(rank_row==0){
      MPI_Send(send_buffer_top,N_loc,MPI_DOUBLE,rank_top,104,MPI_COMM_WORLD);
      MPI_Send(send_buffer_top_v,N_loc,MPI_DOUBLE,rank_top,104,MPI_COMM_WORLD);
    }
    else{
      MPI_Sendrecv(send_buffer_top,N_loc,MPI_DOUBLE,rank_top,104,recv_buffer_bot,N_loc,MPI_DOUBLE,rank_bot,104,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Sendrecv(send_buffer_top_v,N_loc,MPI_DOUBLE,rank_top,104,recv_buffer_bot_v,N_loc,MPI_DOUBLE,rank_bot,104,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }   

// END SENDING AND RECIEVING //

    
    // *** unpack data from recv buffers to ghosts
    if(rank_col>0){
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(0,i,N_loc); //left ghost
	u[id_ghost] = recv_buffer_lft[i]; //load receive buffer to left ghost
  v[id_ghost] = recv_buffer_lft_v[i];
      }      
    }

    if(rank_col<q-1){
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(N_loc+1,i,N_loc); //right ghost
	u[id_ghost] = recv_buffer_rgt[i]; //load receive buffer to right ghost
  v[id_ghost] = recv_buffer_rgt_v[i];
      }      
    }

    if(rank_row>0){
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(i,0,N_loc); //bottom ghost
	u[id_ghost] = recv_buffer_bot[i]; //load receive buffer to bottom ghost
  v[id_ghost] = recv_buffer_bot_v[i];
      }      
    }
    
    if(rank_row<q-1){
      for(i=1;i<N_loc+1;i++){
	id_ghost = ID_2D(i,N_loc+1,N_loc); //top ghost
	u[id_ghost] = recv_buffer_top[i]; //load receive buffer to top ghost
  v[id_ghost] = recv_buffer_top_v[i];
      }      
    }

    // COMMUNICATION COMPLETE
    
    
   if(rank_col>0){ //compute left edge
      for(j=1;j<N_loc+1;j++){
	i=1;
	
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top
		u_new[id] = u[id] +  (dt/pow(h,2))*(f[id] + u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);	
	v_new[id] = v[id] +  (dt/pow(h,2))*g[id];
      }
    }
 
  if(rank_col<q-1){ //compute right edge
      for(j=1;j<N_loc+1;j++){
	i=N_loc;
	
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top
		u_new[id] = u[id] +  (dt/pow(h,2))*(f[id] + u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);
	v_new[id] = v[id] +  (dt/pow(h,2))*g[id];
      }
    }
    
   if(rank_row>0){ //compute bottom edge
      for(i=1;i<N_loc+1;i++){
	j=1;
	
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top
		u_new[id] = u[id] +  (dt/pow(h,2))*(f[id] + u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);
	v_new[id] = v[id] +  (dt/pow(h,2))*g[id];
      }
    }
    
   if(rank_row<q-1){ //compute top edge
      for(i=1;i<N_loc+1;i++){
	j=N_loc;
 
	id = ID_2D(i,j,N_loc);
	id_lft = ID_2D(i-1,j,N_loc); //index left
	id_rgt = ID_2D(i+1,j,N_loc); //index righ
	id_bot = ID_2D(i,j-1,N_loc); //index bottom
	id_top = ID_2D(i,j+1,N_loc); //index top
		u_new[id] = u[id] +  (dt/pow(h,2))*(f[id] + u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);
	v_new[id] = v[id] +  (dt/pow(h,2))*g[id];
      }
    }
  
   //update solution 
    for(j=1;j<N+1;j++){
      for(i=1;i<N+1;i++){
	id = ID_2D(i,j,N);
	u[id] = u_new[id];
	v[id] = v_new[id];
      }
    }
    for(int j=1;j<N+1;j++){
    for(int i=1;i<N+1;i++){
      id = ID_2D(i,j,N);
      f[id] = (1.0/eps)*(u[id])*(1.0 - u[id])*(u[id] - ((v[id] + b)/a));
      g[id] = u[id] - v[id];
    }
  }
    }
  
  write_results(u, N, N_loc, nproc, irank, h, L);
  
    
#ifdef DEBUG // print the content of array u including ghosts
  for(j=0; j<N+2; j++){
    for(i=0; i<N+2; i++){
      id = ID_2D(i,j,N);
      printf("%f\t",u[id]);
    }
    printf("\n");
  }
#endif
  
  free(u);
  free(v);
  free(u_new);
  free(v_new);
  free(g);
  free(f);
  free(x);
  free(y);
 
  return 0;
}

// *** Function definitions ***
void getArgs(int *N, int argc, char *argv[], int irank, MPI_Comm comm){

  if(irank==0){
    if ( argc != 2 ) /* argc should be 2 for correct execution */
      {
	//If not correct number of arguments, assume n=1000
	printf("Incorrect number of arguments. Usage: ./poisson N \n");
	MPI_Finalize ( );
	exit ( 1 ); // exit with error code 1
      }
    else
      {
	//Use input argument
	*N = atoi(argv[1]);
      }
  }
  MPI_Bcast(N,1,MPI_INT,0,comm);
}

void write_results(double *u, int N, int N_loc, int nproc, int irank, double h, double L)
 {

  int i,j, id;

  double *u_local       = malloc((N_loc)*(N_loc)*sizeof(double));
  double *u_global      = malloc((N)*(N)*sizeof(double));
  double *u_write      = malloc((N)*(N)*sizeof(double));

  double x,y;
  
  //pack the data for gather (to avoid sending ghosts)
  int id_loc = 0;
  for(j=1;j<N_loc+1;j++){
    for(i=1;i<N_loc+1;i++){
      id = ID_2D(i,j,N_loc);
      u_local[id_loc] = u[id];
      id_loc++;
    }
  }
  
  MPI_Gather(u_local,id_loc,MPI_DOUBLE,u_global,id_loc,MPI_DOUBLE,0,MPI_COMM_WORLD);

  //unpack data so that it is in nice array format

  int id_write, id_global;
  int p_row, p_col;
  int q = sqrt(nproc);

  if(irank==0){
    for(int p=0; p<nproc;p++){
      p_row = p/q;
      p_col = p%q;
      for(j=0;j<N_loc;j++){
	for(i=0;i<N_loc;i++){
	  id_global = p*N_loc*N_loc + j*N_loc + i;
	  id_write  = p_row*N_loc*N_loc*q + j*N_loc*q + p_col*N_loc + i;
	  u_write[id_write] = u_global[id_global];
	}
      }
    }

#ifdef DEBUG
    //writing to screen
    printf("Checking unpacked array again\n");
    for(j=0;j<N;j++){
      for(i=0;i<N;i++){
	id = j*N+i;
	printf("%f ",u_write[id]);
      }
      printf("\n");
    }
#endif

    
    //write to file
    FILE *f;
    f = fopen("results_task2.dat","w"); //open file
    fprintf(f,"x, y, u\n");

    for(j=0; j<N; j++){
      for(i=0; i<N; i++){
	id = j*N + i;
	x = i*h - L; //I am a bit lazy here with not gathering x and y
	y = j*h - L;
	fprintf(f,"%f, %f, %f\n",x,y,u_write[id]);
      }
    }
    fclose(f);
  }
  free(u_local); free(u_global); free(u_write);
}