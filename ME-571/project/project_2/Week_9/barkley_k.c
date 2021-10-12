/*
The code solves a 2D Poisson equation in a square domain using finite diference method

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

// write results to file
void write_results(double *u, int N, int N_loc, int nproc, int irank, double h, double L);

// *** Main function ***
int main(int argc, char * argv[]){

  int N, i, j, id, id_ghost,M,Tfinal;
  int id_lft, id_rgt, id_top, id_bot;
  
  double dt;
  
  double eps = 0.02, a =0.75, b = 0.01;

  int irank, nproc;

  //initialize MPI                                                       
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&irank);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
 
  //read command line arguments
  getArgs(&N, argc, argv, irank, MPI_COMM_WORLD);

  // define domain size (assume both x and y directions are the same)
  double L = 20.0;

  //compute the grid spacing
  double h = 2*L/(N-1);

  //compute local problem size                                           
  int N_loc = N/nproc;
  
  //allocate arrays
  double *u       = malloc((N+2)*(N_loc+2)*sizeof(double));
  double *v       = malloc((N+2)*(N_loc+2)*sizeof(double));
  double *f      = malloc((N+2)*(N_loc+2)*sizeof(double));
  double *g      = malloc((N+2)*(N_loc+2)*sizeof(double));
  double *u_new   = malloc((N+2)*(N_loc+2)*sizeof(double));
  double *v_new   = malloc((N+2)*(N_loc+2)*sizeof(double));
  double *x       = malloc((N+2)*sizeof(double));
  double *y       = malloc((N+2)*sizeof(double));

  //allocate bufers                                                     
  double *send_bufer_rgt = malloc(N*sizeof(double));
  double *send_bufer_lft = malloc(N*sizeof(double));

  double *recv_bufer_rgt = malloc(N*sizeof(double));
  double *recv_bufer_lft = malloc(N*sizeof(double));

  //initialize x and y
  for(i=1;i<N_loc+1;i++){
    x[i] = -L + (i-1)*h + irank*N_loc*h;
  }
  for(i=1;i<N;i++){
    y[i] = -L + (i-1)*h;
  }
  
  //initialize array u
  for(j=1;j<N+1;j++){
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

  for(int j=1;j<N+1;j++){
      for(int i=1;i<N_loc+1;i++){
          id = ID_2D(i,j,N_loc);
          f[id] = (1/eps)*(u[id])*(1 - u[id])*(u[id] - ((v[id] + b)/a));
          g[id] = u[id] - v[id];
      }
    }

  //start time loop
  double wtime, wtime_global;
  wtime = MPI_Wtime();
  
  //tfinal
  Tfinal = 40; 
  
  //time spacing
  dt = 0.025;
  
  //# time steps
  M = (int)Tfinal/dt;
  
  // Time loop
  for(int k=1;k<=M;k++){
    
    //initialize BC                                                      
      for(i=1;i<N_loc+1;i++){
          id = ID_2D(i,0,N_loc); //bottom ghost                                            
          int id_bn =ID_2D(i,1,N_loc);
          u[id] = u[id_bn];
      }

      for(i=1;i<N_loc+1;i++){
          id = ID_2D(i,N+1,N_loc); //top ghost                                         
          int id_tn =ID_2D(i,N,N_loc);
          u[id] = u[id_tn];
      }

      for(i=1;i<N+1;i++){
          id = ID_2D(0,i,N_loc); //left ghost                                              
          int id_ln = ID_2D(1,i,N_loc);
          u[id] = u[id_ln];
      }

      for(i=1;i<N+1;i++){
          id = ID_2D(N_loc+1,i,N_loc); //right ghost                                       
          int id_rn =ID_2D(N_loc,i,N_loc);
          u[id] = u[id_rn];
      }

    //go over interior points 2:N_loc
    for(j=2;j<N;j++){
      for(i=2;i<N_loc;i++){
          id = ID_2D(i,j,N_loc);
          id_lft = ID_2D(i-1,j,N_loc); //index left
          id_rgt = ID_2D(i+1,j,N_loc); //index righ
          id_bot = ID_2D(i,j-1,N_loc); //index bottom
          id_top = ID_2D(i,j+1,N_loc); //index top

          u_new[id] = u[id] + dt*(f[id])+ (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);
          
          v_new[id] = v[id] +  dt*g[id];
      }
    }
    
    //enforce BC by copying values from ghosts
      for(i=1;i<N_loc+1;i++){
        j = 1;
          
          id = ID_2D(i,j,N_loc);
          id_lft = ID_2D(i-1,j,N_loc); //index left
          id_rgt = ID_2D(i+1,j,N_loc); //index righ
          id_bot = ID_2D(i,j-1,N_loc); //index bottom
          id_top = ID_2D(i,j+1,N_loc); //index top

          u_new[id] = u[id] + dt*(f[id])+ (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);
          
          v_new[id] = v[id] +  dt*g[id];
      }

      for(i=1;i<N_loc+1;i++){
        j = N;

          id = ID_2D(i,j,N_loc);
          id_lft = ID_2D(i-1,j,N_loc); //index left
          id_rgt = ID_2D(i+1,j,N_loc); //index righ
          id_bot = ID_2D(i,j-1,N_loc); //index bottom
          id_top = ID_2D(i,j+1,N_loc); //index top

          u_new[id] = u[id] + dt*(f[id])+ (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);
          
          v_new[id] = v[id] +  dt*g[id];
      }
    
    if(irank == 0){
      for(j=1;j<N+1;j++){
        i = 1;

          id = ID_2D(i,j,N_loc);
          id_lft = ID_2D(i-1,j,N_loc); //index left
          id_rgt = ID_2D(i+1,j,N_loc); //index righ
          id_bot = ID_2D(i,j-1,N_loc); //index bottom
          id_top = ID_2D(i,j+1,N_loc); //index top

          u_new[id] = u[id] + dt*(f[id])+ (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);
          
          v_new[id] = v[id] +  dt*g[id];
      }
    }

    if(irank==nproc-1){
      for(j=1;j<N+1;j++){
        i = N_loc;
         
        id = ID_2D(i,j,N_loc);
        id_lft = ID_2D(i-1,j,N_loc); //index left
        id_rgt = ID_2D(i+1,j,N_loc); //index righ
        id_bot = ID_2D(i,j-1,N_loc); //index bottom
        id_top = ID_2D(i,j+1,N_loc); //index top

        u_new[id] = u[id] + dt*(f[id])+ (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);
        
        v_new[id] = v[id] +  dt*g[id];
      }
    }
                                                 
    if(nproc>1)
    {
    // *** pack data to send bufers from ghosts                       
    for(i=1;i<N+1;i++){
        id = ID_2D(0,i,N_loc);
        send_bufer_lft[i] = u[id];
      }

      for(i=1;i<N+1;i++){
        id = ID_2D(N_loc+1,i,N_loc);
        send_bufer_rgt[i] = u[id];
      }
    
    int rank_lft = irank - 1;
    int rank_rgt = irank + 1;

    //All ranks except those in first column send to the left            
    //move Data to the left                                              
    if(irank==0){
      MPI_Recv(recv_bufer_rgt,N,MPI_DOUBLE,rank_rgt,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    else if(irank==nproc-1){
      MPI_Send(send_bufer_lft,N,MPI_DOUBLE,rank_lft,101,MPI_COMM_WORLD);
      }
    else{
      MPI_Sendrecv(send_bufer_lft,N,MPI_DOUBLE,rank_lft,101,recv_bufer_rgt,N,MPI_DOUBLE,rank_rgt,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }

    //All ranks except those in last column send to the right            
    //Move Data to the right                                             
    if(irank==nproc-1){
      MPI_Recv(recv_bufer_lft,N,MPI_DOUBLE,rank_lft,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    else if(irank==0){
      MPI_Send(send_bufer_rgt,N,MPI_DOUBLE,rank_rgt,102,MPI_COMM_WORLD);
      }
    else{
      MPI_Sendrecv(send_bufer_rgt,N,MPI_DOUBLE,rank_rgt,102,recv_bufer_lft,N,MPI_DOUBLE,rank_lft,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }

    // *** unpack data from recv bufers to ghosts                       
   if(irank>0){
      for(i=1;i<N+1;i++){
        id_ghost = ID_2D(0,i,N_loc); //left ghost                        
        u[id_ghost] = recv_bufer_lft[i]; //load receive bufer to left ghost                                                                    
      }
   }

   if(irank<nproc-1){
      for(i=1;i<N+1;i++){
        id_ghost = ID_2D(N_loc+1,i,N_loc); //right ghost                 
        u[id_ghost] = recv_bufer_rgt[i]; //load receive bufer to right ghost                                                                   
      }
   }

    // COMMUNICATION COMPLETE         
    //complete computation at the rank edges (except the domain boundary conditions)
    if(irank>0){ //compute left edge                                  
      for(j=1;j<N+1;j++){
        i=1;

        id = ID_2D(i,j,N_loc);
        id_lft = ID_2D(i-1,j,N_loc); //index left                        
        id_rgt = ID_2D(i+1,j,N_loc); //index righ                        
        id_bot = ID_2D(i,j-1,N_loc); //index bottom                      
        id_top = ID_2D(i,j+1,N_loc); //index top                         
	
        u_new[id] = u[id] + dt*f[id] + (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);

        v_new[id] = v[id] +  dt*g[id];
      }
    }

  //compute right edge 
  if(irank<nproc-1){                              
      for(j=1;j<N+1;j++){
        i=N_loc;

        id = ID_2D(i,j,N_loc);
        id_lft = ID_2D(i-1,j,N_loc); //index left                        
        id_rgt = ID_2D(i+1,j,N_loc); //index righ                        
        id_bot = ID_2D(i,j-1,N_loc); //index bottom                      
        id_top = ID_2D(i,j+1,N_loc); //index top                         

        u_new[id] = u[id] + dt*f[id] + (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);

        v_new[id] = v[id] +  dt*g[id];
      }
  }

     //compute bottom edge                                
      for(i=1;i<N_loc+1;i++){
	    j=1;

        id = ID_2D(i,j,N_loc);
        id_lft = ID_2D(i-1,j,N_loc); //index left                        
        id_rgt = ID_2D(i+1,j,N_loc); //index righ                        
        id_bot = ID_2D(i,j-1,N_loc); //index bottom                      
        id_top = ID_2D(i,j+1,N_loc); //index top                         
       
        u_new[id] = u[id] + dt*f[id] + (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);

        v_new[id] = v[id] +  dt*g[id];
      }
    
  //compute top edge                                
      for(i=1;i<N_loc+1;i++){
      	j = N;

        id = ID_2D(i,j,N_loc);
        id_lft = ID_2D(i-1,j,N_loc); //index left                        
        id_rgt = ID_2D(i+1,j,N_loc); //index righ                        
        id_bot = ID_2D(i,j-1,N_loc); //index bottom                      
        id_top = ID_2D(i,j+1,N_loc); //index top                         
       
        u_new[id] = u[id] + dt*f[id] + (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);

        v_new[id] = v[id] +  dt*g[id];
      }
    
    }
    // computation done here                                             

    //update solution 
    for(j=1;j<N+1;j++){
      for(i=1;i<N_loc+1;i++){
        id = ID_2D(i,j,N_loc);
        u[id] = u_new[id];
        v[id] = v_new[id];
      }
    }

    for(int j=1;j<N+1;j++){
      for(int i=1;i<N_loc+1;i++){
          id = ID_2D(i,j,N_loc);
          f[id] = (1/eps)*(u[id])*(1 - u[id])*(u[id] - ((v[id] + b)/a));
          g[id] = u[id] - v[id];
      }
    }
  }

  wtime = MPI_Wtime()-wtime;

  MPI_Reduce(&wtime, &wtime_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(irank == 0){
    printf("%d,%d,%f\n",nproc,N,wtime_global);
  }
    write_results(u, N, N_loc, nproc, irank, h, L);
  
  free(u);
  free(v);
  free(u_new);
  free(v_new);
  free(g);
  free(f);
  free(x);
  free(y);

  MPI_Finalize(); 
  return 0;
}


// *** Function definitions ***
void getArgs(int *N, int argc, char *argv[], int irank, MPI_Comm comm)
{

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

void write_results(double *u, int N, int N_loc, int nproc, int irank, double h, double L){

  int i,j, id;

  double *u_local       = malloc((N)*(N_loc)*sizeof(double));
  double *u_global      = malloc((N)*(N)*sizeof(double));
  double *u_write      = malloc((N)*(N)*sizeof(double));

  double x,y;
  
  //pack the data for gather (to avoid sending ghosts)
  int id_loc = 0;
  for(j=1;j<N+1;j++){
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
  int q = nproc;

  if(irank==0){
    
    for(int p=0; p<nproc;p++){
      p_row = p/nproc;
      p_col = p%q;
      for(j=0;j<N;j++){
    	for(i=0;i<N_loc;i++){
    	  id_global = p*N*N_loc + j*N_loc + i;
        id_write  = p_row*N_loc*N*q + j*N_loc*q + p_col*N_loc + i;
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
    f = fopen("results_barkley_task3.dat","w"); //open file
    fprintf(f,"x, y, u\n");

    for(j=0; j<N; j++){
      for(i=0; i<N; i++){
	id = j*N + i;
	x = -L + i*h;
	y = -L + j*h;
	
	fprintf(f,"%f, %f, %f\n",x,y,u_write[id]);
      }
    }
    fclose(f);
  }
  free(u_local); free(u_global); free(u_write);
  
}
