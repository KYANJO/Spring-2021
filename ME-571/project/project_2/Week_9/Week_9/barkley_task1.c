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

#define ID_2D(i,j,nx) ((j)*(nx+2)+(i))

// *** Function interfaces ***

// getArgs_mpi reads input parameters Nfrom command line
void getArgs(int *N, int argc, char *argv[]);

//functions f and g
double* f(int N, double* u, double* v);

double* g(int N, double* u, double* v);

// write results to file
void write_results(double *u, double *x, double *y, int N);


// *** Main function ***
int main(int argc, char * argv[]){

  int N, i, j, id, id_ghost;
  int id_lft, id_rgt, id_top, id_bot;

  //read command line arguments
  getArgs(&N, argc, argv);

  // define domain size (assume both x and y directions are the same)
  double L = 20.0;
  // double L = 150.0;

  //final time
  int tfinal = 40;
  
  //compute the grid spacing
  double h = 2*L/(N-1);
  
  //allocate arrays
  double *u       = malloc((N+2)*(N+2)*sizeof(double));
  double *v       = malloc((N+2)*(N+2)*sizeof(double));
  double *ff      = malloc((N+2)*(N+2)*sizeof(double));
  double *gg      = malloc((N+2)*(N+2)*sizeof(double));
  double *u_new   = malloc((N+2)*(N+2)*sizeof(double));
  double *v_new   = malloc((N+2)*(N+2)*sizeof(double));
  double *x       = malloc((N+2)*sizeof(double));
  double *y       = malloc((N+2)*sizeof(double));

  //initialize x and y
  for(i=1;i<N+1;i++){
    x[i] = -L + (i-1)*h;
    y[i] = -L + (i-1)*h;//y coordinates are the same as x, but need not be
  }
  
  //initialize array u
  for(j=1;j<N+1;j++){
      for(i=1;i<N+1;i++){
	id = ID_2D(i,j,N);
	
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

  ff = f(N, u, v);
  gg = g(N, u,v);


  double dt = 0.001;
  int ntime = (int)tfinal/dt;
  for(int time=1;time<=ntime;time++){
    

    //initialize BC in ghost cells                                                                                                                                                  
    for(i=1;i<N+1;i++){
      id = ID_2D(i,0,N); //bottom ghost                                                                                                                                             
      int id_bn =ID_2D(i,1,N);
      u[id] = u[id_bn];

      id = ID_2D(i,N+1,N); //top ghost                                                                                                                                              
      int id_tn =ID_2D(i,N,N);
      u[id] = u[id_tn];

      id = ID_2D(0,i,N); //left ghost                                                                                                                                               
      int id_ln = ID_2D(1,i,N);
      u[id] = u[id_ln];

      id = ID_2D(N+1,i,N); //right ghost                                                                                                                                            
      int id_rn =ID_2D(N,i,N);
      u[id] = u[id_rn];
    }

    //go over interior points 2:N-1
    for(j=1;j<N+1;j++){
      for(i=1;i<N+1;i++){
        id = ID_2D(i,j,N);
        id_lft = ID_2D(i-1,j,N); //index left
        id_rgt = ID_2D(i+1,j,N); //index righ
        id_bot = ID_2D(i,j-1,N); //index bottom
        id_top = ID_2D(i,j+1,N); //index top

        u_new[id] = u[id] + dt*ff[id] + (dt/pow(h,2))*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - 4*u[id]);
        
        
        v_new[id] = v[id] +  dt*gg[id];
      }
    }
    
    //enforce Dirichlet BC by copying values from ghosts
    for(i=1;i<N+1;i++){
      id = ID_2D(i,1,N); //bottom boundary
      id_ghost = ID_2D(i,0,N); //bottom ghost
      u_new[id] = u[id_ghost];
      v_new[id] = v[id_ghost];

      id = ID_2D(i,N,N); //top boundary
      id_ghost = ID_2D(i,N+1,N); //top ghost
      u_new[id] = u[id_ghost];
      v_new[id] = v[id_ghost];

      id = ID_2D(1,i,N); //left boundary
      id_ghost = ID_2D(0,i,N); //left ghost
      u_new[id] = u[id_ghost];    
      v_new[id] = v[id_ghost];

      id = ID_2D(N,i,N); //right boundary
      id_ghost = ID_2D(N+1,i,N); //right ghost
      u_new[id] = u[id_ghost];
      v_new[id] = v[id_ghost];
    }
    
    //update solution 
    for(j=1;j<N+1;j++){
      for(i=1;i<N+1;i++){
        id = ID_2D(i,j,N);
        u[id] = u_new[id];
        v[id] = v_new[id];
      }
    }
    ff = f(N,u,v);
    gg = g(N,u,v);
  }
  
  write_results(u, x, y, N);
  
    
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
  free(gg);
  free(ff);
  free(x);
  free(y);
 
  return 0;
}


// *** Function definitions ***
void getArgs(int *N, int argc, char *argv[])
{

    if ( argc != 2 ) /* argc should be 2 for correct execution */
      {
	//If not correct number of arguments, assume n=1000
	printf("Incorrect number of arguments. Usage: ./poisson N \n");
       
       }
    else
      {
	//Use input argument
	*N = atoi(argv[1]);	
      }
}

double* f(int N, double* u, double* v){
  
  double *ff = malloc((N+2)*(N+2)*sizeof(double));
  double epslon = 0.02, a =0.75, b = 0.01;
  //double epslon = 0.02, a =0.75, b = 0.02;
  int id;

  for(int j=1;j<N+1;j++){
    for(int i=1;i<N+1;i++){
      id = ID_2D(i,j,N);
      ff[id] = (1.0/epslon)*(u[id])*(1.0 - u[id])*(u[id] - ((v[id] + b)/a));

    }
  }
  return ff;
}

double* g(int N, double* u, double* v){

  double *gg = malloc((N+2)*(N+2)*sizeof(double));
  int id;

  for(int j=1;j<N+1;j++){
    for(int i=1;i<N+1;i++){
      id = ID_2D(i,j,N);
      gg[id] = u[id] - v[id];
	
    }
  }
  return gg;
}

void write_results(double *u, double *x, double *y, int N){
  int i,j, id;

  FILE *f;
  f = fopen("results_barkley_serial.dat","a"); //open file
  fprintf(f,"x, y, u\n");
  
  for(j=1; j<N+1; j++){
    for(i=1; i<N+1; i++){
      id = ID_2D(i,j,N);
      fprintf(f,"%f, %f, %f\n",x[i],y[j],u[id]);
    }
  }
  fclose(f);
  
}
