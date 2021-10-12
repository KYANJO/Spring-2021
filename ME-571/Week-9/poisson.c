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

//compute integral L2 error
double L2_error(int npts, double dx, double L, double* data, double* reference);

// write results to file
void write_results(double *u, double *u_exact, double *x, double *y, int N);



// *** Main function ***
int main(int argc, char * argv[]){

  int N, i, j, id, id_ghost;
  int id_lft, id_rgt, id_top, id_bot;

  double eps = 1e-8; //error tolerance for Jacobi method
  
  //read command line arguments
  getArgs(&N, argc, argv);

  // define domain size (assume both x and y directions are the same)
  double L = 1.0;
  
  //compute the grid spacing
  double h = L/(N-1);
  
  //allocate arrays
  double *u       = malloc((N+2)*(N+2)*sizeof(double));
  double *u_new   = malloc((N+2)*(N+2)*sizeof(double));
  double *u_exact = malloc((N+2)*(N+2)*sizeof(double));
  double *f       = malloc((N+2)*(N+2)*sizeof(double));
  double *x       = malloc((N+2)*sizeof(double));
  double *y       = malloc((N+2)*sizeof(double));

  //initialize x and y
  for(i=1;i<N+1;i++){
    x[i] = (i-1)*h;
    y[i] = (i-1)*h;//y coordinates are the same as x, but need not be
  }
  
  //initialize array u
  for(j=1;j<N+1;j++){
      for(i=1;i<N+1;i++){
	id = ID_2D(i,j,N);

	u_exact[id] = sin(2*M_PI*x[i])*sin(2*M_PI*y[j]);

	f[id] = -8*M_PI*M_PI*sin(2*M_PI*x[i])*sin(2*M_PI*y[j]);

	u[id] = 0.0;
      }
  }

  //initialize BC in ghost cells
  for(i=1;i<N+1;i++){
    id = ID_2D(i,0,N); //bottom ghost
    u[id] = 0.0;
    
    id = ID_2D(i,N+1,N); //top ghost
    u[id] = 0.0;

    id = ID_2D(0,i,N); //left ghost
    u[id] = 0.0;

    id = ID_2D(N+1,i,N); //right ghost
    u[id] = 0.0;
  }

  
  //compute initial error
  double error = L2_error(N, h, L, u, u_exact);

  int iter=0;
  //begin Jacobi iteration
  while(error>eps){
    
    //go over interior points 2:N-1
    for(j=2;j<N;j++){
      for(i=2;i<N;i++){
	id = ID_2D(i,j,N);
	id_lft = ID_2D(i-1,j,N); //index left
	id_rgt = ID_2D(i+1,j,N); //index righ
	id_bot = ID_2D(i,j-1,N); //index bottom
	id_top = ID_2D(i,j+1,N); //index top

	u_new[id] = 0.25*(u[id_lft] + u[id_rgt] + u[id_bot] + u[id_top] - h*h*f[id]);
      }
    }
    
    //enforce Dirichlet BC by copying values from ghosts
    for(i=1;i<N+1;i++){
      id = ID_2D(i,1,N); //bottom boundary
      id_ghost = ID_2D(i,0,N); //bottom ghost
      u_new[id] = u[id_ghost];
      
      id = ID_2D(i,N,N); //top boundary
      id_ghost = ID_2D(i,N+1,N); //top ghost
      u_new[id] = u[id_ghost];
      
      id = ID_2D(1,i,N); //left boundary
      id_ghost = ID_2D(0,i,N); //left ghost
      u_new[id] = u[id_ghost];    
      
      id = ID_2D(N,i,N); //right boundary
      id_ghost = ID_2D(N+1,i,N); //right ghost
      u_new[id] = u[id_ghost];
    }
    
    //compute error
    error = L2_error(N, h, L, u_new, u);
    iter++;
    //printf("iter = %d, error = %e\n",iter,error);
    
    //update solution 
    for(j=1;j<N+1;j++){
      for(i=1;i<N+1;i++){
	id = ID_2D(i,j,N);
	u[id] = u_new[id];
      }
    }
        
  }
  
  printf("N = %d, iter = %d, eps = %e, error = %e\n",N,iter,eps,error);

  write_results(u, u_exact, x, y, N);
  
    
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
  free(u_new);
  free(u_exact);
  free(f);
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


double L2_error(int N, double dx, double L, double* data, double* reference){

  double error=0.0; 
    int id;
    
    for(int j=1;j<N+1;j++){
      for(int i=1;i<N+1;i++){
	id = ID_2D(i,j,N);
	error += pow((data[id]-reference[id]),2);
      }
    }

    return sqrt(error)/N;
}

void write_results(double *u, double *u_exact, double *x, double *y, int N){
  int i,j, id;

  FILE *f;
  f = fopen("results_poisson_serial.dat","a"); //open file
  fprintf(f,"x, y, u, u_exact\n");
  
  for(j=1; j<N+1; j++){
    for(i=1; i<N+1; i++){
      id = ID_2D(i,j,N);
      fprintf(f,"%f, %f, %f, %f\n",x[i],y[j],u[id],u_exact[id]);
    }
  }
  fclose(f);
  
}
