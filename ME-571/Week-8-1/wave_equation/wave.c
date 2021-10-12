#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>


//function headers (definitions at the end of the file)

// getArgs_mpi reads input parameters npts, dt, time_final 
// from command line
void getArgs(int *npts, double *dt, double *time_final, int argc, char *argv[]);

// computes the parameter alpha and checks whether user provided stable dt
double get_alpha( double c, double dt, double dx);

//reference solutions
double exact_solution(double x, double t, double c);
double exact_solution_derivative(double x, double t, double c);

//utility functions:
double maxValue(double myArray[], int size); 
double minValue(double myArray[], int size); 
double L2_error(int npts, double dx, double* data, double* reference);
void print_results(double *u, double *x, int npts, double time, double c);


/********************************************************************************
 WAVE solves a wave equation u_tt = c^2 u_xx in 1D 
 using 2nd order finite difference formula for both space and time derivative
 *******************************************************************************/

int main(int argc, char * argv[]){

  double *u0; // u^{n-1} solution at previous time-step
  double *u1; // u^n     solution at current time-step
  double *u2; // u^{n+1} solution at next time-step
  double *x;  //         spatial locations of points

  double *u_init; //     initial solution
  double *u_t_init; //   derivative of initial solution

  double dt;  //         time-step 
  double time_final; //  final time for the simulation
  double time; //        current time of the simulation
  double dx;  //         spatial resolution between points

  double c = 1.0; //     wave speed (set to 1.0 by default)
  double alpha; //       CFL non-dimensional parameter parameter alpha = c*dt/dx 

  int npts; // global and local number of points

  int i; 

  //******
  // INITIALIZE SIMULATION
  //******

  //get command-line arguments
  getArgs(&npts, &dt, &time_final, argc,argv);

  //set dx based on the number of points
  dx = 1.0/(npts-1); 

  // compute parameter alpha = c*dt/dx
  alpha = get_alpha(c,dt,dx);
  
  //set-up coordinates and reference solutions
  x = malloc(sizeof(double)*npts);
  u0 = malloc(sizeof(double)*npts);
  u1 = malloc(sizeof(double)*npts);
  u2 = malloc(sizeof(double)*npts);

  u_init = malloc(sizeof(double)*npts);
  u_t_init = malloc(sizeof(double)*npts);

  //******  
  // INITIAL CONDITIONS
  //******

  for(i=0;i<npts;i++){
    // create x coordinates
    x[i] = i*dx;
    
    //use exact solution as initial condition for u 
    u_init[i] = exact_solution(x[i],0,c); 

    //store initial time deriative of initial condition
    u_t_init[i] = exact_solution_derivative(x[i],0,c); 

    //set initial condition
    u1[i] = u_init[i]; 

    //set initial condition for u0 based on initial dudt
    u0[i] = u1[i] - dt*u_t_init[i]; 
  }  

#ifdef DEBUG
  printf(" xmin = %f, xmax = %f\n",minValue(x,npts),maxValue(x,npts));
#endif  

    
    /**********
     BEGIN SIMULATION 
    **********/

    //initialize time variable
    time = 0;

    //compute alpha^2
    double alpha2 = alpha*alpha;

    double ghost_right;
    double ghost_left;
    
    //begin time loop
    while(time<time_final){
      
      //compute interior points
      for(i=1;i<npts-1;i++)
	u2[i] = alpha2*(u1[i-1] - 2*u1[i] + u1[i+1]) + 2*u1[i] - u0[i];

      
      //load left ghost
      ghost_left = u1[npts-1];
      
      //compute left end points
      i=0;
      u2[i] = alpha2*(ghost_left - 2*u1[i] + u1[i+1]) + 2*u1[i] - u0[i];    
      
      //load right ghost
      ghost_right = u1[0];
      
      //compute right end points
      i = npts-1;
      u2[i] = alpha2*(u1[i-1] - 2*u1[i] + ghost_right) + 2*u1[i] - u0[i];	
      
      //increment time
      time = time+dt;
      
      //update u0,u1
      for(i=0;i<npts;i++){
	u0[i] = u1[i];
	u1[i] = u2[i];
      }
      
      
    }//end time loop
   
  /**************
    END SIMULATION 
  **************/

  // Computation of the L2 error, assuming we have made a full revolution 
  
  double u_error = L2_error(npts, dx, u1, u_init);
  printf("%d\t%e\n",npts,u_error);

  //write output data to file
  print_results(u1,x,npts,time,c);

  
  //Deallocation
  free(u0);
  free(u1);
  free(u2);
  free(x);


  return 0;
} 

void getArgs(int *npts, double *dt, double *time_final, int argc, char *argv[])
{

    if ( argc != 4 ) /* argc should be 2 for correct execution */
      {
	//If not correct number of arguments, assume n=1000
	printf("Incorrect number of arguments. Usage: ./wave npts dt time_final \n");
       
       }
    else
      {
	//Use input argument
	*npts = atoi(argv[1]);	
	*dt = atof(argv[2]);
	*time_final = atof(argv[3]);
      }
}

double get_alpha( double c, double dt, double dx)
{
  double alpha = c*dt/dx;

    if ( 1.0 <= fabs ( alpha ) )
    {
      int irank;
      MPI_Comm_rank(MPI_COMM_WORLD,&irank);
      if ( irank == 0 )
	{
	  printf ("\n" );
	  printf ("  CFL condition not met!");
	  printf ("  c = %g\n", c );
	  printf ("  dt = %g\n", dt );
	  printf ("  dx = %g\n", dx );
	  printf ("  alpha = d*dt/dx %g - abs value exceeds 1\n", alpha );
	  printf ("  Computation will not be stable!\n" );
	}
      MPI_Finalize ( );
      exit ( 1 );
    }

    return alpha;

}

double exact_solution(double x, double t, double c)
{
  double pi = acos(-1.0);
  
  double u = sin(2.0 * pi*(x - c*t));

  return u;
}

double exact_solution_derivative(double x, double t, double c)
{
  double pi = acos(-1.0);

  double dudt = -2.0*pi*c*cos (2.0*pi*(x-c*t));

  return dudt;

}

double maxValue(double myArray[], int size) {

  double maxV = myArray[0];

  for (int i = 1; i < size; ++i) {
    if ( myArray[i] > maxV ) {
      maxV = myArray[i];
    }
  }
  return maxV;
}

double minValue(double myArray[], int size) {

  double minV = myArray[0];

  for (int i = 1; i < size; ++i) {
    if ( myArray[i] < minV ) {
      minV = myArray[i];
    }
  }
  return minV;
}

double L2_error(int npts, double dx, double* data, double* reference){
    double error=0.0; 

    for(int j=0;j<npts;j++){
      error += dx*pow((data[j]-reference[j]),2);
    }

    return sqrt(error);
}

void print_results(double *u, double *x, int npts, double time, double c)
  {

    int i;

    FILE *f;
    f = fopen("results_wave_serial.dat","a"); //open file
    fprintf(f,"x\t u\t u_exact\n");
    for(i=0; i<npts; i++){
      fprintf(f,"%f\t%f\t%f\n",x[i],u[i],exact_solution(x[i],time,c));
    }
    fclose(f);


  }
