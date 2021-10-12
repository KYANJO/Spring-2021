#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define ID_2D(i,j,nx) ((j)*(nx+2)+(i))

void initialData(double *ip, const double ival, int size);
void matmullocal(double *A, double *B, double *C, const int N);
void getArgs(int *N, int argc, char *argv[], int irank, MPI_Comm comm);
void checkResult(double *hostRef, double *gpuRef, const int N);
void printMatrix(double *C, const int nx, const int ny);

// *** Main function ***
int main(int argc, char * argv[])
{

    int n,n_loc,i,j,k,aid,bid,id,ida,idb;
    int irank, nproc;

    int root_rank = 0;
    
    MPI_Init(&argc,&argv); //initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&irank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

     //compute row and column index in the rank grid                        
    int q = (int)sqrt(nproc); //assume nproc is a square number 
    int rank_row = irank/q;
    int rank_col = irank%q;

     //read command line arguments
    getArgs(&n, argc, argv, irank, MPI_COMM_WORLD);

    //compute the local problem size
    n_loc = n/q;

    int nxy = n_loc * n_loc;
    int nBytes = nxy * sizeof(double);

    // malloc host memory
    double *A, *B, *C, *C_serial;
    A = (double *)malloc(nBytes);
    B = (double *)malloc(nBytes);
    C = (double *)malloc(nBytes);
    C_serial = (double *)malloc(nBytes);


    // initialize matrices A and B
    initialData(A,2.0f,nxy);
    initialData(B,0.5f, nxy);

    //memset(C, 0, nBytes);
    //memset(C_serial, 0, nBytes);

    //Comparission with the serial 
    matmullocal(A,B,C_serial,n);
     if (irank == root_rank){
      printMatrix(C_serial, n, n);
     }

    int rank_lft = irank - 1;
    int rank_rgt = irank + 1;
    int rank_top = irank + q;
    int rank_bot = irank - q;
    //initial alignment of sub-matrices
    for(i=0; i<=q-1; i++)
    {
      for(j=0; j<=q-1; j++)
      {
        id = ID_2D(i,j,q);
        aid = ID_2D(i,(j-1+q)%q,q);
        bid = ID_2D((i-j+q)%q,j,q);
      
        MPI_Sendrecv(&A[id],1,MPI_DOUBLE,rank_lft,101,&A[aid],1,MPI_DOUBLE,rank_rgt,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Sendrecv(&B[id],1,MPI_DOUBLE,rank_top,102,&B[bid],1,MPI_DOUBLE,rank_bot,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
    } 
    
    //int rank_rgt = irank + 1;
    //int rank_bot = irank - q;
    /* double *A_rec, *B_rec;
    A_rec = (double *)malloc(nBytes);
    B_rec = (double *)malloc(nBytes);
    
    if (rank_col != q-1)
    {
      MPI_Recv(A_rec,nxy,MPI_DOUBLE,rank_rgt,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    if (rank_row != 0)
    {
      MPI_Recv(B_rec,nxy,MPI_DOUBLE,rank_bot,102,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    } */
    //Each rank multiplies its local A and B sub-matrix, and increments its local
   matmullocal(A,B,C,n_loc);

     for(k=1; k<=q-1; k++)
    {
      //for(i = 0,j=0; i<=q-1,j<=q-1; i++,j++)
      for(i=0; i<=q-1; i++)
      {
        for(j=0; j<=q-1; j++)
        {
          id = ID_2D(i,j,q);
          ida = ID_2D(i,j-1,q);
          idb = ID_2D(i-1,j,q);
          A[id] = A[ida];
          B[id] = B[idb];
        }
      }
      matmullocal(A,B,C,n_loc);
    }

    if (irank == root_rank){
      printMatrix(C, n, n);

      //check results
      //checkResult(C_serial,C,n);
    } 
    
    free(A);
    free(B);
    //free(A_rec);
    //free(B_rec);
    free(C);
    free(C_serial);

    MPI_Finalize(); 
    return 0;
}


//Generate random matrices
void initialData(double *ip, const double ival, int size)
 {  
    int i;
     for (i = 0; i < size; i++)
     {
         ip[i] = (double)(rand() & 0xFF) / 100.0f;
     }
 
     return;
 }

//local Matrix product
 void matmullocal(double *A, double *B, double *C, const int N)
{
  int   id, ida, idb;
  double cc;
  int iy,ix,k;
    for (iy = 0; iy < N; iy++)
    {
        for (ix = 0; ix < N; ix++)
        {
            
            cc = 0;
            for (k = 0; k < N; k++){
                ida = iy*N + k;
                idb = k *N + ix;
                    cc += A[ida]*B[idb];
	        }
	  id = iy*N+ix;
	  C[id] = cc;
        }
    }

    return;
}

void printMatrix(double *C, const int nx, const int ny)
{
    double *ic = C;

    for (int iy = 0; iy < ny; iy++)
    {
        for (int ix = 0; ix < nx; ix++)
        {
            printf("%f ", ic[ix]);

        }

        ic += nx;
        printf("\n");
    }

    return;
}

void checkResult(double *hostRef, double *gpuRef, const int N)
{
    double epsilon = 1.0E-6;
    _Bool match = 1;
    for (int i = 0; i < N; i++)
    {
      if (abs(hostRef[i] - gpuRef[i])/abs(hostRef[i]) > epsilon)
        {
            match = 0;
            printf("host %f gpu %f, err = %e\n", hostRef[i], gpuRef[i], abs(hostRef[i]-gpuRef[i])/abs(hostRef[i]));
            break;
        }
    }
    if (match)
        printf("Arrays match.\n\n");
    else
        printf("Arrays do not match.\n\n");
}



//get inline arguments
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

