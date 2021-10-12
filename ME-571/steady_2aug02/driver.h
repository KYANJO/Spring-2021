#ifndef _DRIVER_
#define _DRIVER_

/* Set the precision of Real variables (float or double)   */

#define PRECISION double
typedef PRECISION Real;

#define FIFE_SCALING 0    /* if 1 use fife scaling */
#define GRAPHICS 0        /* if 1 use compile in ezplot graphics */

/* 
 * Memory allocation macros 
 * ------------------------ */

#define  VECTOR(n)  (Real *)malloc((unsigned)(n)*sizeof(Real))
#define DVECTOR(n)  (double *)malloc((unsigned)(n)*sizeof(double))
#define FVECTOR(n)  (float *)malloc((unsigned)(n)*sizeof(float))
#define IVECTOR(n)  (int *) malloc((unsigned)(n)*sizeof(int))

/* 
 * I always define min, max, and make sure M_PI is defined 
 * ------------------------------------------------------- */

#define min(a,b)      ((a)>(b) ? (b) : (a))
#define max(a,b)      ((a)>(b) ? (a) : (b))
#ifndef M_PI
#define M_PI	      3.14159265358979323846
#endif

#define INDEX(j,k) ( (k) + (j)*(tinc) )

extern int nr, ntheta, tinc, ntot;
extern int r_fix, theta_fix;
extern double one_o_eps, one_o_a, b_o_a, tau2;


/* 
 * Prototypes for public functions defined in driver.c 
 * --------------------------------------------------- */
double R (int j);
double Th (int k);

/* Numerical Recipes routines with double precision data */
void realft (double *data, int n, int isign);
/* Lapack routines */
double dnrm2_ (int *ntot, double *rhs, int *one);
void dgbtrf_ (int *m_rows, int *n_cols, int *kl, int *ku, double *AB,
	      int *ldab, int *ipiv, int *info);
void dgbtrs_ (char *trans, int *n_cols, int *kl, int *ku, int *nrhs,
	      double *AB, int *ldab, int *ipiv, double *bvec,
	      int *ldb, int *info);

#endif /* ifndef _DRIVER_ */
