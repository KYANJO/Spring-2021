/* Include Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "driver.h"
#include "operators.h"
#include "solver.h"

#define BANDED 1	/* 1 for banded (un-banded only used for testing) */
#define FAST   1	/* 1 for fast but tricky matrix build */

#if BANDED			/* banded storage */
#define AB_INDEX  (kl + ku + i_f77 - j_f77 + (j_f77 - 1) * ldab)
#else /* standard storage */
#define AB_INDEX  ((i_f77 - 1) + (j_f77 - 1) * ldab)
#endif

/* Matrix elements set in operators.c */
extern double *Ld, *Ll, *Lu;
extern double *omega_column;

/* 
 * Private variables 
 * ----------------- */
static double *AB;
static int info, kl, ku, ldab, n_cols, *ipiv;
static int j_f77_fix, istart_fix, istop_fix;

/* 
 * Private functions 
 * ----------------- */
static int Convert_indx (int i);

/*===========================================================================*/

void
Solve_A (double *u)
{
  int i_f77, ldb, nrhs;
  char trans = 'N';
  double *bvec;
  double coef;

  ldb = n_cols;
  nrhs = 2;
  bvec = VECTOR (nrhs * ldb);

  for (i_f77 = 1; i_f77 <= ldb; i_f77++)
    bvec[i_f77 - 1] = u[Convert_indx (i_f77)];

  /* Since omega is variable we will use Sherman-Morrison */
  for (i_f77 = 1; i_f77 <= ldb; i_f77++) {
    bvec[ldb + i_f77 - 1] = omega_column[Convert_indx (i_f77)];
  }
  for (i_f77 = istart_fix; i_f77 <= istop_fix; i_f77++) {
    bvec[ldb + i_f77 - 1] -= 1.;
  }

  /* Solve using LAPACK */
  dgbtrs_ (&trans, &n_cols, &kl, &ku, &nrhs, AB, &ldab, ipiv,
	   bvec, &ldb, &info);

  if (info != 0) {
    printf ("info = %d\n", info);
  }

  for (i_f77 = 1; i_f77 <= ldb; i_f77++) {
    u[Convert_indx (i_f77)] = bvec[i_f77 - 1];
  }

  /* Here is the Sherman-Morrison correction */
  coef = bvec[j_f77_fix - 1] / (1. + bvec[ldb + j_f77_fix - 1]);
  for (i_f77 = 1; i_f77 <= ldb; i_f77++) {
    u[Convert_indx (i_f77)] -= coef * bvec[ldb + i_f77 - 1];
  }

  B_op (u);

  free (bvec);
}

/*===========================================================================*/

void
Build_A (void)
{
  int m_rows, i_f77, j_f77, istart, istop, block, j, k;
  double *u1;

  kl = ku = ntheta;
  m_rows = n_cols = ntheta * nr + 1;
#if BANDED
  ldab = 2 * kl + ku + 1;
#else
  ldab = m_rows;
#endif

  if (AB == NULL) {
    /* only allocate AB and ipiv once */
    AB = VECTOR (ldab * n_cols);
    if (AB == NULL) {
      printf ("Matrix A memory failure\n");
      exit (1);
    }
    ipiv = IVECTOR (n_cols);
  }
  u1 = VECTOR (ntot);

  /* Zero AB */
  for (j = 0; j < ldab * n_cols; j++) {
    AB[j] = 0.;
  }

#if FAST
  /* Fast but tricky build */

  /* First fill in each r=constant block with terms due to theta
   * derivatives and kinetics (the greens function means that the
   * linearized kinetics fills out each r-block).
   * 
   * The values are found by setting one unit value in each r-block of
   * u1 and calling To_Build () (defined in operators.c.  It computes
   * theta derivatives and linearized kinetics.)  */

  /* r=0 block */

  i_f77 = j_f77 = 1;
  for (j = 0; j < ntot; j++) {
    u1[j] = 0.;
  }
  /* set unit value (all others zeroed above) */
  u1[Convert_indx (j_f77)] = 1.;
  To_Build (u1, u1);
  /* set matrix element */
  AB[AB_INDEX] = u1[Convert_indx (i_f77)];

  /* r>0 blocks */

  /* do all other blocks at same time because they do not couple */
  for (k = 0; k < ntheta; k++) {

    for (j = 0; j < ntot; j++) {
      u1[j] = 0.;
    }

    /* set one unit value in each r-block */
    for (block = 2; block <= nr + 1; block++) {
      j_f77 = k + (block - 2) * ntheta + 2;
      u1[Convert_indx (j_f77)] = 1.;
    }

    To_Build (u1, u1);

    /* set matrix elements */
    for (block = 2; block <= nr + 1; block++) {
      j_f77 = k + (block - 2) * ntheta + 2;
      istart = (block - 2) * ntheta + 2;
      istop = istart + ntheta - 1;
      for (i_f77 = istart; i_f77 <= istop; i_f77++) {
	AB[AB_INDEX] = u1[Convert_indx (i_f77)];
      }
    }
  }

  /* Now add the r-derivative contributions.  We know these directly
   * from the matrix elements Ld, Lu, Ll. */

  /* Diagonal elements */
  j_f77 = i_f77 = 1;
  AB[AB_INDEX] += Ld[0];
  for (j = 1; j <= nr; j++) {
    for (k = 0; k < ntheta; k++) {
      j_f77 = k + (j - 1) * ntheta + 2;
      i_f77 = j_f77;
      AB[AB_INDEX] += Ld[j];
    }
  }
  /* Lower diagonal elements */
  for (j = 1; j <= nr; j++) {
    for (k = 0; k < ntheta; k++) {
      i_f77 = k + (j - 1) * ntheta + 2;
      j_f77 = max (1, i_f77 - ntheta);
      AB[AB_INDEX] += Ll[j - 1];
    }
  }
  /* Upper diagonal elements */
  for (j = 1; j <= nr; j++) {
    for (k = 0; k < ntheta; k++) {
      j_f77 = k + (j - 1) * ntheta + 2;
      i_f77 = max (1, j_f77 - ntheta);
      AB[AB_INDEX] += Lu[j - 1];
      /* note the j-1 in the above line is a bit strange but correct
       * due to the relationship between the fortran and c indexing */
    }
  }

  /* End of fast build code */

#else

  /*  Slow but reliable build */

  for (j_f77 = 1; j_f77 <= n_cols; j_f77++) {
    /* Set u1 = 0 except element j_f77 which is 1 */
    for (j = 0; j < ntot; j++) {
      u1[j] = 0.;
    }
    u1[Convert_indx (j_f77)] = 1.;

    /* call to B_op is necessary for proper treatment of center */
    if (j_f77 == 1) {
      B_op (u1);
    }

    Eval_Au (u1, u1);

    istart = max (1, j_f77 - ku);
    istop = min (m_rows, j_f77 + kl);
    for (i_f77 = istart; i_f77 <= istop; i_f77++) {
      AB[AB_INDEX] = u1[Convert_indx (i_f77)];
    }
  }

  /*  End of slow build */

#endif

  /* Because omega is a variable */
  j_f77 = 2 + ntheta * (r_fix - 1) + theta_fix;
  istart_fix = max (1, j_f77 - ku);
  istop_fix = min (m_rows, j_f77 + kl);
  for (i_f77 = istart_fix; i_f77 <= istop_fix; i_f77++) {
    AB[AB_INDEX] = 1.;
  }
  j_f77_fix = j_f77;		/* save for use in solve */

  /*  printf ("finished making AB \n"); */

#if 0
  /* check that matrix is really banded */
  printf ("Checking that matrix is banded \n");
  for (j_f77 = 1; j_f77 <= n_cols; j_f77++) {
    istart = max (1, j_f77 - ku);
    istop = min (m_rows, j_f77 + kl);
    for (i_f77 = 1; i_f77 < istart; i_f77++) {
      if (u1[Convert_indx (i_f77)] != 0.) {
	printf ("element not zero: i,j,u = %d, %d, %g\n",
		i_f77, j_f77, u1[Convert_indx (i_f77)]);
      }
    }
    for (i_f77 = istop + 1; i_f77 <= n_cols; i_f77++) {
      if (u1[Convert_indx (i_f77)] != 0.) {
	printf ("element not zero: i,j,u = %d, %d, %g\n",
		i_f77, j_f77, u1[Convert_indx (i_f77)]);
      }
    }
  }
#endif

#if 0
  /* output matrix */
  for (i_f77 = 1; i_f77 <= ldab; i_f77++) {
    for (j_f77 = 1; j_f77 <= n_cols; j_f77++) {
      printf ("%-8.3f", AB[(i_f77 - 1) + (j_f77 - 1) * ldab]);
    }
    printf ("\n");
  }
#endif

  /* Factorize the matrix */
  dgbtrf_ (&m_rows, &n_cols, &kl, &ku, AB, &ldab, ipiv, &info);

  if (info != 0)
    printf ("info = %d in Build_A() \n", info);

  free (u1);
}

/*===========================================================================*/

static int
Convert_indx (int i_f77)
{
  /* Given fortran "row" index i_f77, return index for u-type arrays
   * by first computing j and k (r and theta indicies). */

  int j, k, return_val;

  if (i_f77 == 1) {
    return (0);
  }

  j = 1 + (i_f77 - 2) / ntheta;
  k = (i_f77 - 2) - (j - 1) * ntheta;
  return_val = INDEX (j, k);

  /* bounds check */
#if 0
  if ((return_val < 0) || (return_val >= ntot)) {
    printf ("bad return_val in Convert_indx. return_val = %d\n", return_val);
    exit (1);
  }
#endif

  return (return_val);
}

/*===========================================================================*/
