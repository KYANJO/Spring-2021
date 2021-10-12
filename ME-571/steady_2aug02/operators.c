/* Include Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "driver.h"
#include "operators.h"

#define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )

#if FIFE_SCALING
#define U_KINETICS(u,uth) ( one_o_eps*one_o_eps*(u)*(1.-(u))*((u)-(uth)) )
#else
#define U_KINETICS(u,uth) ( one_o_eps*(u)*(1.-(u))*((u)-(uth)) )
#endif

/* Testing */
#define TEST_DERIVATIVES 0
#define TEST_KINETICS 0
#define AXISYMMETRIC 0

double *Ld, *Ll, *Lu;
double *omega_column;

/* 
 * Private variables 
 * ----------------- */
static double *dcoef, *gcoef, *ocoef;
static double *f_u, *f_v;

/* 
 * Private functions 
 * ----------------- */

static void Eval_Operators (char type, double *u, double *u_prime);
static void R_Derivatives (double *u, double *u_prime);
static void Theta_Derivatives (double *u, double *u_prime);
static void Linear_Kinetics (double *u, double *u_prime);
static void Kinetics (double *u, double *u_prime);
static void Domega (double *u, double *u_prime);

/*===========================================================================*/

void
Eval_RHS (double *u, double *rhs)
{
  /* u and rhs may be the same */
  Eval_Operators ('N', u, rhs);
}

/*===========================================================================*/

void
Eval_Au (double *u, double *A_on_u)
{
  int j;
  double omega;

  /* We do this since omega is a variable */
  omega = u[INDEX (r_fix, theta_fix)];
  u[INDEX (r_fix, theta_fix)] = 0.;

  Eval_Operators ('L', u, A_on_u);

  /* We do this since omega is a variable */
  for (j = 0; j < ntot; j++) {
    A_on_u[j] += omega * omega_column[j];
  }
  u[INDEX (r_fix, theta_fix)] = omega;
}

/*===========================================================================*/

void
To_Build (double *u, double *u_prime)
{
  double *a_t, *a_k;
  int j;

  a_t = VECTOR (ntot);
  a_k = VECTOR (ntot);
  Theta_Derivatives (u, a_t);
  Linear_Kinetics (u, a_k);
  for (j = 0; j < ntot; j++) {
    u_prime[j] = a_t[j] + a_k[j];
  }
  B_op (u_prime);
  free (a_t);
  free (a_k);
}

/*===========================================================================*/

static void
Eval_Operators (char type, double *u, double *u_prime)
{
  /* u and u_prime may be the same */

  double *a_r, *a_t, *a_k;
  int j;

  a_r = VECTOR (ntot);
  a_t = VECTOR (ntot);
  a_k = VECTOR (ntot);

  R_Derivatives (u, a_r);
  Theta_Derivatives (u, a_t);

  if (type == 'L') {
    Linear_Kinetics (u, a_k);
  }
  else if (type == 'N') {
    Kinetics (u, a_k);
  }
  else {
    printf ("unknown type in Eval_Operators(): type=%c\n", type);
  }

#if TEST_KINETICS
  for (j = 0; j < ntot; j++) {
    a_t[j] = 0.;
  }
  for (j = 0; j < ntot; j++) {
    a_r[j] = 0.;
  }
#endif

#if TEST_DERIVATIVES
  for (j = 0; j < ntot; j++) {
    a_k[j] = 0.;
  }
#endif

  for (j = 0; j < ntot; j++) {
    u_prime[j] = a_r[j] + a_t[j] + a_k[j];
  }

  B_op (u_prime);

  free (a_r);
  free (a_t);
  free (a_k);
}

/*===========================================================================*/

static void
R_Derivatives (double *u, double *u_prime)
{
  /* r derivatives plus Laplacian at center */

  int j, k;

  /* Laplacian at center */
#if AXISYMMETRIC
  u_prime[INDEX (0, 0)] = Ld[0] * u[INDEX (0, 0)] + Lu[0] * u[INDEX (1, 0)];
#else
  /* "average around the circle j=1 and subtract center" */
  u_prime[INDEX (0, 0)] = Ld[0] * u[INDEX (0, 0)];
  for (k = 0; k < ntheta; k++) {
    u_prime[INDEX (0, 0)] += Lu[0] * u[INDEX (1, k)];
  }
#endif

  /* r derivatives at all interior radii */
  for (j = 1; j < nr; j++) {
    for (k = 0; k < ntheta; k++) {
      u_prime[INDEX (j, k)]
	= Ll[j - 1] * u[INDEX (j - 1, k)]
	+ Ld[j] * u[INDEX (j, k)]
	+ Lu[j] * u[INDEX (j + 1, k)];
    }
  }

  /* r derivatives at radius */
  for (k = 0; k < ntheta; k++) {
    u_prime[INDEX (nr, k)]
      = Ll[nr - 1] * u[INDEX (nr - 1, k)]
      + Ld[nr] * u[INDEX (nr, k)];
  }
}

/*===========================================================================*/

static void
Theta_Derivatives (double *u, double *u_prime)
{
  int ntheta_o_2 = ntheta / 2, j, m;
  double real_part, imag_part;

  for (j = 0; j < ntot; j++) {
    u_prime[j] = u[j];
  }

  /* Added Dec 01 */
  for (j = 0; j < ntheta; j++) {
    u_prime[j] = 0.;
  }

  for (j = 1; j <= nr; j++) {
    realft (u_prime + INDEX (j, 0), ntheta_o_2, 1);

    u_prime[INDEX (j, 0)] *= dcoef[INDEX (j, 0)];
    u_prime[INDEX (j, 1)] *= dcoef[INDEX (j, 1)];

    for (m = 1; m < ntheta_o_2; m++) {
      real_part = u_prime[INDEX (j, 2 * m)] * dcoef[INDEX (j, 2 * m)]
	- u_prime[INDEX (j, 2 * m + 1)] * dcoef[INDEX (j, 2 * m + 1)];
      imag_part = u_prime[INDEX (j, 2 * m)] * dcoef[INDEX (j, 2 * m + 1)]
	+ u_prime[INDEX (j, 2 * m + 1)] * dcoef[INDEX (j, 2 * m)];

      u_prime[INDEX (j, 2 * m)] = real_part;
      u_prime[INDEX (j, 2 * m + 1)] = imag_part;
    }

    realft (u_prime + INDEX (j, 0), ntheta_o_2, -1);
  }
}

/*===========================================================================*/

static void
Linear_Kinetics (double *u, double *u_prime)
{
  int j;
  double *v = VECTOR (ntot);

  Green (u, v);

  for (j = 0; j < ntot; j++) {
    u_prime[j] = f_u[j] * u[j] + f_v[j] * v[j];
  }

  free (v);
}

/*===========================================================================*/

static void
Kinetics (double *u, double *u_prime)
{
  int j;
  double *v = VECTOR (ntot), uth;

  Green (u, v);

  for (j = 0; j < ntot; j++) {
    uth = U_THRESHOLD (v[j]);
    u_prime[j] = U_KINETICS (u[j], uth);
  }

  free (v);
}

/*===========================================================================*/

void
Green (double *u, double *u_prime)
{
  int ntheta_o_2 = ntheta / 2, j, m;
  double real_part, imag_part;

  for (j = 0; j < ntot; j++) {
    u_prime[j] = u[j];
  }

  /* Above sets v_center=u_center.  Now all other radii */
  for (j = 1; j <= nr; j++) {
    realft (u_prime + INDEX (j, 0), ntheta_o_2, 1);

    u_prime[INDEX (j, 0)] *= gcoef[0];
    u_prime[INDEX (j, 1)] *= gcoef[1];

    for (m = 1; m < ntheta_o_2; m++) {
      real_part = u_prime[INDEX (j, 2 * m)] * gcoef[2 * m]
	- u_prime[INDEX (j, 2 * m + 1)] * gcoef[2 * m + 1];
      imag_part = u_prime[INDEX (j, 2 * m)] * gcoef[2 * m + 1]
	+ u_prime[INDEX (j, 2 * m + 1)] * gcoef[2 * m];

      u_prime[INDEX (j, 2 * m)] = real_part;
      u_prime[INDEX (j, 2 * m + 1)] = imag_part;
    }

    realft (u_prime + INDEX (j, 0), ntheta_o_2, -1);
  }
}

/*===========================================================================*/

void
Domega (double *u, double *u_prime)
{
  int ntheta_o_2 = ntheta / 2, j, m;
  double real_part, imag_part;
  double *u_t = VECTOR (ntot);	/* contribution from d/d theta */
  double *u_k = VECTOR (ntot);	/* contribution from kinetics through
				 * Greens function */

  for (j = 0; j < ntot; j++) {
    u_t[j] = u[j];
    u_k[j] = u[j];
  }

  /* Set values at center = 0 */
  u_t[0] = 0.;
  u_k[0] = 0.;

  /* Now all other radii */
  for (j = 1; j <= nr; j++) {
    realft (u_t + INDEX (j, 0), ntheta_o_2, 1);
    realft (u_k + INDEX (j, 0), ntheta_o_2, 1);

    u_t[INDEX (j, 0)] = 0.;
    u_t[INDEX (j, 1)] = 0.;

    u_k[INDEX (j, 0)] *= ocoef[INDEX (j, 0)];
    u_k[INDEX (j, 1)] *= ocoef[INDEX (j, 1)];

    /* u_t is simply d/d theta of u */
    for (m = 1; m < ntheta_o_2; m++) {
      real_part = u_t[INDEX (j, 2 * m + 1)] * m / (double) ntheta_o_2;
      imag_part = -u_t[INDEX (j, 2 * m)] * m / (double) ntheta_o_2;
      u_t[INDEX (j, 2 * m)] = real_part;
      u_t[INDEX (j, 2 * m + 1)] = imag_part;
    }

    for (m = 1; m < ntheta_o_2; m++) {
      real_part = u_k[INDEX (j, 2 * m)] * ocoef[INDEX (j, 2 * m)]
	- u_k[INDEX (j, 2 * m + 1)] * ocoef[INDEX (j, 2 * m + 1)];
      imag_part = u_k[INDEX (j, 2 * m)] * ocoef[INDEX (j, 2 * m + 1)]
	+ u_k[INDEX (j, 2 * m + 1)] * ocoef[INDEX (j, 2 * m)];

      u_k[INDEX (j, 2 * m)] = real_part;
      u_k[INDEX (j, 2 * m + 1)] = imag_part;
    }

    realft (u_t + INDEX (j, 0), ntheta_o_2, -1);
    realft (u_k + INDEX (j, 0), ntheta_o_2, -1);
  }

  for (j = 0; j < ntot; j++) {
    u_prime[j] = u_t[j] + f_v[j] * u_k[j];
  }
  B_op (u_prime);

  free (u_t);
  free (u_k);
}

/*===========================================================================*/

void
B_op (double *u)
{
  int j, k;

  /* Set "theta=2pi" line to "theta=0" line */
  for (j = 1; j <= nr; j++) {
    u[INDEX (j, ntheta)] = u[INDEX (j, 0)];
  }

  /* Set all "extra" center (points r=0, theta>0) to center (r=0, theta=0) */
  for (k = 1; k <= ntheta; k++) {
    u[INDEX (0, k)] = u[INDEX (0, 0)];
  }
}

/*===========================================================================*/

void
Make_Operators (double *u_base, double omega_base, double dr)
{
  int j, m, ntheta_o_2 = ntheta / 2;
  double one_o_dr2 = 1. / (dr * dr);
  double one_o_2dr = 1. / (2. * dr);
  double real_part, imag_part, denominator;
  double *v_base, uth;

  /* Set L matrix elements: r-derivatives */

  for (j = 1; j < nr; j++) {
    Ll[j - 1] = one_o_dr2 - one_o_2dr / R (j);
    Ld[j] = -2. * one_o_dr2;
    Lu[j] = one_o_dr2 + one_o_2dr / R (j);
  }

  /* Take into account BCs */
#if AXISYMMETRIC
  Ld[0] = -4. * one_o_dr2;	/* these two used only for */
  Lu[0] = 4. * one_o_dr2;	/* axisymmetric case */
#else
  Ld[0] = -4. * one_o_dr2;
  Lu[0] = 4. * one_o_dr2 / ntheta;
#endif

  Ld[nr] = -2. * one_o_dr2;
  Ll[nr - 1] = 2. * one_o_dr2;


  /* Set Fourier coefficients: theta derivatives */

  for (j = 1; j <= nr; j++) {
    for (m = 0; m <= ntheta_o_2; m++) {

      /*      real_part = -m * m / (R (j) * R (j)); */
      real_part = -m * m * (tau2 + 1. / (R (j) * R (j)));
      imag_part = m * omega_base;
      real_part /= (double) ntheta_o_2;
      imag_part /= -(double) ntheta_o_2;
      /* the minus sign above is because NumRec fft gives fourier 
       * coefs of e^{-i m theta} */

      if (m == 0) {
	dcoef[INDEX (j, 0)] = real_part;
      }
      else if (m == ntheta_o_2) {
	dcoef[INDEX (j, 1)] = real_part;
      }
      else {
	dcoef[INDEX (j, 2 * m)] = real_part;
	dcoef[INDEX (j, 2 * m + 1)] = imag_part;
      }
    }
  }

#if FIFE_SCALING
  omega_base *= one_o_eps;
#endif

  /* Set Fourier coefficients: d/domega terms */

  for (j = 1; j <= nr; j++) {
    for (m = 0; m <= ntheta_o_2; m++) {

      denominator = pow ((1. + m * m * omega_base * omega_base), 2.);
      real_part = -2. * m * m * omega_base / denominator;
      imag_part = m * (1. - m * m * omega_base * omega_base) / denominator;

      real_part /= (double) ntheta_o_2;
      imag_part /= -(double) ntheta_o_2;
      /* the minus sign above is because NumRec fft gives fourier 
       * coefs of e^{-i m theta} */

#if FIFE_SCALING
      real_part *= one_o_eps;
      imag_part *= one_o_eps;
#endif

      if (m == 0) {
	ocoef[INDEX (j, 0)] = real_part;
      }
      else if (m == ntheta_o_2) {
	ocoef[INDEX (j, 1)] = real_part;
      }
      else {
	ocoef[INDEX (j, 2 * m)] = real_part;
	ocoef[INDEX (j, 2 * m + 1)] = imag_part;
      }
    }
  }

  /* Set coefficients for Greens function */

  for (m = 0; m <= ntheta_o_2; m++) {

    real_part = 1. / (1. + m * m * omega_base * omega_base);
    imag_part = m * omega_base / (1. + m * m * omega_base * omega_base);
    real_part /= (double) ntheta_o_2;
    imag_part /= -(double) ntheta_o_2;
    /* the minus sign above is because NumRec fft gives fourier 
     * coefs of e^{-i m theta} */

    if (m == 0) {
      gcoef[0] = real_part;
    }
    else if (m == ntheta_o_2) {
      gcoef[1] = real_part;
    }
    else {
      gcoef[2 * m] = real_part;
      gcoef[2 * m + 1] = imag_part;
    }
  }

  /* Set elements of Linear_Kinetics (note call to Green below must
   * come after gcoef is set above) */

  v_base = VECTOR (ntot);
  Green (u_base, v_base);
  for (j = 0; j < ntot; j++) {
    uth = U_THRESHOLD (v_base[j]);
    f_u[j] = one_o_eps * ((1. - u_base[j]) * (u_base[j] - uth)
			  - u_base[j] * (u_base[j] - uth)
			  + u_base[j] * (1. - u_base[j]));
    f_v[j] = one_o_eps * u_base[j] * (1. - u_base[j]) * (-one_o_a);

#if FIFE_SCALING
    f_u[j] *= one_o_eps;
    f_v[j] *= one_o_eps;
#endif
  }
  free (v_base);

  /* Compute d_domega terms (must come after f_v is set above) */
  Domega (u_base, omega_column);
}

/*===========================================================================*/

void
Operators_ini (void)
{
  /* Memory allocation */

  dcoef = VECTOR (ntot);
  gcoef = VECTOR (ntheta);
  ocoef = VECTOR (ntot);
  Ld = VECTOR (nr + 1);
  Ll = VECTOR (nr);
  Lu = VECTOR (nr);
  f_u = VECTOR (ntot);
  f_v = VECTOR (ntot);
  omega_column = VECTOR (ntot);

  /* Could also set L matrix coefficients because there are
   * independent of u_base and omega_base, but set them in Make_Ops
   * with other coefficients. */
}

/*===========================================================================*/

void
Operators_fin (void)
{
  free (dcoef);
  free (gcoef);
  free (ocoef);
  free (Ld);
  free (Ll);
  free (Lu);
  free (f_u);
  free (f_v);
  free (omega_column);
}

/*===========================================================================*/
