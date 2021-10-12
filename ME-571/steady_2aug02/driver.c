/*---------------------------------------------------------------------------*
 *                                                                           
 *                                   Steady
 *     The code computes steady spiral and straight twist scroll solutions
 *                                                                           
 *                     Copyright (C) 2002  Dwight Barkley
 *                   	 barkley@maths.warwick.ac.uk
 *                                                                           
 *                                Version 1.0
 *                                                                           
 *---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "driver.h"
#include "solver.h"
#include "operators.h"

#if GRAPHICS
#include "ezplot.h"
#endif

/* Set model and numerical parameter values */

#define A_PARAM     1.0
#define B_PARAM     0.1
#define ONE_O_EPS   50.0
#define TAU         0.0
#define RADIUS      10.0
#define NR          100		/* 0 for read at run time */
#define NTHETA      128		/* 0 for read at run time */

/* Set the continuation parameter PARAMETER. Normally one of: 
 * a_param, b_param, one_o_eps, radius, tau.  
 *
 * Continuation ends when the continuation parameter reaches the value
 * PARAM_END. STEP_INI, STEP_INI, STEP_INI should all be positive and
 * set the initial, minimum, and maximum parameter step size.
 *
 * CONTINUOUS 1 allows for the step size to vary continuously to keep
 * the number of Newton iterations about 3.  Else the step size is
 * fixed except that is may be decreased in the event of failure of
 * Newton iterations.
 * */

#define PARAMETER   a_param
#define PARAM_END   1.0
#define STEP_INI    0.001       /* > 0 */
#define STEP_MIN    0.0001	/* > 0 */
#define STEP_MAX    0.5		/* > 0 */
#define CONTINUOUS  1

/* Parameters for Newton iteration. Success if norm of rhs <
 * TOL. Failure is norm exceeds MAX_NORM or more the MAX_ITS
 * iterations are taken. */

#define TOL         1.e-3
#define MAX_NORM    1000.
#define MAX_ITS     10


/* Globals */
int nr, ntheta, tinc, ntot;
int r_fix, theta_fix;
double one_o_eps, one_o_a, b_o_a, tau2;

/* 
 * Private variables 
 * ----------------- */
static double *u, *rhs, omega, dr, dtheta;
static double a_param, b_param, tau, radius, param_step;
static FILE *omega_file;
static int sign, one = 1;

/* 
 * Private functions 
 * ----------------- */
static double Newton_ini (void);
static double Newton (void);
static void Initialize (void);
static void Finish_up (void);
static void Read_ic (void);
static void Write_fc (void);

/*===========================================================================*/

int
main (void)
{
  int j, nstep, first_time = 1;
  double *u_save, omega_save, norm;
  double new_param_step, new_u, new_omega;

  Initialize ();

  /* Save solution for recovery from a bad Newton step */
  u_save = VECTOR (ntot);
  for (j = 0; j < ntot; j++) {
    u_save[j] = u[j];
  }
  omega_save = omega;

#if GRAPHICS
  /* Plot initial condition */
  EZplot_polar (u, dr, dtheta, nr, ntheta);
  EZevent ();
#endif

  /* Start Continuation Loop */

  while (sign * PARAMETER <= sign * ((PARAM_END) + 0.5 * param_step)) {

    /* Compute these in case they involve PARAMETER */
    tau2 = tau * tau;
    dr = radius / nr;
    one_o_a = 1. / a_param;
    b_o_a = b_param * one_o_a;

    /* Initialize Newton iteration */
    nstep = 0;
    norm = Newton_ini ();
    printf ("   step = %d", nstep);
    printf ("   norm = %g", norm);
    printf ("   omega = %g \n", omega);

    /* Start Newton iterations */
    while ((nstep < MAX_ITS) && (norm > TOL) && (norm < MAX_NORM)) {
      nstep++;
      norm = Newton ();
      printf ("   step = %d", nstep);
      printf ("   norm = %g", norm);
      printf ("   omega = %g \n", omega);

#if GRAPHICS
      /* Plot new solution */
      EZplot_polar (u, dr, dtheta, nr, ntheta);
      if (EZevent ()) {
	Finish_up ();
	exit (1);
      }
#endif
    }
    /* Newton iterations have finished, either successfully or not */

    if (norm < TOL) {
      /* Success */
      printf ("PARAMETER = %g, omega = %g\n\n", PARAMETER, omega);
      fprintf (omega_file, "%g  %g \n", PARAMETER, omega);

      if (PARAMETER == PARAM_END)
	break;

      if (CONTINUOUS)
	new_param_step = param_step * min (3. / nstep, 1.5);
      else
	new_param_step = param_step;
    }
    else {
      /* Failure.  */
      if (first_time) {
	/* Could not even get started so just exit */
	exit (1);
      }

      /* Recover saved values and adjust param_step */
      for (j = 0; j < ntot; j++) {
	u[j] = u_save[j];
      }
      omega = omega_save;
      PARAMETER -= param_step;

      if (fabs (param_step) == STEP_MIN) {
	/* nothing else to try */
	Finish_up ();
	exit (1);
      }
      else {
	/* Reduce parameter step and try again */
	new_param_step = param_step / 2.;
      }
    }

    first_time = 0;

    /* Adjust new_param_step */
    new_param_step = sign * max (fabs (new_param_step), STEP_MIN);
    if (fabs (PARAMETER - (PARAM_END)) < fabs (new_param_step)) {
      new_param_step = (PARAM_END) - PARAMETER;
    }
    printf ("  param_step = %g\n", new_param_step);

    /* Extrapolate and save u, omega for next time or in case of failure. */
    for (j = 0; j < ntot; j++) {
      new_u = u[j] + new_param_step / param_step * (u[j] - u_save[j]);
      u_save[j] = u[j];
      u[j] = new_u;
    }
    new_omega = omega + new_param_step / param_step * (omega - omega_save);
    omega_save = omega;
    omega = new_omega;
    param_step = new_param_step;

    PARAMETER += param_step;
    if (fabs ((PARAMETER - (PARAM_END)) / param_step) < 0.1)
      PARAMETER = PARAM_END;
  }
  /* End Continuation Loop */

  Finish_up ();

  return (0);
}

/*===========================================================================*/

static double
Newton_ini (void)
{
  double norm;

  Make_Operators (u, omega, dr);
  Eval_RHS (u, rhs);
  norm = dnrm2_ (&ntot, rhs, &one);
  return (norm);
}

/*===========================================================================*/

static double
Newton (void)
{
  int j;
  double domega, norm;

  /* Take one Newton step. Note, virtually all the work of this
   * code in direct-solve mode occurs with the call to Build_A. 
   * */

  Build_A ();
  Solve_A (rhs);

  domega = rhs[INDEX (r_fix, theta_fix)];
  rhs[INDEX (r_fix, theta_fix)] = 0.;

  omega = omega - domega;
  for (j = 0; j < ntot; j++) {
    u[j] = u[j] - rhs[j];
  }
  Make_Operators (u, omega, dr);
  Eval_RHS (u, rhs);
  norm = dnrm2_ (&ntot, rhs, &one);

  return (norm);
}

/*===========================================================================*/

static void
Initialize (void)
{
  /* Section for setting and writing problem parameters 
   * -------------------------------------------------- */

  /* set global integers */
#if NR
  nr = NR;
#else
  printf ("nr = ? ");
  scanf ("%d", &nr);
#endif

#if NTHETA
  ntheta = NTHETA;
#else
  printf ("ntheta = ? ");
  scanf ("%d", &ntheta);
#endif

  tinc = ntheta + 1;
  ntot = (ntheta + 1) * (nr + 1);

  /* Set model and numerical parameters */

  a_param = A_PARAM;
  b_param = B_PARAM;
  one_o_eps = ONE_O_EPS;
  radius = RADIUS;
  tau = TAU;

  one_o_a = 1. / a_param;
  b_o_a = b_param * one_o_a;
  tau2 = tau * tau;
  dr = radius / nr;
  dtheta = 2 * M_PI / ntheta;

  if ((PARAM_END - PARAMETER) >= 0.) {
    sign = 1;
  }
  else {
    sign = -1;
  }
  param_step = sign * STEP_INI;

  /* Allocate memory for solution variable(s) */
  u = VECTOR (ntot);
  rhs = VECTOR (ntot);

  /* Read initial conditions */
  Read_ic ();
  B_op (u);

  /* Set point for fixed u (I prefer r_fix=nr/2 and theta_fix=0) */
  r_fix = nr / 2;
#if 0 
  { /* this code is useful if theta_fix=0 is not appropriate */
    int k;
    for (k = 0; k < ntheta; k++) {
      if ((u[INDEX (r_fix, k)] - 0.5) <= 0 &&
	  (u[INDEX (r_fix, k + 1)] - 0.5) > 0) {
	theta_fix = k;
	break;
      }
    }
  }
#endif
  theta_fix = 0; 
  u[INDEX (r_fix, theta_fix)] = 0.5;

  /* Initialize any other sections of the code */
  Operators_ini ();
#if GRAPHICS
  /* EZplot initialization */
  EZplot_ini (-radius - 1, radius + 1, -radius - 1, radius + 1);
#endif

  /* Write to screen */
  printf ("\n\nModel Parameters: \n");
  printf ("a     = %g\n", a_param);
  printf ("b     = %g\n", b_param);
  printf ("eps   = 1/%g = %g\n", one_o_eps, 1.0 / one_o_eps);
  printf ("tau2  = %g \n", tau2);
  printf ("Radius = %g\n", radius);

  printf ("\nNumerical Parameters: \n");
  printf ("N_r = %d,  N_theta = %d, dr = %g, dtheta = %g \n",
	  nr, ntheta, dr, dtheta);
  printf ("r_fix = %d,   theta_fix  = %d,  u_fix = %g \n",
	  r_fix, theta_fix, u[INDEX (r_fix, theta_fix)]);
  printf ("maximum iterations = %d \n\n", MAX_ITS);

  printf ("Continuing to PARAMETER = %g\n\n", PARAM_END);

  /* initialize output file for omega */
  omega_file = fopen ("omega.dat", "w");
  fprintf (omega_file, "# Model Parameters:  ");
  fprintf (omega_file, "a = %g  ", a_param);
  fprintf (omega_file, "b = %g  ", b_param);
  fprintf (omega_file, "eps = 1/%g = %g  ", one_o_eps, 1.0 / one_o_eps);
  fprintf (omega_file, "tau2 = %g ", tau2);
  fprintf (omega_file, "Radius = %g", radius);
  fprintf (omega_file, "\n# Numerical Parameters:  ");
  fprintf (omega_file, "N_r = %d,  N_theta = %d \n", nr, ntheta);
  fprintf (omega_file, "# PARAMETER    omega   \n");
}

/*===========================================================================*/

static void
Finish_up (void)
{
  /* Write final results */
  Write_fc ();
}

/*===========================================================================*/

#define NEXT_LINE(fp) while(getc(fp)!='\n');	/* macro used to skip to end 
						 * * of input line */

static void
Read_ic (void)
{
  /* Reads initial condition file (ic.dat) */

  double u_in, v_in, u_old;
  int nr_ic, ntheta_ic, j, k;
  char f_type;
  FILE *fp;

  if ((fp = fopen ("ic.dat", "r")) == NULL) {
    printf ("no ic.dat\n");
    exit (1);
  }

  /* Read nr_ic etc following = sign on second line of file */
  NEXT_LINE (fp);
  while ('=' != getc (fp));
  fscanf (fp, "%d, %d", &nr_ic, &ntheta_ic);

  /* Skip to 9th line */
  for (j = 0; j < 7; j++) {
    NEXT_LINE (fp);
  }

  while ('=' != getc (fp));
  fscanf (fp, "%lg", &omega);
  NEXT_LINE (fp);

  f_type = getc (fp);
  NEXT_LINE (fp);
  if ((f_type != 'B') && (f_type != 'A')) {
    printf ("\n ic.dat exists but of unrecognized type Binary or Ascii \n");
    exit (1);
  }

  if (f_type == 'B') {
    /* Binary data file */
    if ((nr_ic != nr) || (ntheta_ic != ntheta)) {
      printf ("\n ic.dat has nr, ntheta = %d, %d \n\n", nr_ic, ntheta_ic);
      exit (1);
    }
    fread (u, sizeof (double), ntot, fp);
  }
  else {
    /* Ascii data file */

    /* ic size = current size */
    if ((nr_ic == nr) && (ntheta_ic == ntheta)) {
      for (j = 0; j < ntot; j++) {
	fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	u[j] = u_in;
      }
    }

    /* ic with half current ntheta */
    else if ((nr_ic == nr) && (2 * ntheta_ic == ntheta)) {
      u_old = 0.;
      for (j = 0; j <= nr; j++) {
	for (k = 0; k < ntheta_ic; k++) {
	  fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	  u[INDEX (j, 2 * k)] = (u_in + u_old) / 2.;
	  u[INDEX (j, 2 * k + 1)] = u_in;
	  u_old = u_in;
	}
	fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	u_old = u_in;
      }
    }

    /* ic with twice current ntheta */
    else if ((nr_ic == nr) && (ntheta_ic == 2 * ntheta)) {
      for (j = 0; j <= nr; j++) {
	for (k = 0; k < ntheta; k++) {
	  fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	  fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	  u[INDEX (j, k)] = u_in;
	}
	fscanf (fp, "%lg %lg\n", &u_in, &v_in);
      }
    }

    /* ic with half current nr */
    else if ((2 * nr_ic == nr) && (ntheta_ic == ntheta)) {
      for (j = 0; j <= nr_ic; j++) {
	for (k = 0; k <= ntheta_ic; k++) {
	  fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	  u[INDEX (2 * j, k)] = u_in;
	  if (j != 0)
	    u[INDEX (2 * j - 1, k)] =
	      0.5 * (u[INDEX (2 * j, k)] + u[INDEX (2 * j - 2, k)]);
	}
      }
    }

    /* ic with twice current nr */
    else if ((nr_ic == 2 * nr) && (ntheta_ic == ntheta)) {
      for (j = 0; j <= nr; j++) {
	for (k = 0; k <= ntheta; k++) {
	  fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	  u[INDEX (j, k)] = u_in;
	}
	if (j != nr) {
	  for (k = 0; k <= ntheta; k++) {
	    fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	  }
	}
      }
    }

    /* ic with nr <  current nr */
    else if ((nr > nr_ic) && (ntheta_ic == ntheta)) {
      for (j = 0; j <= nr_ic; j++) {
	for (k = 0; k <= ntheta_ic; k++) {
	  fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	  u[INDEX (j, k)] = u_in;
	}
      }
      for (j = nr_ic + 1; j <= nr; j++) {
	for (k = 0; k <= ntheta_ic; k++) {
	  u[INDEX (j, k)] = u[INDEX (nr_ic, k)];
	}
      }
    }

    else {
      printf ("\n ic.dat has nr, ntheta = %d, %d \n\n", nr_ic, ntheta_ic);
      exit (1);
    }
  }
}

/*===========================================================================*/

static void
Write_fc (void)
{
  /* Write final condition file */

  FILE *fp;
  double *v = VECTOR (ntot);
  time_t tt1;
  int j;

  time (&tt1);
  fp = fopen ("fc.dat", "w");

  fprintf (fp, "Model Parameters: a, b, 1/eps, tau, radius = ");
  fprintf (fp, "%g, %g, %g, %g, %g \n", a_param, b_param, one_o_eps,
	   sqrt (tau2), radius);
  fprintf (fp, "Numerical Parameters: nr, ntheta = ");
  fprintf (fp, "%d, %d \n", nr, ntheta);
  fprintf (fp, "File written: %s", ctime (&tt1));

  if (FIFE_SCALING)
    fprintf (fp, "Comments: Fife scaling \n");
  else
    fprintf (fp, "Comments: \n");

  fprintf (fp, "\n");

  fprintf (fp, "\n");
  fprintf (fp, "\n");
  fprintf (fp, "\n");
  fprintf (fp, "omega = %g\n", omega);

  Green (u, v);

  /* Write ascii data */
  fprintf (fp, "Ascii values of u and v follow\n");
  for (j = 0; j < ntot; j++) {
    fprintf (fp, "%.12g %.12g\n", u[j], v[j]);
  }

#if 0				/* code to change rotation direction */
  {
    int k;
    for (j = 0; j <= nr; j++) {
      for (k = 0; k <= ntheta; k++)
	fprintf (fp, "%.12g %.8g\n", u[INDEX (j, ntheta - k)],
		 v[INDEX (j, ntheta - k)]);
    }
  }
#endif

  fclose (fp);
  free (v);
}

/*===========================================================================*/

double
R (int j)
{
  return (j * dr);
}

/*===========================================================================*/

double
Th (int k)
{
  return (k * dtheta);
}

/*===========================================================================*/

void
MAIN__ (void)
{				/* dummy for fort77 */
}

/*===========================================================================*/
