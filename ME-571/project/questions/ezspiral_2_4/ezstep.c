/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                                EZSTEP.C                                   */
/*                  Time-stepping Subroutines for EZ-SPIRAL                  */
/*                                                                           */
/*                                Version 2                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
 * $Revision: 2.3 $
 * $Date: 1997/11/13 23:18:21 $
 */
/*---------------------------------------------------------------------------*/

#include "ezspiral.h"
#include "ezstep.h"

/* External variables */
extern precision  **u, **v, **u_lap[2], **v_lap[2], 
                  delta, dt, NDu, NDv, one_o_a, b_o_a, dt_o_eps, dt_o_2eps; 
extern int        ndim, k1, k2;

/* Functions local to ezstep */
static void       nbc();
static void       pbc();

/* constants */
const precision zero   = 0.;
const precision half   = 0.5;
const precision one    = 1.;
const precision two    = 2.;
const precision three  = 3.;
const precision four   = 4.;
const precision twenty = 20.;

/*---------------------------------------------------------------------------*/

void step()

     /* 
      *  Routine for taking one time step. 
      *  See ezstep.h for definitions of various macros used here.
      */
{
  precision u_thresh;
  register int i, j;
  
  /* Interchange k1 and k2 */
  k1 = 1-k1;
  k2 = 1-k2;
  
  /* Main Loop (ALMOST ALL WORK DONE HERE !!!) */
  for(i=1;i<=ndim;i++) {
    for(j=1;j<=ndim;j++) {
      if(u[i][j]<delta) { 
	u[i][j] = U_DIFFUSION;
	v[i][j] = V_KINETICS(0,v[i][j]) + V_DIFFUSION;
      }
      else { 
	u_thresh = U_THRESHOLD(v[i][j]);
	v[i][j]  = V_KINETICS(u[i][j],v[i][j])  + V_DIFFUSION;
	u[i][j]  = U_KINETICS(u[i][j],u_thresh) + U_DIFFUSION;
	ADD_TO_U_LAPLACIAN;
      }
      ADD_TO_V_LAPLACIAN;
      ZERO_USED_LAPLACIANS;
    }
  }

  /* Impose boundary conditions */
  if(NBC) {
    nbc();
  } 
  else {
    pbc();
  }
}    
/*---------------------------------------------------------------------------*/

void step_ini()

     /* 
      *  Routine for initializing Laplacians. 
      */
{
  int i, j;
  
  /* Set initial k1 and k2 */
  k1 = 0;
  k2 = 1;
  
  /* Compute Laplacians of initial fields */
  for(i=1;i<=ndim;i++) {
    for(j=1;j<=ndim;j++) {
      ADD_TO_U_LAPLACIAN;
      ADD_TO_V_LAPLACIAN;
      ZERO_USED_LAPLACIANS;
    }
  }
  
  /* Impose boundary conditions */
  if(NBC) {
    nbc();
  } 
  else {
    pbc();
  }
}
/*---------------------------------------------------------------------------*/

void nbc()

     /* 
      *  Imposes Neumann boundary conditions. 
      */
{
  int ib;

  for(ib=1;ib<=ndim;ib++) {
#if NINEPOINT  /* Nine-point Laplacian formulas */
    u_lap[k2][ib][1]    += four*u[ib][2];
    u_lap[k2][ib][ndim] += four*u[ib][ndim-1];
    u_lap[k2][1][ib]    += four*u[2][ib];
    u_lap[k2][ndim][ib] += four*u[ndim-1][ib];
    if(ib==1) {
      u_lap[k2][1][1]    += three*u[2][2];
      u_lap[k2][1][ndim] += three*u[2][ndim-1];
      u_lap[k2][ndim][1] += three*u[ndim-1][2];
    }
    else if (ib==ndim) {
      u_lap[k2][ndim][ndim] += three*u[ndim-1][ndim-1];
    }
    else {
      u_lap[k2][ib][1]    += u[ib-1][2];
      u_lap[k2][ib][1]    += u[ib+1][2];
      u_lap[k2][ib][ndim] += u[ib-1][ndim-1];
      u_lap[k2][ib][ndim] += u[ib+1][ndim-1];
      u_lap[k2][1][ib]    += u[2][ib-1];
      u_lap[k2][1][ib]    += u[2][ib+1];
      u_lap[k2][ndim][ib] += u[ndim-1][ib-1];
      u_lap[k2][ndim][ib] += u[ndim-1][ib+1];
    }
#  if V_DIFF_ON  
    v_lap[k2][ib][1]    += four*v[ib][2];
    v_lap[k2][ib][ndim] += four*v[ib][ndim-1];
    v_lap[k2][1][ib]    += four*v[2][ib];
    v_lap[k2][ndim][ib] += four*v[ndim-1][ib];
    if(ib==1) {
      v_lap[k2][1][1]    += three*v[2][2];
      v_lap[k2][1][ndim] += three*v[2][ndim-1];
      v_lap[k2][ndim][1] += three*v[ndim-1][2];
    }
    else if (ib==ndim) {
      v_lap[k2][ndim][ndim] += three*v[ndim-1][ndim-1];
    }
    else {
      v_lap[k2][ib][1]    += v[ib-1][2];
      v_lap[k2][ib][1]    += v[ib+1][2];
      v_lap[k2][ib][ndim] += v[ib-1][ndim-1];
      v_lap[k2][ib][ndim] += v[ib+1][ndim-1];
      v_lap[k2][1][ib]    += v[2][ib-1];
      v_lap[k2][1][ib]    += v[2][ib+1];
      v_lap[k2][ndim][ib] += v[ndim-1][ib-1];
      v_lap[k2][ndim][ib] += v[ndim-1][ib+1];
    }
#  endif  /* end  #if V_DIFF_ON */
#else  /* Five-point Laplacian formulas */
    u_lap[k2][ib][1]    += u[ib][2];  
    u_lap[k2][ib][ndim] += u[ib][ndim-1]; 
    u_lap[k2][1][ib]    += u[2][ib]; 
    u_lap[k2][ndim][ib] += u[ndim-1][ib]; 
#  if V_DIFF_ON  
    v_lap[k2][ib][1]    += v[ib][2];  
    v_lap[k2][ib][ndim] += v[ib][ndim-1]; 
    v_lap[k2][1][ib]    += v[2][ib]; 
    v_lap[k2][ndim][ib] += v[ndim-1][ib]; 
#  endif  /* end  #if V_DIFF_ON */
#endif
  }
}
/*---------------------------------------------------------------------------*/

void pbc()

     /* 
      *  Imposes Periodic boundary conditions. 
      */
{
  int ib;

  for(ib=1;ib<=ndim;ib++) {
#if NINEPOINT  /* Nine-point Laplacian formulas */
    u_lap[k2][ib][1]    += four*u[ib][ndim];  
    u_lap[k2][ib][ndim] += four*u[ib][1]; 
    u_lap[k2][1][ib]    += four*u[ndim][ib]; 
    u_lap[k2][ndim][ib] += four*u[1][ib]; 
    if(ib==1) {
      u_lap[k2][1][1]    += u[2][ndim];
      u_lap[k2][1][1]    += u[ndim][2];
      u_lap[k2][1][1]    += u[ndim][ndim];
      u_lap[k2][ndim][1] += u[1][2];
      u_lap[k2][ndim][1] += u[1][ndim];
      u_lap[k2][1][ndim] += u[ndim][1];
      u_lap[k2][1][ndim] += u[2][1];
    }
    else if (ib==ndim) {
      u_lap[k2][ndim][1]    += u[ndim-1][ndim];
      u_lap[k2][1][ndim]    += u[ndim][ndim-1];
      u_lap[k2][ndim][ndim] += u[ndim-1][1];
      u_lap[k2][ndim][ndim] += u[1][ndim-1];
      u_lap[k2][ndim][ndim] += u[1][1];
    }
    else {
      u_lap[k2][ib][1]    += u[ib-1][ndim];
      u_lap[k2][ib][1]    += u[ib+1][ndim];
      u_lap[k2][ib][ndim] += u[ib-1][1];
      u_lap[k2][ib][ndim] += u[ib+1][1];
      u_lap[k2][1][ib]    += u[ndim][ib-1];
      u_lap[k2][1][ib]    += u[ndim][ib+1];
      u_lap[k2][ndim][ib] += u[1][ib-1];
      u_lap[k2][ndim][ib] += u[1][ib+1];
    }
#  if V_DIFF_ON  
    v_lap[k2][ib][1]    += four*v[ib][ndim];  
    v_lap[k2][ib][ndim] += four*v[ib][1]; 
    v_lap[k2][1][ib]    += four*v[ndim][ib]; 
    v_lap[k2][ndim][ib] += four*v[1][ib]; 
    if(ib==1) {
      v_lap[k2][1][1]    += v[2][ndim];
      v_lap[k2][1][1]    += v[ndim][2];
      v_lap[k2][1][1]    += v[ndim][ndim];
      v_lap[k2][ndim][1] += v[1][2];
      v_lap[k2][ndim][1] += v[1][ndim];
      v_lap[k2][1][ndim] += v[ndim][1];
      v_lap[k2][1][ndim] += v[2][1];
    }
    else if (ib==ndim) {
      v_lap[k2][ndim][1]    += v[ndim-1][ndim];
      v_lap[k2][1][ndim]    += v[ndim][ndim-1];
      v_lap[k2][ndim][ndim] += v[ndim-1][1];
      v_lap[k2][ndim][ndim] += v[1][ndim-1];
      v_lap[k2][ndim][ndim] += v[1][1];
    }
    else {
      v_lap[k2][ib][1]    += v[ib-1][ndim];
      v_lap[k2][ib][1]    += v[ib+1][ndim];
      v_lap[k2][ib][ndim] += v[ib-1][1];
      v_lap[k2][ib][ndim] += v[ib+1][1];
      v_lap[k2][1][ib]    += v[ndim][ib-1];
      v_lap[k2][1][ib]    += v[ndim][ib+1];
      v_lap[k2][ndim][ib] += v[1][ib-1];
      v_lap[k2][ndim][ib] += v[1][ib+1];
    }
#  endif  /* end  #if V_DIFF_ON */
#else  /* Five-point Laplacian formulas */
    u_lap[k2][ib][1]    += u[ib][ndim];  
    u_lap[k2][ib][ndim] += u[ib][1]; 
    u_lap[k2][1][ib]    += u[ndim][ib]; 
    u_lap[k2][ndim][ib] += u[1][ib]; 
#  if V_DIFF_ON  
    v_lap[k2][ib][1]    += v[ib][ndim];  
    v_lap[k2][ib][ndim] += v[ib][1]; 
    v_lap[k2][1][ib]    += v[ndim][ib]; 
    v_lap[k2][ndim][ib] += v[1][ib]; 
#  endif  /* end  #if V_DIFF_ON */
#endif
  }
}
/*---------------------------------------------------------------------------*/

