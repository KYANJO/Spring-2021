/*
 *  Define the macros U_THRESHOLD(v) and G(u,v) according to the model you 
 *  wish to simulate. In principle these can be almost anything you want.  
 *  Three examples are included.
 * 
 * $Revision: 2.3 $
 * $Date: 1997/11/13 23:17:22 $
 *---------------------------------------------------------------------------*/

#if 1  /* Standard Model */
#  define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
#  define G(u,v)          ( (u)-(v) )
#endif

#if 0  /* Simple Model for Turbulent Spirals */
#  define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
#  define G(u,v)          ( (u)*(u)*(u) - (v) )
#endif

#if 0  /*  Bar and Eiswirth Model (Phys. Rev. E V. 48, p. 1635, 1993).
	*  I do not include their case u>1 because it is not relevant.
	*/
#  define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
#  define G(u,v)          (-(v) + ( ((u)<1./3.) ? zero : \
                                    (1-6.75*(u)*(1-(u))*(1-(u))) ) )  
#endif

/*---------------------------------------------------------------------------*/
/*          In general you should not need to change anything below          */
/*---------------------------------------------------------------------------*/

/*  I rely on the compiler (gcc) to take care of most optimization.  
 *  The only place where I give the compiler help is with the nine-point 
 *  Laplacian formulas because apparently gcc does not do the obvious 
 *  thing.  I do not unroll the main loops in step() nor do I write explicit   
 *  pointer operations because gcc seems to handle all the indirections very 
 *  efficiently.
 */

/*---------------------------------------------------------------------------*/
/* 
 *  Defining the kinetics macros: U_KINETICS(u,uth) and V_KINETICS(u,v) 
 */

#if EXPLICIT   /* Explicit u-kinetics */

  #define U_KINETICS(u,uth) ( (u)+dt_o_eps*(u)*(one-(u))*((u)-(uth)) )

#else          /* Implicit u-kinetics */

/* The original (Physica 49D) scheme can be obtained by setting F2m =
   F2p = one */

  #define F1m(u,uth) ( dt_o_eps*(one-(u))*((u)-(uth)) )
  #define F1p(u,uth) ( dt_o_eps*(u)*((u)-(uth)) )
  #define F2m(u,uth) ( one + dt_o_2eps*((uth)-(u)*(u)) )
  #define F2p(u,uth) ( one + dt_o_2eps*(two*(u) -(uth)-(u)*(u)) )

  #define U_KINETICS(u,uth) (                                      \
          ( (u) < (uth) ) ?                                        \
          (u)/(one-F1m(u,uth)*F2m(u,uth) ) :                       \
          ((u)+F1p(u,uth)*F2p(u,uth))/(one+F1p(u,uth)*F2p(u,uth)) )

#endif

#define V_KINETICS(u,v) ( (v)+dt*G(u,v) )

/*---------------------------------------------------------------------------*/
/* 
 *  Defining the diffusion macros: U_DIFFUSION,  V_DIFFUSION,  
 *                                 ADD_TO_U_LAPLACIAN, ADD_TO_V_LAPLACIAN, 
 *                                 and ZERO_USED_LAPLACIANS. 
 */

#define U_DIFFUSION   ( NDu * u_lap[k1][i][j] ) 
#if V_DIFF_ON  
#  define V_DIFFUSION ( NDv * v_lap[k1][i][j] )
#else   
#  define V_DIFFUSION zero
#endif

#if NINEPOINT
/* Nine-point Laplacian formulas */
#  define ADD_TO_U_LAPLACIAN           \
     {precision stupid_cc = four*u[i][j]; \
      u_lap[k2][i][j] -= twenty*u[i][j];  \
      u_lap[k2][i+1][j] += stupid_cc;  \
      u_lap[k2][i-1][j] += stupid_cc;  \
      u_lap[k2][i][j+1] += stupid_cc;  \
      u_lap[k2][i][j-1] += stupid_cc;  \
      u_lap[k2][i+1][j+1] += u[i][j];  \
      u_lap[k2][i-1][j+1] += u[i][j];  \
      u_lap[k2][i+1][j-1] += u[i][j];  \
      u_lap[k2][i-1][j-1] += u[i][j];  \
     }
#  if V_DIFF_ON 
#    define ADD_TO_V_LAPLACIAN            \
       {precision stupid_cc = four*v[i][j]; \
        v_lap[k2][i][j] -= twenty*v[i][j];  \
        v_lap[k2][i+1][j] += stupid_cc;   \
        v_lap[k2][i-1][j] += stupid_cc;   \
        v_lap[k2][i][j+1] += stupid_cc;   \
        v_lap[k2][i][j-1] += stupid_cc;   \
        v_lap[k2][i+1][j+1] += v[i][j];   \
        v_lap[k2][i-1][j+1] += v[i][j];   \
        v_lap[k2][i+1][j-1] += v[i][j];   \
        v_lap[k2][i-1][j-1] += v[i][j];   \
       }
#  else
#    define ADD_TO_V_LAPLACIAN 
#  endif
#else
/* Five-point Laplacian formulas */
#  define ADD_TO_U_LAPLACIAN       \
     u_lap[k2][i][j] -= four*u[i][j]; \
     u_lap[k2][i+1][j] += u[i][j]; \
     u_lap[k2][i-1][j] += u[i][j]; \
     u_lap[k2][i][j+1] += u[i][j]; \
     u_lap[k2][i][j-1] += u[i][j];  
#  if V_DIFF_ON 
#    define ADD_TO_V_LAPLACIAN       \
       v_lap[k2][i][j] -= four*v[i][j]; \
       v_lap[k2][i+1][j] += v[i][j]; \
       v_lap[k2][i-1][j] += v[i][j]; \
       v_lap[k2][i][j+1] += v[i][j]; \
       v_lap[k2][i][j-1] += v[i][j]; 
#  else
#    define ADD_TO_V_LAPLACIAN 
#  endif
#endif


#if V_DIFF_ON 
#  define ZERO_USED_LAPLACIANS \
    u_lap[k1][i][j] = zero;       \
    v_lap[k1][i][j] = zero; 
#else  
#  define ZERO_USED_LAPLACIANS \
    u_lap[k1][i][j] = zero; 
#endif

/*---------------------------------------------------------------------------*/
