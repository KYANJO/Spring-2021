/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                                EZ-SPIRAL                                  */
/*                    A Code for Simulating Spiral Waves                     */
/*                                                                           */
/*           Copyright (C) 1992, 1994, 1997 Dwight Barkley                   */
/*                                          barkley@maths.warwick.ac.uk      */
/*                                                                           */
/*                                Version 2                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
 * $Revision: 2.4 $
 * $Date: 1997/11/13 23:22:00 $
 */
/*---------------------------------------------------------------------------*/

#include "ezspiral.h"

/* Global variables used throughout the EZ-Spiral Code */
precision  **u, **v,                    /* u and v fields */
           **u_lap[2], **v_lap[2],      /* u and v Laplacian fields */
           *tip_x, *tip_y,              /* spiral tip arrays */
           a_param, b_param,            /* model parameters */
           one_o_eps, length, Dv,       /* model parameters */
           ts, delta,                   /* numerical parameters */
           grid_h, dt,                  /* numerical parameters */
           NDu, NDv, one_o_a, b_o_a,    /* useful parameter combinations */
           dt_o_eps, dt_o_2eps;         /* useful parameter combinations */
int        ndim,                        /* # of grid points per direction */
           u_plot, v_plot,              /* flags indicating fields to plot */
           verbose,                     /* verbosity level */
           k1, k2, ntips;               /* useful integers */

/* Global variables and functions for this file only */
static int      nits, itip, iplot, iwrt, hist_x, hist_y;
static FILE     *history, *path;
static void     initialize();
static void     finish_up ();
static void     write_dat (int iter);
static void     write_fc  (char *filename);
static void     read_ic   (char *filename);
static void     next_line (FILE *filep);

/*---------------------------------------------------------------------------*/

void main()

     /* 
      *  Main routine. Takes nits time steps of the model. 
      */
{
  int iter;

  initialize();

  for(iter=0;iter<nits;iter++) {
    if(keyboard_chk())         break;
    if(iwrt && (iter%iwrt)==0) write_dat(iter);
    if(itip && (iter%itip)==0) tip_finder();
    if((iter%iplot)==0)        plot();
    step();
  }

  finish_up();
}
/*---------------------------------------------------------------------------*/

void initialize()

     /* 
      *  Opens files, reads data, initializes graphics and time stepping. 
      */
{
  double p_in;
  FILE *fp;
  int i;

  /* Read parameters from task.dat */
  if((fp=fopen("task.dat","r"))==NULL) {
    if(verbose) fprintf(stderr,"Cannot open task file: task.dat \n");
    exit(1);
  }
  else {
    fscanf(fp,"%lg",&p_in); next_line(fp); a_param=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); b_param=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); one_o_eps=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); length=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); Dv=p_in;
                            next_line(fp);
    fscanf(fp,"%d", &ndim); next_line(fp); 
    fscanf(fp,"%lg",&p_in); next_line(fp); ts=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); delta=p_in;
                            next_line(fp);
    fscanf(fp,"%d",&nits);  next_line(fp);
    fscanf(fp,"%d",&iplot); next_line(fp);
    fscanf(fp,"%d",&itip);  next_line(fp);
    fscanf(fp,"%d",&iwrt);  next_line(fp);
    fscanf(fp,"%d,%d",&hist_x,&hist_y); next_line(fp);
                            next_line(fp);
    fscanf(fp,"%d",&u_plot);next_line(fp);
    fscanf(fp,"%d",&v_plot);next_line(fp);
    fscanf(fp,"%d",&verbose);next_line(fp);
    fclose(fp);
  }

#define STABILITY_LIMIT (0.25*grid_h*grid_h)

#if NBC  
  /* Neumann  boundary conditions are used */
  grid_h   = length/(ndim-1);
#else    
  /* Periodic boundary conditions are used */
  grid_h   = length/ndim;
#endif
  dt       = ts*STABILITY_LIMIT;
  NDu      = dt/(grid_h*grid_h);
  NDv      = Dv*NDu;
#if NINEPOINT 
  /* There is a 6 in the denominator of the 9-point Laplacian formula */
  NDu      = NDu/6.;
  NDv      = NDv/6.;
#endif
  one_o_a   = 1.0/a_param; 
  b_o_a     = b_param/a_param; 
  dt_o_eps  = dt*one_o_eps; 
  dt_o_2eps = dt_o_eps/2.;

  if(verbose) {
    printf("\n\nModel Parameters: \n");
    printf("a     = %g\n", a_param);
    printf("b     = %g\n", b_param);
    printf("eps   = 1/%g = %g\n", one_o_eps, 1.0/one_o_eps);
    printf("L     = %g\n", length);
    printf("Dv    = %g\n", Dv);
    printf("\nNumerical Parameters: \n");
    printf("ndim = %-10d ts     = %-10g delta  = %-10g\n", 
	   ndim, ts, delta);
    printf("dt   = %-10g dt/eps = %-10g grid_h = %-10g\n", 
	   dt, dt_o_eps, grid_h);
    printf("dt/(grid_h*grid_h) = %-10g ", dt/(grid_h*grid_h) );
#if NINEPOINT
    printf("Using 9 point Laplacians.");
#endif
    printf("\n\nNumber of time steps = %d\n", nits);
    if(u_plot||v_plot) printf("    time steps per plot = %d\n", iplot);
    if(itip)   printf("    time steps per locating tip = %d\n", itip);
    if(iwrt)   printf("    time steps per write to data files= %d\n", iwrt);
    if(hist_x) printf("    history point = (%d, %d)\n", hist_x, hist_y ); 
  }

  /* Perform some tests */
  if( V_DIFF_ON && Dv==0.) {
    fprintf(stderr,"***** V_DIFF_ON is 1 and Dv == 0. ******\n");
    exit(1);
  }
  if(!V_DIFF_ON && Dv!=0.) {
    fprintf(stderr,"***** V_DIFF_ON is 0 and Dv != 0. ******\n");
    exit(1);
  }
  if(hist_x<0 || hist_x>ndim || hist_y<0 || hist_y>ndim) {
    fprintf(stderr,"***** history point out of range ******\n");
    exit(1);
  }
  if(ts > 1.0 ) {
    fprintf(stderr,"***** ts > 1 (the diffusion stability limit) ******\n");
    exit(1);
  }

  /*  Allocate memory for u[ndim+2][ndim+2], v[ndim+2][ndim+2], 
   *  u_lap[2][ndim+2][ndim+2], and if V_DIFF_ON, for v_lap[2][ndim+2][ndim+2]. 
   *  u and v could be of size [ndim][ndim]; however, because u_lap and v_lap 
   *  MUST be of size [2][ndim+2][ndim+2], I think it best that u and v also
   *  be of size [ndim+2][ndim+2].  The memory for each field is allocated 
   *  with the (ndim+2)*(ndim+2) locations contiguous.
   */
  u = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  v = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  u[0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  v[0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  for(i=1;i<=ndim+1;i++){
    u[i] = u[0]+i*(ndim+2);
    v[i] = v[0]+i*(ndim+2);
  }

  u_lap[0] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  u_lap[1] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  u_lap[0][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  u_lap[1][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  for(i=1;i<=ndim+1;i++){
    u_lap[0][i] = u_lap[0][0]+i*(ndim+2);
    u_lap[1][i] = u_lap[1][0]+i*(ndim+2);
  }

#if V_DIFF_ON
  v_lap[0] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  v_lap[1] = (precision **) malloc((unsigned)(ndim+2)*sizeof(precision*));
  v_lap[0][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  v_lap[1][0] = (precision *) calloc((ndim+2)*(ndim+2), sizeof(precision));
  for(i=1;i<=ndim+1;i++){
    v_lap[0][i] = v_lap[0][0]+i*(ndim+2);
    v_lap[1][i] = v_lap[1][0]+i*(ndim+2);
  }
#endif

  ntips = 0;
  if(itip) {  /* memory for tip path. */
    tip_x=(precision *)malloc((unsigned)(nits+1)*sizeof(precision));
    tip_y=(precision *)malloc((unsigned)(nits+1)*sizeof(precision));
  }
  if(iwrt) {  /* open history and tip-path files. */
    if(hist_x!=0) history = fopen("history.dat", "w");
    if(itip  !=0) path    = fopen("path.dat",    "w");
  }

  /* Read initial conditions and initialize time stepping and graphics */
  read_ic("ic.dat"); 
  step_ini();
  plot_ini();
}
/*---------------------------------------------------------------------------*/

void finish_up()

     /* 
      *  Writes final conditions, possibly writes last history and/or 
      *  tip-path data, and closes data files.
      */
{
  write_fc("fc.dat");
  if(iwrt) {
    if((nits%iwrt)==0) write_dat(nits); 
    if(hist_x!=0) fclose(history);
    if(itip  !=0) fclose(path);
  }
}
/*---------------------------------------------------------------------------*/

void write_dat(int iter)

     /* 
      *  Writes history and tip-path data. 
      */
{
  if(verbose>1) printf("iter = %d", iter);
  if(hist_x!=0) {
    fprintf(history, "%.5f %.5f %.5f \n", dt*(precision)iter, 
	    u[hist_x][hist_y], v[hist_x][hist_y]); 
    if(verbose>1) printf("; history (u,v) = %.5f, %.5f", 
			 u[hist_x][hist_y], v[hist_x][hist_y]); 
  }
  if(ntips>0) {
    fprintf(path, "%.5f %.5f %.5f \n", dt*(precision)iter, 
	    grid_h*(tip_x[ntips-1]-1), grid_h*(tip_y[ntips-1]-1) );
    if(verbose>1) printf("; tip (x,y) = %.5f %.5f", 
	    grid_h*(tip_x[ntips-1]-1), grid_h*(tip_y[ntips-1]-1) );
  }
  if(verbose>1) printf("\n"); 
}
/*---------------------------------------------------------------------------*/

void write_fc (char *filename)

     /* 
      *  Write final conditions to a file. 
      */
{
  FILE *fp;
  time_t tt1;
  int i,j;

  time(&tt1); 
  fp = fopen(filename, "w");
  fprintf(fp,"Model Parameters: a, b, 1/eps, L, Dv = %g, %g, %g, %g, %g\n",
	  a_param, b_param, one_o_eps, length, Dv);
  fprintf(fp,"Numerical Parameters: ndim, ts, delta = %d, %g, %g \n", 
	  ndim, ts, delta);
  fprintf(fp,"File written: %s",ctime(&tt1));
  fprintf(fp,"Comments:\n");
  fprintf(fp,"Values of u and v follow\n");
  for(i=1;i<=ndim;i++) {
    for(j=1;j<=ndim;j++) {
      fprintf (fp, "%g %g\n", (float)u[i][j], (float)v[i][j]);
    }
  }
  fclose(fp);
}
/*---------------------------------------------------------------------------*/

void read_ic (char *filename)
     
     /* 
      *  Reads initial condition file if it exists, 
      *  otherwise generates new initial condition. 
      */
{
  precision **u_tmp, **v_tmp, ratio;
  double u_in, v_in;
  int ndim_ic, i_tmp, j_tmp, dummy, i, j;
  FILE *fp;
  
  if((fp=fopen(filename,"r"))!=NULL) { /* if file can be opened then read it */
    next_line(fp);
    while( (dummy=getc(fp)) != '='); fscanf (fp, "%d", &ndim_ic); 
    next_line(fp); 
    next_line(fp);
    next_line(fp);
    next_line(fp);
    if(verbose) printf("Reading ic.dat with ndim = %d... \n\n", ndim_ic);

    /* Allocate temporary memory (appropriate for in ndim_ic) */
    u_tmp =((precision **) malloc((unsigned)(ndim_ic)*sizeof(precision*)))-1;
    v_tmp =((precision **) malloc((unsigned)(ndim_ic)*sizeof(precision*)))-1;
    for(i=1;i<=ndim_ic;i++){
      u_tmp[i]=((precision *) malloc((unsigned)(ndim_ic)*sizeof(precision)))-1;
      v_tmp[i]=((precision *) malloc((unsigned)(ndim_ic)*sizeof(precision)))-1;
    }

    /* Read into temporary memory */
    for(i=1;i<=ndim_ic;i++) {
      for(j=1;j<=ndim_ic;j++) {
	fscanf (fp, "%lg %lg\n", &u_in, &v_in);
	u_tmp[i][j] = u_in; v_tmp[i][j] = v_in;
      }
    }
    
    /* Copy into u and v */
    ratio = (ndim_ic-1.0)/(ndim-1.0);
    for(i=1;i<=ndim;i++) {
      i_tmp = 1+(int)(ratio*(i-1));
      for(j=1;j<=ndim;j++) {
	j_tmp = 1+(int)(ratio*(j-1));
	u[i][j] = u_tmp[i_tmp][j_tmp];
	v[i][j] = v_tmp[i_tmp][j_tmp];
      }
    }
    
    /* Free temporary memory */
    for(i=ndim_ic;i>=1;i--) free((char*)(u_tmp[i]+1));
    for(i=ndim_ic;i>=1;i--) free((char*)(v_tmp[i]+1));
    free((char*) (u_tmp+1));
    free((char*) (v_tmp+1));
  }
  else {     
    /* If file cannot be opened then generate new initial condition */
    if(verbose) printf("No file ic.dat. Generating initial condition. \n\n");
    for(i=1;i<=ndim;i++) {
      for(j=1;j<=ndim;j++) {
	if(i<(10+ndim/2)) u[i][j] = 0.; else u[i][j] = 1.0; 
	if(j<ndim/2)      v[i][j] = 0.; else v[i][j] = a_param/2.; 
      }
    }
  }
}
/*---------------------------------------------------------------------------*/

void next_line(FILE *filep)
     
     /* 
      *  Skips to next input line. 
      */
{
  int dummy;
  while( (dummy=getc(filep)) != '\n');
}
/*---------------------------------------------------------------------------*/
