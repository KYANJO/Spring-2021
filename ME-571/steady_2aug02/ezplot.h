/* ------------------------------------------------------------------------- */
/* ezplot.h 
 *
 * Copyright (C) 1998-2002 Dwight Barkley                 
 *
 * ------------------------------------------------------------------------- */

#ifndef _EZPLOT2D_H_
#define _EZPLOT2D_H_

#ifndef _XFUNCPROTOBEGIN
#ifdef __cplusplus		/* for C++ V2.0 */
#define _XFUNCPROTOBEGIN extern "C" {	/* do not leave open across includes */
#define _XFUNCPROTOEND }
#else
#define _XFUNCPROTOBEGIN
#define _XFUNCPROTOEND
#endif
#endif /* _XFUNCPROTOBEGIN */

_XFUNCPROTOBEGIN

/* =========================== Public functions ============================ */
void EZplot_ini (double xmin, double xmax, double umin, double umax);

/* Must be called before any other EZplot functions.  Initializes window and
 * sets plot limits: xmin <= x <= xmax, umin <= u <= umax. */

/* ========================================================================= */

void EZplot_simple1 (double *u, int npoints);

/* Simple plotting function.  Plots the curve (j, u[j]) for j=0 to
 * npoints-1.  EZplot_ini() should be called with xmin=0. and
 * xmax=(double)npoints-1.  umin and umax should be set appropriately. */

/* ========================================================================= */

void EZplot_line1 (double *x, double *u, int npoints, char color);

/* Plot the curve (x[j], u[j]) where x and u are array with npoints
 * elements.  Note, the plot will not appear until EZplot_display() is
 * called. The character color sets the color of the curve.  See note on
 * colors below. */

/* ========================================================================= */

void EZplot_points1 (double *x, double *u, int npoints, char color);

/* Plot the curve (x[j], u[j]) where x and u are array with npoints
 * elements.  Note, the plot will not appear until EZplot_display() is
 * called. The character color sets the color of the curve.  See note on
 * colors below. */

/* ========================================================================= */

void EZplot_clear (char color);

/* Clear the plot region to specified color.  Note, the plot will not be
 * cleared until EZplot_display() is called. See note on colors below. */

/* ========================================================================= */

void EZplot_display (void);

/* Display the plot in the X-window. */

/* ========================================================================= */

int EZevent (void);

/* Simple event function that responds to characters typed on the keyboard
 * (the pointer must be in the X-window).  The functions are: 
 *    "p" pauses the simulation
 *    space resumes the simulation
 *    "q" causes EZevent to return a value 1
 *    escape exits the program immediately
 * */

/* ========================================================================= */

#ifdef MY_KEY_EVENT_HANDLE
int My_Key_Event_Handle (char Key_character);
#endif

void EZ_Pause (void);
void EZ_Resume (void);
int  EZ_Quit (void);
void EZ_Kill (void);

/* ========================================================================= */

/* Colors: The following colors are allowed with EZplot and correspond to the
 * following char specification in function calls.
 *
 * color   char
 * -----   ----
 * red     'r'
 * green   'g'
 * blue    'b'
 * white   'w'
 * black   ' '
 */

/* ============================= 2D functions ============================== */

void EZplot2d_ini (double min, double max, double umin, double umax);

/* Must be called before any other 2d functions.  Initializes window and
 * sets plot limits: min <= x <= max, min <= y <= max.  Sets contour
 * levels based on umin and umax which must satisfy: umin<0 and umax>0. 
 *
 * Contour colors are (from smallest to largest) green, blue, red white
 * with the blue/red border being at u=0. */

/* ========================================================================= */

void EZplot_simple2 (double *u, int nj);

/* Simple 2D contour plotting function.  Plots u(j,k) for j,k =0 to
 * nj-1.  EZplot2d_ini() should be called with min=0 and max=(double)nj-1.
 * umin and umax should be set appropriately. */

/* ========================================================================= */

void EZplot_polar (double *u, double dr, double dtheta, int nr, int ntheta);

/* Simple 2D contour plotting function in polar coordinates.
 * EZplot2d_ini() should be called with min=umin=-radius and
 * max=umax=radius */

/* ========================================================================= */

void EZplot_polar_line1 (double *r, double *theta, int nj, char color);

/* Simple curve plotting function in polar coordinates.
 * EZplot2d_ini() should be called with min=umin=-radius and
 * max=umax=radius */

/* ========================================================================= */


_XFUNCPROTOEND
#endif /* _EZPLOT2D_H_ */
