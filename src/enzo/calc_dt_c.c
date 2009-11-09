#ifdef CONFIG_PFLOAT_16
/* calc_dt.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include <stdio.h>
#include <math.h>

#define ENZO_PYTHON_IMPORTED
#include "macros_and_parameters.h"
#undef ENZO_PYTHON_IMPORTED

/* Table of constant values */

static int c__9 = 9;
static int c__1 = 1;
static int c__4 = 4;
static int c__3 = 3;

/* ======================================================================= */
/* ////////////////////////  SUBROUTINE CALC_DT  \\\\\\\\\\\\\\\\\\\\\\\\c */

/* Subroutine */ 
int calc_dt_c(int *rank, int *idim, int *jdim, 
	    int *kdim, int *i1, int *i2, int *j1, int *j2, 
	    int *k1, int *k2, int *ihydro, float *c2, FLOAT *dx, FLOAT *dy, 
	    FLOAT *dz, float *vgx, float *vgy, float *vgz, float *gamma, 
	    int * ipfree, float *aye, float *d__, float *p, float *u, float *v, 
	    float *w, float *dt, float *dtviscous)
{
  /* System generated locals */
  int d_dim1, d_dim2, d_offset, p_dim1, p_dim2, p_offset, u_dim1, 
    u_dim2, u_offset, v_dim1, v_dim2, v_offset, w_dim1, w_dim2, 
    w_offset, i__1, i__2, i__3;
  float r__1, r__2, r__3, r__4, r__5;

  /* Local variables */
  static int i__, j, k;
  static float cs, dt1;

/*  COMPUTES TIME STEP FOR NEXT CYCLE */

/*     written by: Greg Bryan */
/*     date:       February, 1996 */
/*     modified1: */

/*  PURPOSE:  Computes the new timestep using the Courant condition. */
/*            (For rank < 3, the unused fields and cell widths may be */
/*             null) */

/*  INPUTS: */
/*    rank    - rank of fields */
/*    i,j,dim - declared dimensions of fields */
/*    i,j,k1  - start index of active region in fields (0 based) */
/*    i,j,k2  - end index of active region in fields (0 based) */
/*    ihydro  - Hydro method (2 - Zeus), used for viscosity computation */
/*    C2      - coefficient of quadratic artificial viscosity */
/*    dx,y,z  - cell widths along each dimension */
/*    vgx,y,z - grid bulk velocity */
/*    gamma   - ratio of specific heats */
/*    ipfree  - pressure free flag (1 = on, 0 = off) */
/*    aye     - expansion factor (or 1 if not using comoving coordinates) */
/*    d,p     - density and pressure fields */
/*    u,v,w   - velocity fields (x,y,z) */
/*    dtviscous - viscous time for stability (if used) */

/*  OUTPUTS: */
/*    dt      - minimum allowed dt (without Courant safety factor) */

/*  LOCALS: */

/* ----------------------------------------------------------------------- */


/*     Arguments */


/*     Locals */


/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////// */
/* ======================================================================= */

/*     Set initial timestep to a large number */

    /* Parameter adjustments */
    --dx;
    --dy;
    w_dim1 = *idim;
    w_dim2 = *jdim;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idim;
    v_dim2 = *jdim;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    u_dim1 = *idim;
    u_dim2 = *jdim;
    u_offset = 1 + u_dim1 * (1 + u_dim2);
    u -= u_offset;
    p_dim1 = *idim;
    p_dim2 = *jdim;
    p_offset = 1 + p_dim1 * (1 + p_dim2);
    p -= p_offset;
    d_dim1 = *idim;
    d_dim2 = *jdim;
    d_offset = 1 + d_dim1 * (1 + d_dim2);
    d__ -= d_offset;
    --dz;

    /* Function Body */
    *dt = 1e20f;

/*     one-dimensional version */

    if (*rank == 1) {
	i__1 = *i2 + 1;
	for (i__ = *i1 + 1; i__ <= i__1; ++i__) {
/* Computing max */
	    r__1 = sqrt(*gamma * p[i__ + (p_dim2 + 1) * p_dim1] / d__[i__ + (
		    d_dim2 + 1) * d_dim1]);
	    cs = max(r__1,1e-20f);
	    if (*ipfree == 1) {
		cs = 1e-20f;
	    }
/* Computing min */
	    r__2 = *dt, r__3 = dx[i__] * *aye / (cs + (r__1 = u[i__ + (u_dim2 
		    + 1) * u_dim1] - *vgx, dabs(r__1)));
	    *dt = min(r__2,r__3);
	}
	if (*ihydro == 2) {
	    i__1 = *i2 + 1;
	    for (i__ = *i1 + 1; i__ <= i__1; ++i__) {
/* Computing min */
/* Computing max */
		r__3 = -u[i__ + 1 + (u_dim2 + 1) * u_dim1] + u[i__ + (u_dim2 
			+ 1) * u_dim1];
		r__1 = *dtviscous, r__2 = dx[i__] * *aye / (*c2 * 4.f * max(
			r__3,1e-20f));
		*dtviscous = min(r__1,r__2);
	    }
	}
    }

/*     two-dimensional version */

    if (*rank == 2) {
	i__1 = *j2 + 1;
	for (j = *j1 + 1; j <= i__1; ++j) {
	    i__2 = *i2 + 1;
	    for (i__ = *i1 + 1; i__ <= i__2; ++i__) {
/* Computing max */
		r__1 = sqrt(*gamma * p[i__ + (j + p_dim2) * p_dim1] / d__[i__ 
			+ (j + d_dim2) * d_dim1]);
		cs = max(r__1,1e-20f);
		if (*ipfree == 1) {
		    cs = 1e-20f;
		}
/* Computing min */
		r__3 = *dt, r__4 = dx[i__] * *aye / (cs + (r__1 = u[i__ + (j 
			+ u_dim2) * u_dim1] - *vgx, dabs(r__1))), r__3 = min(
			r__3,r__4), r__4 = dy[j] * *aye / (cs + (r__2 = v[i__ 
			+ (j + v_dim2) * v_dim1] - *vgy, dabs(r__2)));
		*dt = min(r__3,r__4);
	    }

	    if (*ihydro == 2) {
		i__2 = *i2 + 1;
		for (i__ = *i1 + 1; i__ <= i__2; ++i__) {
/* Computing min */
/* Computing max */
		    r__3 = -u[i__ + 1 + (j + u_dim2) * u_dim1] + u[i__ + (j + 
			    u_dim2) * u_dim1];
/* Computing max */
		    r__4 = -v[i__ + (j + 1 + v_dim2) * v_dim1] + v[i__ + (j + 
			    v_dim2) * v_dim1];
		    r__1 = *dtviscous, r__2 = dx[i__] * *aye / (*c2 * 4.f * 
			    max(r__3,1e-20f)), r__1 = min(r__1,r__2), r__2 = 
			    dy[j] * *aye / (*c2 * 4.f * max(r__4,1e-20f));
		    *dtviscous = min(r__1,r__2);
		}
	    }
	}
    }

/*     three-dimensional version */

    if (*rank == 3) {
	i__1 = *k2 + 1;
	for (k = *k1 + 1; k <= i__1; ++k) {
	    i__2 = *j2 + 1;
	    for (j = *j1 + 1; j <= i__2; ++j) {
		i__3 = *i2 + 1;
		for (i__ = *i1 + 1; i__ <= i__3; ++i__) {
		    if (d__[i__ + (j + k * d_dim2) * d_dim1] != d__[i__ + (j 
			    + k * d_dim2) * d_dim1] || p[i__ + (j + k * 
			    p_dim2) * p_dim1] != p[i__ + (j + k * p_dim2) * 
			    p_dim1]) {

		      fprintf(stdout, "calc_dt d,p,i,j,k= %g %g %i %i %i ",
			      d__[i__ + (j + k * d_dim2) * d_dim1],
			        p[i__ + (j + k * d_dim2) * d_dim1], i__,j,k);
			      
		    }
/* Computing max */
		    r__1 = sqrt(*gamma * p[i__ + (j + k * p_dim2) * p_dim1] / 
			    d__[i__ + (j + k * d_dim2) * d_dim1]);
		    cs = max(r__1,1e-20f);
		    if (*ipfree == 1) {
			cs = 1e-20f;
		    }
/* Computing min */
		    r__4 = dx[i__] * *aye / (cs + (r__1 = u[i__ + (j + k * 
			    u_dim2) * u_dim1] - *vgx, dabs(r__1))), r__5 = dy[
			    j] * *aye / (cs + (r__2 = v[i__ + (j + k * v_dim2)
			     * v_dim1] - *vgy, dabs(r__2))), r__4 = min(r__4,
			    r__5), r__5 = dz[k] * *aye / (cs + (r__3 = w[i__ 
			    + (j + k * w_dim2) * w_dim1] - *vgz, dabs(r__3)));
		    dt1 = min(r__4,r__5);
		    *dt = min(*dt,dt1);
/*                  if (dt1 .lt. 1.0e-5) write(6,1000) dt1,d(i,j,k), */
/*     &                  p(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k) */
/* L1000: */
		}

		if (*ihydro == 2) {
		    i__3 = *i2 + 1;
		    for (i__ = *i1 + 1; i__ <= i__3; ++i__) {
/* Computing min */
/* Computing max */
			r__3 = -u[i__ + 1 + (j + k * u_dim2) * u_dim1] + u[
				i__ + (j + k * u_dim2) * u_dim1];
/* Computing max */
			r__4 = -v[i__ + (j + 1 + k * v_dim2) * v_dim1] + v[
				i__ + (j + k * v_dim2) * v_dim1];
/* Computing max */
			r__5 = -w[i__ + (j + (k + 1) * w_dim2) * w_dim1] + w[
				i__ + (j + k * w_dim2) * w_dim1];
			r__1 = dx[i__] * *aye / (*c2 * 4.f * max(r__3,1e-20f)
				), r__2 = dy[j] * *aye / (*c2 * 4.f * max(
				r__4,1e-20f)), r__1 = min(r__1,r__2), r__2 = 
				dz[k] * *aye / (*c2 * 4.f * max(r__5,1e-20f))
				;
			dt1 = min(r__1,r__2);
			*dtviscous = min(*dtviscous,dt1);
		    }
		}
	    }
	}
    }

    return 0;
} /* calc_dt__ */

#endif
