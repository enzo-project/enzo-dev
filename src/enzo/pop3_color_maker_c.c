#ifdef CONFIG_PFLOAT_16
/* pop3_color_maker.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

#define ENZO_PYTHON_IMPORTED
#include "macros_and_parameters.h"
#undef ENZO_PYTHON_IMPORTED

/* ======================================================================= */
/* /////////////////////  SUBROUTINE POP3_COLOR_MAKER \\\\\\\\\\\\\\\\\\\\ */

/* Subroutine */ int pop3_color_maker_c(integer *nx, integer *ny, integer *nz,
	 real *d__, real *dm, real *u, real *v, real *w, real *dt, real *r__, 
	real *dx, FLOAT *t, real *z__, integer *procnum, real *d1, real *x1, 
	real *v1, real *t1, integer *nmax,
    FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
    integer *ibuff, integer *imethod, real *odthresh, integer *level,
    integer *np, FLOAT *xp, FLOAT *yp, FLOAT *zp, real *up, real *vp, 
	real *wp, real *mp, real *tcp, integer *type__, integer *ctype)
{
    /* System generated locals */
    integer d_dim1, d_dim2, d_offset, dm_dim1, dm_dim2, dm_offset, u_dim1, 
	    u_dim2, u_offset, v_dim1, v_dim2, v_offset, w_dim1, w_dim2, 
	    w_offset, r_dim1, r_dim2, r_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static real m1;
    static integer ii;


/*  CREATES PRIMORDIAL STAR PARTICLES */

/*  written by: Chris Loken */
/*  date:       3 March 1997 */
/*  modified1: 2 March 1999 by Brian O''Shea */
/*    2 inputs were added: odthresh and masseff, and the star creation */
/*    code was implemented */
/*  modified2: 18 May 1999 by Brian O''Shea */
/*    1 input added: smthresh, and star particle selection code added */
/*  modified3: 30 August 2005 by John Wise */
/*    2 inputs added: type and ctype to distinugish between particles now */
/*    that we have multiple star particle types, and added BWO's fix on */
/*    the runaway star particle phenomenon by averaging the velocities */
/*    over 5 cells */
/*  modified4: 30 August 2005 by John Wise */
/*    modified star_maker2 for single primordial star formation */
/*  modified5: June, 2009 by Matthew Turk */
/*    modified pop3_maker for much simpler, overdensity star formation */

/*  INPUTS: */

/*    d     - density field */
/*    dm    - dark matter field */
/*    h2d   - molecular hydrogen field */
/*    u,v,w - velocity fields */
/*    r     - refinement field (non-zero if zone is further refined) */
/*    dt    - current timestep */
/*    dx    - zone size (code units) */
/*    t     - current time */
/*    z     - current redshift */
/*    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units */
/*    nx,ny,nz - dimensions of field arrays */
/*    ibuff    - number of buffer zones at each end of grid */
/*    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS) */
/*    odthresh - overdensity threshold (some number * avg. density) */
/*    starmass - primordial star mass (solar masses) */
/*    h2crit   - critical value for primordial star formation */
/*    metalcrit- critical value for primordial star formation */
/*    level - current level of refinement */
/*    procnum - processor number (for output) */

/*  OUTPUTS: */

/*    np   - number of particles created */
/*    x/y/z start - starting position of grid origin */
/*    xp,yp,zp - positions of created particles */
/*    up,vp,wp - velocities of created particles */
/*    mp       - mass of new particles */
/*    tcp      - creation time of particle */
/*    nmax     - particle array size specified by calling routine */
/*    type     - particle types */
/*    ctype    - current type to set the particles to */
/*    justburn     - time-weighted mass of star formation (code units) */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

/*  Arguments */


/*  Locals: */


    /* Parameter adjustments */
    r_dim1 = *nx;
    r_dim2 = *ny;
    r_offset = 1 + r_dim1 * (1 + r_dim2 * 1);
    r__ -= r_offset;
    w_dim1 = *nx;
    w_dim2 = *ny;
    w_offset = 1 + w_dim1 * (1 + w_dim2 * 1);
    w -= w_offset;
    v_dim1 = *nx;
    v_dim2 = *ny;
    v_offset = 1 + v_dim1 * (1 + v_dim2 * 1);
    v -= v_offset;
    u_dim1 = *nx;
    u_dim2 = *ny;
    u_offset = 1 + u_dim1 * (1 + u_dim2 * 1);
    u -= u_offset;
    dm_dim1 = *nx;
    dm_dim2 = *ny;
    dm_offset = 1 + dm_dim1 * (1 + dm_dim2 * 1);
    dm -= dm_offset;
    d_dim1 = *nx;
    d_dim2 = *ny;
    d_offset = 1 + d_dim1 * (1 + d_dim2 * 1);
    d__ -= d_offset;
    --type__;
    --tcp;
    --mp;
    --wp;
    --vp;
    --up;
    --zp;
    --yp;
    --xp;

    /* Function Body */
    ii = *np;

/*     calculate how many solar masses are in d1*x1^3 */
/* Computing 3rd power */
    d__1 = (doublereal) (*x1 * *dx);
    m1 = d__1 * (d__1 * d__1) * (doublereal) (*d1) / 1.989e33;

/*  for each zone, : a primordial "star" particle is created if answers */
/*  to all the following questions are affirmative: */

/*    is this the finest level of refinement ? */
/*    is the density greater than a critical density ? */
/*    is the molecular hydrogen fraction greater than critical value ? */
/*    is the metallicity less than a critical value ? */
/*    is the flow convergent ? */
/*    is the cooling time less than a dynamical time ? */

    i__1 = *nz - *ibuff;
    for (k = *ibuff + 1; k <= i__1; ++k) {
	i__2 = *ny - *ibuff;
	for (j = *ibuff + 1; j <= i__2; ++j) {
	    i__3 = *nx - *ibuff;
	    for (i__ = *ibuff + 1; i__ <= i__3; ++i__) {

/*              1) finest level of refinement? */

/*               write(0,*) d(i,j,k) */
		if (r__[i__ + (j + k * r_dim2) * r_dim1] != 0.f) {
		    goto L10;
		}

/*              2) density greater than threshold */

		if (d__[i__ + (j + k * d_dim2) * d_dim1] * *d1 < *odthresh) {
		    goto L10;
		}

/*              3) divergence negative */
/*                 (the first calculation is face centered for ZEUS, */
/*                  the second is cell-centered for PPM) */


/*              Create a star particle */

		++ii;
		type__[ii] = -(*ctype);
		tcp[ii] = *t;
		xp[ii] = *xstart + ((FLOAT) i__ - .5f) * *dx;
		yp[ii] = *ystart + ((FLOAT) j - .5f) * *dx;
		zp[ii] = *zstart + ((FLOAT) k - .5f) * *dx;

/*              Star velocities averaged over multiple cells to */
/*              avoid "runaway star particle" phenomenon */
/*              imethod = 2 is zeus, otherwise PPM */

		if (*imethod == 2) {
/*                  up(ii) = 0.5*(u(i,j,k)+u(i+1,j,k)) */
/*                  vp(ii) = 0.5*(v(i,j,k)+v(i,j+1,k)) */
/*                  wp(ii) = 0.5*(w(i,j,k)+w(i,j,k+1)) */
		    up[ii] = ((u[i__ + (j + k * u_dim2) * u_dim1] + u[i__ + 1 
			    + (j + k * u_dim2) * u_dim1]) * .5f * d__[i__ + (
			    j + k * d_dim2) * d_dim1] + (u[i__ - 1 + (j + k * 
			    u_dim2) * u_dim1] + u[i__ + (j + k * u_dim2) * 
			    u_dim1]) * .5f * d__[i__ - 1 + (j + k * d_dim2) * 
			    d_dim1] + (u[i__ + 1 + (j + k * u_dim2) * u_dim1] 
			    + u[i__ + 2 + (j + k * u_dim2) * u_dim1]) * .5f * 
			    d__[i__ + 1 + (j + k * d_dim2) * d_dim1] + (u[i__ 
			    + (j + 1 + k * u_dim2) * u_dim1] + u[i__ + 1 + (j 
			    + 1 + k * u_dim2) * u_dim1]) * .5f * d__[i__ + (j 
			    + 1 + k * d_dim2) * d_dim1] + (u[i__ + (j - 1 + k 
			    * u_dim2) * u_dim1] + u[i__ + 1 + (j - 1 + k * 
			    u_dim2) * u_dim1]) * .5f * d__[i__ + (j - 1 + k * 
			    d_dim2) * d_dim1] + (u[i__ + (j + (k + 1) * 
			    u_dim2) * u_dim1] + u[i__ + 1 + (j + (k + 1) * 
			    u_dim2) * u_dim1]) * .5f * d__[i__ + (j + (k + 1) 
			    * d_dim2) * d_dim1] + (u[i__ + (j + (k - 1) * 
			    u_dim2) * u_dim1] + u[i__ + 1 + (j + (k - 1) * 
			    u_dim2) * u_dim1]) * .5f * d__[i__ + (j + (k - 1) 
			    * d_dim2) * d_dim1]) / (d__[i__ + (j + k * d_dim2)
			     * d_dim1] + d__[i__ - 1 + (j + k * d_dim2) * 
			    d_dim1] + d__[i__ + 1 + (j + k * d_dim2) * d_dim1]
			     + d__[i__ + (j - 1 + k * d_dim2) * d_dim1] + d__[
			    i__ + (j + 1 + k * d_dim2) * d_dim1] + d__[i__ + (
			    j + (k - 1) * d_dim2) * d_dim1] + d__[i__ + (j + (
			    k + 1) * d_dim2) * d_dim1]);
		    vp[ii] = ((v[i__ + (j + k * v_dim2) * v_dim1] + v[i__ + 1 
			    + (j + k * v_dim2) * v_dim1]) * .5f * d__[i__ + (
			    j + k * d_dim2) * d_dim1] + (v[i__ - 1 + (j + k * 
			    v_dim2) * v_dim1] + v[i__ + (j + k * v_dim2) * 
			    v_dim1]) * .5f * d__[i__ - 1 + (j + k * d_dim2) * 
			    d_dim1] + (v[i__ + 1 + (j + k * v_dim2) * v_dim1] 
			    + v[i__ + 2 + (j + k * v_dim2) * v_dim1]) * .5f * 
			    d__[i__ + 1 + (j + k * d_dim2) * d_dim1] + (v[i__ 
			    + (j + 1 + k * v_dim2) * v_dim1] + v[i__ + 1 + (j 
			    + 1 + k * v_dim2) * v_dim1]) * .5f * d__[i__ + (j 
			    + 1 + k * d_dim2) * d_dim1] + (v[i__ + (j - 1 + k 
			    * v_dim2) * v_dim1] + v[i__ + 1 + (j - 1 + k * 
			    v_dim2) * v_dim1]) * .5f * d__[i__ + (j - 1 + k * 
			    d_dim2) * d_dim1] + (v[i__ + (j + (k + 1) * 
			    v_dim2) * v_dim1] + v[i__ + 1 + (j + (k + 1) * 
			    v_dim2) * v_dim1]) * .5f * d__[i__ + (j + (k + 1) 
			    * d_dim2) * d_dim1] + (v[i__ + (j + (k - 1) * 
			    v_dim2) * v_dim1] + v[i__ + 1 + (j + (k - 1) * 
			    v_dim2) * v_dim1]) * .5f * d__[i__ + (j + (k - 1) 
			    * d_dim2) * d_dim1]) / (d__[i__ + (j + k * d_dim2)
			     * d_dim1] + d__[i__ - 1 + (j + k * d_dim2) * 
			    d_dim1] + d__[i__ + 1 + (j + k * d_dim2) * d_dim1]
			     + d__[i__ + (j - 1 + k * d_dim2) * d_dim1] + d__[
			    i__ + (j + 1 + k * d_dim2) * d_dim1] + d__[i__ + (
			    j + (k - 1) * d_dim2) * d_dim1] + d__[i__ + (j + (
			    k + 1) * d_dim2) * d_dim1]);
		    wp[ii] = ((w[i__ + (j + k * w_dim2) * w_dim1] + w[i__ + 1 
			    + (j + k * w_dim2) * w_dim1]) * .5f * d__[i__ + (
			    j + k * d_dim2) * d_dim1] + (w[i__ - 1 + (j + k * 
			    w_dim2) * w_dim1] + w[i__ + (j + k * w_dim2) * 
			    w_dim1]) * .5f * d__[i__ - 1 + (j + k * d_dim2) * 
			    d_dim1] + (w[i__ + 1 + (j + k * w_dim2) * w_dim1] 
			    + w[i__ + 2 + (j + k * w_dim2) * w_dim1]) * .5f * 
			    d__[i__ + 1 + (j + k * d_dim2) * d_dim1] + (w[i__ 
			    + (j + 1 + k * w_dim2) * w_dim1] + w[i__ + 1 + (j 
			    + 1 + k * w_dim2) * w_dim1]) * .5f * d__[i__ + (j 
			    + 1 + k * d_dim2) * d_dim1] + (w[i__ + (j - 1 + k 
			    * w_dim2) * w_dim1] + w[i__ + 1 + (j - 1 + k * 
			    w_dim2) * w_dim1]) * .5f * d__[i__ + (j - 1 + k * 
			    d_dim2) * d_dim1] + (w[i__ + (j + (k + 1) * 
			    w_dim2) * w_dim1] + w[i__ + 1 + (j + (k + 1) * 
			    w_dim2) * w_dim1]) * .5f * d__[i__ + (j + (k + 1) 
			    * d_dim2) * d_dim1] + (w[i__ + (j + (k - 1) * 
			    w_dim2) * w_dim1] + w[i__ + 1 + (j + (k - 1) * 
			    w_dim2) * w_dim1]) * .5f * d__[i__ + (j + (k - 1) 
			    * d_dim2) * d_dim1]) / (d__[i__ + (j + k * d_dim2)
			     * d_dim1] + d__[i__ - 1 + (j + k * d_dim2) * 
			    d_dim1] + d__[i__ + 1 + (j + k * d_dim2) * d_dim1]
			     + d__[i__ + (j - 1 + k * d_dim2) * d_dim1] + d__[
			    i__ + (j + 1 + k * d_dim2) * d_dim1] + d__[i__ + (
			    j + (k - 1) * d_dim2) * d_dim1] + d__[i__ + (j + (
			    k + 1) * d_dim2) * d_dim1]);
		} else {
/*                  up(ii) = u(i,j,k) */
/*                  vp(ii) = v(i,j,k) */
/*                  wp(ii) = w(i,j,k) */
		    up[ii] = (u[i__ + (j + k * u_dim2) * u_dim1] * d__[i__ + (
			    j + k * d_dim2) * d_dim1] + u[i__ - 1 + (j + k * 
			    u_dim2) * u_dim1] * d__[i__ - 1 + (j + k * d_dim2)
			     * d_dim1] + u[i__ + 1 + (j + k * u_dim2) * 
			    u_dim1] * d__[i__ + 1 + (j + k * d_dim2) * d_dim1]
			     + u[i__ + (j - 1 + k * u_dim2) * u_dim1] * d__[
			    i__ + (j - 1 + k * d_dim2) * d_dim1] + u[i__ + (j 
			    + 1 + k * u_dim2) * u_dim1] * d__[i__ + (j + 1 + 
			    k * d_dim2) * d_dim1] + u[i__ + (j + (k - 1) * 
			    u_dim2) * u_dim1] * d__[i__ + (j + (k - 1) * 
			    d_dim2) * d_dim1] + u[i__ + (j + (k + 1) * u_dim2)
			     * u_dim1] * d__[i__ + (j + (k + 1) * d_dim2) * 
			    d_dim1]) / (d__[i__ + (j + k * d_dim2) * d_dim1] 
			    + d__[i__ - 1 + (j + k * d_dim2) * d_dim1] + d__[
			    i__ + 1 + (j + k * d_dim2) * d_dim1] + d__[i__ + (
			    j - 1 + k * d_dim2) * d_dim1] + d__[i__ + (j + 1 
			    + k * d_dim2) * d_dim1] + d__[i__ + (j + (k - 1) *
			     d_dim2) * d_dim1] + d__[i__ + (j + (k + 1) * 
			    d_dim2) * d_dim1]);
		    vp[ii] = (v[i__ + (j + k * v_dim2) * v_dim1] * d__[i__ + (
			    j + k * d_dim2) * d_dim1] + v[i__ - 1 + (j + k * 
			    v_dim2) * v_dim1] * d__[i__ - 1 + (j + k * d_dim2)
			     * d_dim1] + v[i__ + 1 + (j + k * v_dim2) * 
			    v_dim1] * d__[i__ + 1 + (j + k * d_dim2) * d_dim1]
			     + v[i__ + (j - 1 + k * v_dim2) * v_dim1] * d__[
			    i__ + (j - 1 + k * d_dim2) * d_dim1] + v[i__ + (j 
			    + 1 + k * v_dim2) * v_dim1] * d__[i__ + (j + 1 + 
			    k * d_dim2) * d_dim1] + v[i__ + (j + (k - 1) * 
			    v_dim2) * v_dim1] * d__[i__ + (j + (k - 1) * 
			    d_dim2) * d_dim1] + v[i__ + (j + (k + 1) * v_dim2)
			     * v_dim1] * d__[i__ + (j + (k + 1) * d_dim2) * 
			    d_dim1]) / (d__[i__ + (j + k * d_dim2) * d_dim1] 
			    + d__[i__ - 1 + (j + k * d_dim2) * d_dim1] + d__[
			    i__ + 1 + (j + k * d_dim2) * d_dim1] + d__[i__ + (
			    j - 1 + k * d_dim2) * d_dim1] + d__[i__ + (j + 1 
			    + k * d_dim2) * d_dim1] + d__[i__ + (j + (k - 1) *
			     d_dim2) * d_dim1] + d__[i__ + (j + (k + 1) * 
			    d_dim2) * d_dim1]);
		    wp[ii] = (w[i__ + (j + k * w_dim2) * w_dim1] * d__[i__ + (
			    j + k * d_dim2) * d_dim1] + w[i__ - 1 + (j + k * 
			    w_dim2) * w_dim1] * d__[i__ - 1 + (j + k * d_dim2)
			     * d_dim1] + w[i__ + 1 + (j + k * w_dim2) * 
			    w_dim1] * d__[i__ + 1 + (j + k * d_dim2) * d_dim1]
			     + w[i__ + (j - 1 + k * w_dim2) * w_dim1] * d__[
			    i__ + (j - 1 + k * d_dim2) * d_dim1] + w[i__ + (j 
			    + 1 + k * w_dim2) * w_dim1] * d__[i__ + (j + 1 + 
			    k * d_dim2) * d_dim1] + w[i__ + (j + (k - 1) * 
			    w_dim2) * w_dim1] * d__[i__ + (j + (k - 1) * 
			    d_dim2) * d_dim1] + w[i__ + (j + (k + 1) * w_dim2)
			     * w_dim1] * d__[i__ + (j + (k + 1) * d_dim2) * 
			    d_dim1]) / (d__[i__ + (j + k * d_dim2) * d_dim1] 
			    + d__[i__ - 1 + (j + k * d_dim2) * d_dim1] + d__[
			    i__ + 1 + (j + k * d_dim2) * d_dim1] + d__[i__ + (
			    j - 1 + k * d_dim2) * d_dim1] + d__[i__ + (j + 1 
			    + k * d_dim2) * d_dim1] + d__[i__ + (j + (k - 1) *
			     d_dim2) * d_dim1] + d__[i__ + (j + (k + 1) * 
			    d_dim2) * d_dim1]);
		}
		mp[ii] = 1e-20f;
/*                write(6,777) d(i,j,k)*d1, xp(ii), yp(ii), zp(ii) */
/*                write(6,7777) i, j, k */
/*                write(6,778) up(ii), vp(ii), wp(ii), mp(ii) */
/*                write(6,779) type(ii), tcp(ii), tdp(ii) */
/* L777: */
/* L7777: */
/* L778: */
/* L779: */

/*               write(7+procnum,1000) bmass*starfraction,tdp(ii),tcp(ii), */
/*     &                       metalf(ii) */
/*               write(7+procnum,1000) level,bmass*starfraction,tcp(ii), */
/*     &                           tdp(ii)*t1,d(i,j,k)*d1,z,metalf(ii) */
/* L1000: */

/*              Do not generate more star particles than available */

		if (ii >= *nmax) {
		    goto L20;
		}

L10:

		;
	    }
	}
    }
L20:

/*      if (ii .ge. nmax) then */
/*         write(6,*) 'pop3_maker: reached max new particle count' */
/*         stop */
/*      endif */
    *np = ii;
/*      if (np .ne. 0) then */
/*         write(6,*) 'pop3_maker: number,time,level: ', np, t, level, */
/*     $        irad */
/*      endif */

    return 0;
} /* pop3_color_maker__ */

#endif
