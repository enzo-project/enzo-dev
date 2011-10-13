#ifdef CONFIG_PFLOAT_16
/* smooth_deposit.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include <stdio.h>
#include <math.h>

#define ENZO_PYTHON_IMPORTED
#include "macros_and_parameters.h"
#undef ENZO_PYTHON_IMPORTED

/* function proto types */
void my_exit(int status);


/* ======================================================================= */
/* //////////////////////  SUBROUTINE SMOOTH_DEPOSIT  \\\\\\\\\\\\\\\\\\\c */
/* Subroutine */ 
int smooth_deposit_c(FLOAT *posx, FLOAT *posy, FLOAT *posz, int *ndim, 
		   int *npositions, float *mass, float *field, FLOAT *leftedge, 
		   int *dim1, int *dim2, int *dim3, float *cellsize, 
		   float *rsmooth)
{

  /* Table of constant values */
  static int c__9 = 9;
  static int c__1 = 1;


    /* System generated locals */
    int field_dim1, field_dim2, field_offset, i__1, i__2, i__3, i__4;
    FLOAT r__1, r__2, r__3;

    /* Local variables */
    static int i__, j, k, n;
    static FLOAT coef;
    static FLOAT xpos, ypos, zpos, rsqr, rsmsqr;

/*  PERFORMS 3D SMOOTHED INTERPOLATION FROM FIELD TO SUMFIELD */

/*  written by: Greg Bryan */
/*  date:       January, 2000 */
/*  modified1: */

/*  PURPOSE: This routine performs a three-dimension, second-order */
/*           interpolation from field to sumfield (without clearing sumfield */
/*           first) at the positions specified. */

/*  INPUTS: */
/*     ndim       - dimensionality */
/*     cellsize   - the cell size of field */
/*     dim1,2,3   - real dimensions of field */
/*     leftedge   - the left edge(s) of field */
/*     npositions - number of particles */
/*     posx,y,z   - particle positions */
/*     sumfield   - 1D field (length npositions) of masses */

/*  OUTPUT ARGUMENTS: */
/*     field      - field to be deposited to */

/*  EXTERNALS: */

/*  LOCALS: */

/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */

/*  argument declarations */


/*  locals */


/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////// */
/* ======================================================================= */

    /* Parameter adjustments */
    --mass;
    --posz;
    --posy;
    --posx;
    --leftedge;
    field_dim1 = *dim1;
    field_dim2 = *dim2;
    field_offset = 1 + field_dim1 * (1 + field_dim2);
    field -= field_offset;

    /* Function Body */
/* Computing 2nd power */
    r__1 = *rsmooth;
    rsmsqr = r__1 * r__1;
/* Computing 3rd power */
    r__1 = *cellsize / *rsmooth;
    coef = r__1 * (r__1 * r__1) * ((FLOAT).9549304651466296);

    if (*ndim != 3) {
      fprintf(stderr, "SMOOTH_DEPOSIT: only ndim=3 supported.");
      exit(c__1);
    }

/*     3D */
    if (*ndim == 3) {

/*        loop over grid */

	i__1 = *dim3;
	for (k = 1; k <= i__1; ++k) {
	    zpos = leftedge[3] + ((FLOAT) k - .5f) * *cellsize;
	    i__2 = *dim2;
	    for (j = 1; j <= i__2; ++j) {
		ypos = leftedge[2] + ((FLOAT) j - .5f) * *cellsize;
		i__3 = *dim1;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    xpos = leftedge[1] + ((FLOAT) i__ - .5f) * *cellsize;

/*                 Loop over particles */

		    i__4 = *npositions;
		    for (n = 1; n <= i__4; ++n) {

/*                     Compute distance from particle to cell center */

/* Computing 2nd power */
			r__1 = posx[n] - xpos;
/* Computing 2nd power */
			r__2 = posy[n] - ypos;
/* Computing 2nd power */
			r__3 = posz[n] - zpos;
			rsqr = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;

			if (rsqr < rsmsqr) {
			    field[i__ + (j + k * field_dim2) * field_dim1] += 
				    mass[n] * coef * (((FLOAT)1.0) - sqrtl(rsqr) / *
				    rsmooth);
			}

		    }

		}
	    }
	}

    }

    return 0;
} /* smooth_deposit__ */

#endif
