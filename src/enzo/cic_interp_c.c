#ifdef CONFIG_PFLOAT_16
/* cic_interp.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

#define ENZO_PYTHON_IMPORTED
#include "macros_and_parameters.h"
#undef ENZO_PYTHON_IMPORTED

/* ======================================================================= */
/* //////////////////////  SUBROUTINE CIC_INTERP  \\\\\\\\\\\\\\\\\\\\\\\c */
/* Subroutine */ 
int cic_interp_c(FLOAT *posx, FLOAT *posy, FLOAT *posz, int *ndim,
	       int *npositions, float *sumfield, float *field,
	       FLOAT *leftedge, int *dim1, int *dim2, int *dim3, 
	       FLOAT *cellsize)
{
    /* System generated locals */
    int field_dim1, field_dim2, field_offset, i__1;
    FLOAT r__1, r__2;

    /* Local variables */
    static int n, i1, j1, k1;
    static FLOAT dx, dy, dz, fact,  xpos, ypos, zpos;
    static FLOAT edge1, edge2, edge3;


/*  PERFORMS 1/2/3D CLOUD-IN-CELL INTERPOLATION FROM FIELD TO SUMFIELD */

/*  written by: Greg Bryan */
/*  date:       January, 1998 */
/*  modified1: */

/*  PURPOSE: This routine performs a three-dimension, second-order */
/*           interpolation from field to sumfield (without clearing sumfield */
/*           first) at the positions specified. */

/*  INPUTS: */
/*     ndim       - dimensionality */
/*     cellsize   - the cell size of field */
/*     dim1,2,3   - real dimensions of field */
/*     field      - field to interpolate from */
/*     leftedge   - the left edge(s) of field */
/*     npositions - number of particles */
/*     posx,y,z   - particle positions */
/*     sumfield   - 1D field (length npositions) to interpolate into */

/*  OUTPUT ARGUMENTS: */
/*     sumfield   - 1D field (length npositions) to interpolate into */

/*  EXTERNALS: */

/*  LOCALS: */

/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */

/*  argument declarations */


/*  locals */


/*  constants (half is slightly more than 1/2 to make sure rounding is */
/*             performed in the correct direction) */


/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////// */
/* ======================================================================= */

    /* Parameter adjustments */
    --sumfield;
    --posz;
    --posy;
    --posx;
    --leftedge;
    field_dim1 = *dim1;
    field_dim2 = *dim2;
    field_offset = 1 + field_dim1 * (1 + field_dim2);
    field -= field_offset;

    /* Function Body */
    fact = (FLOAT)1. / *cellsize;
    edge1 = (FLOAT) (*dim1) - (FLOAT).500;
    edge2 = (FLOAT) (*dim2) - (FLOAT).500;
    edge3 = (FLOAT) (*dim3) - (FLOAT).500;

/*     1D */

    if (*ndim == 1) {

	i__1 = *npositions;
	for (n = 1; n <= i__1; ++n) {

/*           Compute the position of the central cell */

/* Computing min */
/* Computing max */
	    r__2 = (posx[n] - leftedge[1]) * fact;
	    r__1 = max(r__2,(FLOAT).500);
	    xpos = min(r__1,edge1);

/*           Convert this into an integer index */

	    i1 = (int) (xpos + (FLOAT).5);

/*           Compute the weights */

	    dx = (FLOAT) i1 + (FLOAT).5 - xpos;

/*           Interpolate from field into sumfield */

	    sumfield[n] = sumfield[n] + field[i1 + (field_dim2 + 1) * 
		    field_dim1] * dx + field[i1 + 1 + (field_dim2 + 1) * 
		    field_dim1] * ((FLOAT)1. - dx);

	}

    }

/*     2D */

    if (*ndim == 2) {

	i__1 = *npositions;
	for (n = 1; n <= i__1; ++n) {

/*           Compute the position of the central cell */

/* Computing min */
/* Computing max */
	    r__2 = (posx[n] - leftedge[1]) * fact;
	    r__1 = max(r__2,(FLOAT).500);
	    xpos = min(r__1,edge1);
/* Computing min */
/* Computing max */
	    r__2 = (posy[n] - leftedge[2]) * fact;
	    r__1 = max(r__2,(FLOAT).500);
	    ypos = min(r__1,edge2);

/*           Convert this into an integer index */

	    i1 = (int) (xpos + (FLOAT).5);
	    j1 = (int) (ypos + (FLOAT).5);

/*           Compute the weights */

	    dx = (FLOAT) i1 + (FLOAT).5 - xpos;
	    dy = (FLOAT) j1 + (FLOAT).5 - ypos;

/*           Interpolate from field into sumfield */

	    sumfield[n] = sumfield[n] + field[i1 + (j1 + field_dim2) * 
		    field_dim1] * dx * dy + field[i1 + 1 + (j1 + field_dim2) *
		     field_dim1] * ((FLOAT)1. - dx) * dy + field[i1 + (j1 + 1 
		    + field_dim2) * field_dim1] * dx * ((FLOAT)1. - dy) + 
		    field[i1 + 1 + (j1 + 1 + field_dim2) * field_dim1] * ((
		    FLOAT)1. - dx) * ((FLOAT)1. - dy);

	}

    }

/*     3D */

    if (*ndim == 3) {

	i__1 = *npositions;
	for (n = 1; n <= i__1; ++n) {

/*           Compute the position of the central cell */

/* Computing min */
/* Computing max */
	    r__2 = (posx[n] - leftedge[1]) * fact;
	    r__1 = max(r__2,(FLOAT).500);
	    xpos = min(r__1,edge1);
/* Computing min */
/* Computing max */
	    r__2 = (posy[n] - leftedge[2]) * fact;
	    r__1 = max(r__2,(FLOAT).500);
	    ypos = min(r__1,edge2);
/* Computing min */
/* Computing max */
	    r__2 = (posz[n] - leftedge[3]) * fact;
	    r__1 = max(r__2,(FLOAT).500);
	    zpos = min(r__1,edge3);

/*           Convert this into an integer index */

	    i1 = (int) (xpos + (FLOAT).5);
	    j1 = (int) (ypos + (FLOAT).5);
	    k1 = (int) (zpos + (FLOAT).5);

/*           Compute the weights */

	    dx = (FLOAT) i1 + (FLOAT).5 - xpos;
	    dy = (FLOAT) j1 + (FLOAT).5 - ypos;
	    dz = (FLOAT) k1 + (FLOAT).5 - zpos;

/*           Interpolate from field into sumfield */

	    sumfield[n] = sumfield[n] + 
	      field[i1+(j1 + k1 * field_dim2) * field_dim1] * dx * dy * dz + 
	      field[i1+1 + (j1 + k1 * field_dim2) * field_dim1] * ((FLOAT)1. - dx) * dy * dz + 
	      field[i1+(j1 + 1 + k1 * field_dim2) * field_dim1] * dx * ((FLOAT)1. - dy) * dz +
	      field[i1+1 + (j1 + 1 + k1 * field_dim2) * field_dim1] * ((FLOAT)1. - dx) *
	        ((FLOAT)1. - dy) * dz + 
	      field[i1 + (j1 + (k1 + 1) * field_dim2) * field_dim1] * dx * dy * ((FLOAT)1. - dz) +
	      field[i1 + 1 + (j1 + (k1 + 1) * field_dim2) * field_dim1] *
	        ((FLOAT)1. - dx) * dy * ((FLOAT)1. - dz) + 
	      field[i1 + (j1 + 1 + (k1 + 1) * field_dim2) * field_dim1] * dx *
	        ((FLOAT)1. - dy) * ((FLOAT)1. - dz) + 
	      field[i1 + 1 + (j1 + 1 + (k1 + 1) * field_dim2) * field_dim1] *
	        ((FLOAT)1. - dx) * ((FLOAT)1. - dy) * ((FLOAT)1. - dz);
	}
    }

    return 0;
} /* cic_interp__ */


#endif
