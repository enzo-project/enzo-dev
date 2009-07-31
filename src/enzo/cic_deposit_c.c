#ifdef CONFIG_PFLOAT_16
/* cic_deposit.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "macros_and_parameters.h"

/* ======================================================================= */
/* //////////////////////  SUBROUTINE CIC_DEPOSIT  \\\\\\\\\\\\\\\\\\\\\\c */
/* Subroutine */ 
int cic_deposit_c(FLOAT *posx, FLOAT *posy, FLOAT *posz, 
		int *ndim, int *npositions, float *mass, float *field, 
		FLOAT *leftedge, int *dim1, int *dim2, int *dim3, 
		FLOAT *cellsize)
{
    /* System generated locals */
    int field_dim1, field_dim2, field_offset, i__1;
    float r__1, r__2;

    /* Local variables */
    static int n, i1, j1, k1;
    static float dx, dy, dz, fact, xpos, ypos, zpos;
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
    fact = 1.f / *cellsize;
    edge1 = (float) (*dim1) - 0.5001;
    edge2 = (float) (*dim2) - 0.5001;
    edge3 = (float) (*dim3) - 0.5001;

/*     1D */

    if (*ndim == 1) {

	i__1 = *npositions;
	for (n = 1; n <= i__1; ++n) {

/*           Compute the position of the central cell */

/* Computing MIN */
/* Computing MAX */
	    r__2 = (posx[n] - leftedge[1]) * fact;
	    r__1 = MAX(r__2,0.5001);
	    xpos = MIN(r__1,edge1);

/*           Convert this into an int index */

	    i1 = (int) (xpos + .5f);

/*           Compute the weights */

	    dx = (float) i1 + .5f - xpos;

/*           Interpolate from field into sumfield */

	    field[i1 + (field_dim2 + 1) * field_dim1] += mass[n] * dx;
	    field[i1 + 1 + (field_dim2 + 1) * field_dim1] += mass[n] * (1.f - 
		    dx);

	}

    }

/*     2D */

    if (*ndim == 2) {

	i__1 = *npositions;
	for (n = 1; n <= i__1; ++n) {

/*           Compute the position of the central cell */

/* Computing MIN */
/* Computing MAX */
	    r__2 = (posx[n] - leftedge[1]) * fact;
	    r__1 = MAX(r__2,0.5001);
	    xpos = MIN(r__1,edge1);
/* Computing MIN */
/* Computing MAX */
	    r__2 = (posy[n] - leftedge[2]) * fact;
	    r__1 = MAX(r__2,0.5001);
	    ypos = MIN(r__1,edge2);

/*           Convert this into an int index */

	    i1 = (int) (xpos + .5f);
	    j1 = (int) (ypos + .5f);

/*           Compute the weights */

	    dx = (float) i1 + .5f - xpos;
	    dy = (float) j1 + .5f - ypos;

/*           Interpolate from field into sumfield */

	    field[i1 + (j1 + field_dim2) * field_dim1] += mass[n] * dx * dy;
	    field[i1 + 1 + (j1 + field_dim2) * field_dim1] += mass[n] * (1.f 
		    - dx) * dy;
	    field[i1 + (j1 + 1 + field_dim2) * field_dim1] += mass[n] * dx * (
		    1.f - dy);
	    field[i1 + 1 + (j1 + 1 + field_dim2) * field_dim1] += mass[n] * (
		    1.f - dx) * (1.f - dy);

	}

    }

/*     3D */

    if (*ndim == 3) {

	i__1 = *npositions;
	for (n = 1; n <= i__1; ++n) {

/*           Compute the position of the central cell */

/* Computing MIN */
/* Computing MAX */
	    r__2 = (posx[n] - leftedge[1]) * fact;
	    r__1 = MAX(r__2,0.5001);
	    xpos = MIN(r__1,edge1);
/* Computing MIN */
/* Computing MAX */
	    r__2 = (posy[n] - leftedge[2]) * fact;
	    r__1 = MAX(r__2,0.5001);
	    ypos = MIN(r__1,edge2);
/* Computing MIN */
/* Computing MAX */
	    r__2 = (posz[n] - leftedge[3]) * fact;
	    r__1 = MAX(r__2,0.5001);
	    zpos = MIN(r__1,edge3);

/*           Convert this into an integer index */

	    i1 = (int) (xpos + .5f);
	    j1 = (int) (ypos + .5f);
	    k1 = (int) (zpos + .5f);

/*           Compute the weights */

	    dx = (float) i1 + .5f - xpos;
	    dy = (float) j1 + .5f - ypos;
	    dz = (float) k1 + .5f - zpos;

/*           Interpolate from field into sumfield */

	    field[i1 + (j1 + k1 * field_dim2) * field_dim1] += mass[n] * dx * 
		    dy * dz;
	    field[i1 + 1 + (j1 + k1 * field_dim2) * field_dim1] += mass[n] * (
		    1.f - dx) * dy * dz;
	    field[i1 + (j1 + 1 + k1 * field_dim2) * field_dim1] += mass[n] * 
		    dx * (1.f - dy) * dz;
	    field[i1 + 1 + (j1 + 1 + k1 * field_dim2) * field_dim1] += mass[n]
		     * (1.f - dx) * (1.f - dy) * dz;
	    field[i1 + (j1 + (k1 + 1) * field_dim2) * field_dim1] += mass[n] *
		     dx * dy * (1.f - dz);
	    field[i1 + 1 + (j1 + (k1 + 1) * field_dim2) * field_dim1] += mass[
		    n] * (1.f - dx) * dy * (1.f - dz);
	    field[i1 + (j1 + 1 + (k1 + 1) * field_dim2) * field_dim1] += mass[
		    n] * dx * (1.f - dy) * (1.f - dz);
	    field[i1 + 1 + (j1 + 1 + (k1 + 1) * field_dim2) * field_dim1] += 
		    mass[n] * (1.f - dx) * (1.f - dy) * (1.f - dz);

	}

    }

    return 0;
} /* cic_deposit__ */

#endif
