#ifdef CONFIG_PFLOAT_16
/* cic_flag.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

#define ENZO_PYTHON_IMPORTED
#include "macros_and_parameters.h"
#undef ENZO_PYTHON_IMPORTED

/* ======================================================================= */
/* //////////////////////  SUBROUTINE CIC_FLAG  \\\\\\\\\\\\\\\\\\\\\\c */
/* Subroutine */ 
int cic_flag_c(FLOAT *posx, FLOAT *posy, FLOAT *posz, int *ndim, int *npositions, 
	     int *itype, int *ffield, FLOAT *leftedge, int *dim1, int *dim2, 
	       int *dim3, float *cellsize, int *imatch1, int *imatch2,
	       int *buffersize)
{
    /* System generated locals */
    int ffield_dim1, ffield_dim2, ffield_offset, i__1;
    FLOAT r__1, r__2;

    /* Local variables */
    static int n, i1, j1, k1;
    static FLOAT fact,  xpos, ypos, zpos;
    static FLOAT edge1, edge2, edge3;

/*  PERFORMS 1/2/3D CLOUD-IN-CELL MARKING OF FLAGGING FIELD */

/*  written by: Greg Bryan */
/*  date:       September, 2003 */
/*  modified1: */

/*  PURPOSE: This routine performs a 1,2 or 3 dimension setting of a */
/*           flagging field for must-refine particles. */

/*  INPUTS: */
/*     ndim       - dimensionality */
/*     cellsize   - the cell size of field */
/*     dim1,2,3   - real dimensions of field */
/*     leftedge   - the left edge(s) of field */
/*     npositions - number of particles */
/*     posx,y,z   - particle positions */
/*     itype      - 1D field (length npositions) of types */
/*     imatch     - integer indicating type of particle to match */

/*  OUTPUT ARGUMENTS: */
/*     ffield      - field to be deposited to */

/*  EXTERNALS: */

/*  LOCALS: */

/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */

/*  argument declarations */


/*  locals */


/* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////// */
/* ======================================================================= */

    /* Parameter adjustments */
    --itype;
    --posz;
    --posy;
    --posx;
    --leftedge;
    ffield_dim1 = *dim1;
    ffield_dim2 = *dim2;
    ffield_offset = 1 + ffield_dim1 * (1 + ffield_dim2);
    ffield -= ffield_offset;

    /* Function Body */
    fact = 1.f / *cellsize;
    edge1 = (FLOAT) (*dim1) - 0.5001;
    edge2 = (FLOAT) (*dim2) - 0.5001;
    edge3 = (FLOAT) (*dim3) - 0.5001;

/*     1D */

    if (*ndim == 1) {

	i__1 = *npositions;
	for (n = 1; n <= i__1; ++n) {

/*           only do this for must-refine particles */

	    if (itype[n] == *imatch1 || itype[n] == *imatch2) {

/*           Compute the position of the central cell */

/* Computing min */
/* Computing max */
		r__2 = (posx[n] - leftedge[1]) * fact;
		r__1 = max(r__2,0.5001);
		xpos = min(r__1,edge1);

/*           Convert this into an integer index */

		i1 = (int) (xpos + .5f);

/*           set flagging field */

		++ffield[i1 + (ffield_dim2 + 1) * ffield_dim1];
		++ffield[i1 + 1 + (ffield_dim2 + 1) * ffield_dim1];

	    }

	}

    }

/*     2D */

    if (*ndim == 2) {

	i__1 = *npositions;
	for (n = 1; n <= i__1; ++n) {

/*           only do this for must-refine particles */

	    if (itype[n] == *imatch1 || itype[n] == *imatch2) {

/*           Compute the position of the central cell */

/* Computing min */
/* Computing max */
		r__2 = (posx[n] - leftedge[1]) * fact;
		r__1 = max(r__2,0.5001);
		xpos = min(r__1,edge1);
/* Computing min */
/* Computing max */
		r__2 = (posy[n] - leftedge[2]) * fact;
		r__1 = max(r__2,0.5001);
		ypos = min(r__1,edge2);

/*           Convert this into an integer index */

		i1 = (int) (xpos + .5f);
		j1 = (int) (ypos + .5f);

/*           Interpolate from field into sumfield */

		++ffield[i1 + (j1 + ffield_dim2) * ffield_dim1];
		++ffield[i1 + 1 + (j1 + ffield_dim2) * ffield_dim1];
		++ffield[i1 + (j1 + 1 + ffield_dim2) * ffield_dim1];
		++ffield[i1 + 1 + (j1 + 1 + ffield_dim2) * ffield_dim1];

	    }

	}

    }

/*     3D */

    if (*ndim == 3) {

	i__1 = *npositions;
	for (n = 1; n <= i__1; ++n) {

/*           only do this for must-refine particles */

	    if (itype[n] == *imatch1 || itype[n] == *imatch2) {

/*           Compute the position of the central cell */

/* Computing min */
/* Computing max */
		r__2 = (posx[n] - leftedge[1]) * fact;
		r__1 = max(r__2,0.5001);
		xpos = min(r__1,edge1);
/* Computing min */
/* Computing max */
		r__2 = (posy[n] - leftedge[2]) * fact;
		r__1 = max(r__2,0.5001);
		ypos = min(r__1,edge2);
/* Computing min */
/* Computing max */
		r__2 = (posz[n] - leftedge[3]) * fact;
		r__1 = max(r__2,0.5001);
		zpos = min(r__1,edge3);

/*           Convert this into an integer index */

		i1 = (int) (xpos + .5f);
		j1 = (int) (ypos + .5f);
		k1 = (int) (zpos + .5f);

/*           Set flagging field */

		++ffield[i1 + (j1 + k1 * ffield_dim2) * ffield_dim1];
		++ffield[i1 + 1 + (j1 + k1 * ffield_dim2) * ffield_dim1];
		++ffield[i1 + (j1 + 1 + k1 * ffield_dim2) * ffield_dim1];
		++ffield[i1 + 1 + (j1 + 1 + k1 * ffield_dim2) * ffield_dim1];
		++ffield[i1 + (j1 + (k1 + 1) * ffield_dim2) * ffield_dim1];
		++ffield[i1 + 1 + (j1 + (k1 + 1) * ffield_dim2) * ffield_dim1]
			;
		++ffield[i1 + (j1 + 1 + (k1 + 1) * ffield_dim2) * ffield_dim1]
			;
		++ffield[i1 + 1 + (j1 + 1 + (k1 + 1) * ffield_dim2) * 
			ffield_dim1];

	    }

	}

    }

    return 0;
} /* cic_flag__ */

#endif
