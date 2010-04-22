/***********************************************************************

  CREATES MBH PARTICLES

  written by: Ji-hoon Kim
  date:       April, 2010
  modified1:

  PURPOSE: Insert MBH particles at the location predetermined
           by external input file MBHInsertLocationFilename

  INPUTS:
    d     - density field
    u,v,w - velocity fields
    r     - refinement field (non-zero if zone is further refined)
    dt    - current timestep
    dx    - zone size (code units)
    t     - current time
    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
    nx,ny,nz - dimensions of field arrays
    ibuff    - number of buffer zones at each end of grid
    level    - current level of refinement
    nmax     - particle array size specified by calling routine
    x/y/z start - starting position of grid origin

  OUTPUTS:
    np   - number of particles created
    xp,yp,zp - positions of created particles
    up,vp,wp - velocities of created particles
    mp       - mass of new particles
    tcp   - particle creation time (0 for dm particles)
    tdp   - particle dynamic time (not used here)

***********************************************************************/
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"

int mbh_maker(int *nx, int *ny, int *nz, int *size, float *d, float *u, 
	      float *v, float *w, float *dt, float *r, float *dx, FLOAT *t, 
	      float *d1, float *x1, float *v1, float *t1, 
	      int *nmax, FLOAT *xstart, FLOAT *ystart, 
	      FLOAT *zstart, int *ibuff, 
	      int *level, int *np, FLOAT *xp, FLOAT *yp, 
	      FLOAT *zp, float *up, float *vp, float *wp, float *mp, 
	      float *tcp, float *tdp, int *type, int *ctype)
{

  int i, j, k, ii, index;
  FILE *fptr;
  double dummy_double[3], msolar = 1.989e33;
  float dummy_float[3];
  char line[MAX_LINE_LENGTH];

  /* Set new particle index to number of created star particles */

  ii = *np;


  /* Open the file and read in the MBH masses and locations */

  if ((fptr = fopen(MBHInsertLocationFilename, "r")) == NULL) {

    fprintf(stderr, "mbh_maker: Error opening file %s\n", MBHInsertLocationFilename);
    return FAIL;

  } else {

    /* Now read the files line by line */

    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) 
      if (line[0] != '#') {

	/* order: MBH mass (in Ms), MBH position[3], MBH creation time */

	if (sscanf(line, " %"FSYM"  %"PSYM"  %"PSYM"  %"PSYM"  %"FSYM, 
		   &dummy_float[0], &dummy_double[0], &dummy_double[1],
		   &dummy_double[2], &dummy_float[1]) != 5) {
	  fprintf(stderr, "mbh_maker: File structure wrong in %s\n", MBHInsertLocationFilename);
	  return FAIL;
	}

	/* Find the indices in the grid */

	i = (float) ((dummy_double[0] - *xstart) / (*dx));
	j = (float) ((dummy_double[1] - *ystart) / (*dx));
	k = (float) ((dummy_double[2] - *zstart) / (*dx));

	/* Only if the indices are all good-looking, and if this grid is at the 
	   finest level in the hierarchy, then insert the particles; otherwise,
	   this particle should be created in one of the other grids! */

	if ( (i >= *ibuff) && (i < (*nx - *ibuff)) &&
	     (j >= *ibuff) && (j < (*ny - *ibuff)) &&
	     (k >= *ibuff) && (k < (*nz - *ibuff)) ) {

	  index = (k * (*ny) + j) * (*nx) + i;

//	  fprintf(stdout, "mbh_maker: i, j, k, index, r[index] = %d, %d, %d, %d, %g\n", 
//		  i, j, k, index, r[index]);  

	  if (r[index] == 0.0) {
	  
	    // MBH mass and position
	    mp[ii] = (float) ( dummy_float[0] * msolar / 
			       ((double)(*d1) * pow((double)((*x1)*(*dx)), 3)) );
	    xp[ii] = dummy_double[0]; 
	    yp[ii] = dummy_double[1];
	    zp[ii] = dummy_double[2];
	    
	    // if creation time is not properly given, take the current time
	    if (dummy_float[1] <= 0.0)
	      tcp[ii] = (float) *t;
	    else
	      tcp[ii] = dummy_float[1];
	    
	    // MBH other attributes 
	    up[ii] = u[index];
	    vp[ii] = v[index];
	    wp[ii] = w[index];
	    tdp[ii] = MBHMinDynamicalTime;
	    
	    // when the creation time is still in the future, type is negative
	    if (tcp[ii] <= (*t))
	      type[ii] = *ctype;
	    else 
	      type[ii] = -(*ctype);

	    // take out the mass from the cell; if you then get a negative density,
	    // probably your choice of the MBH position was horribly wrong - NOW UNUSED!

	    d[index] -= mp[ii];
	    if (d[index] <= 0) {
	      fprintf(stderr, "mbh_maker: check the MBH position; not dense enough for your MBH mass. (mp[]=%g, d[]=%g)\n",
		      mp[ii], d[index]);
	      return FAIL;
	    }

	    fprintf(stdout, "mbh_maker: a mbh inserted at (%g,%g,%g) with v=(%g,%g,%g), m=%g, tc=%g, type=%d\n",
		    xp[ii], yp[ii], zp[ii], up[ii], vp[ii], wp[ii], mp[ii], tcp[ii], type[ii]);

	    // increase the counter
	    ii++;
	    
	  } //ENDIF r[index] == 0.0

	} //ENDIF indices correct

      } //ENDIF line[0] != '#'

  } //ENDIF fptr


  if ((ii-*np) > 0)
    fprintf(stdout, "mbh_maker: %"ISYM" new mbh particle(s)\n", ii-*np);

  if ((ii-*np) >= *nmax) {
    fprintf(stdout, "mbh_maker: reached max new particle count");
    return FAIL;
  }

  *np = ii;
  return SUCCESS;

}
