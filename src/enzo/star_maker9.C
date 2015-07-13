/***********************************************************************
  CREATES AND ADDS TO SINK PARTICLES - based on star_maker8

  written by: Elizabeth Harper-Clark
  date:       Sept 2010

  INPUTS:
    d     - density field
    u,v,w - velocity fields
    r     - refinement field (non-zero if zone is further refined)
    dt    - current timestep
    dx    - zone size (code units)
    t     - current time
    z     - current redshift
    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
    nx,ny,nz - dimensions of field arrays
    ibuff    - number of buffer zones at each end of grid
    imethod  - Hydro method (0/1 -- PPM DE/LR, 2 - ZEUS)
    massthresh - maximum allowed mass in a cell
    level - current level of refinement
    procnum - processor number (for output)
    npold - number of sink particles already on grid
    xp-zpold - current sink particle positions
    up-wpold - current sink particle velocities
    mpold - curernt sink particle masses
    tcp   - star particle creation time (0 for dm particles)
    tdp   - star particle dynamic time (not used here)
    nmax     - particle array size specified by calling routine
    x/y/z start - starting position of grid origin
    jlrefine - Jeans length refinement (<0 for none)
    temp - temperature

  OUTPUTS:
    np   - number of particles created
    xp,yp,zp - positions of created particles
    up,vp,wp - velocities of created particles
    mp       - mass of new particles
***********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"
// #include "CommunicationUtilities.h"

#define USE

/* function prototypes */
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

int  GetUnits(float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MassUnits, FLOAT Time);

int star_maker9(int *nx, int *ny, int *nz, int *size, float *d, float *te, float *ge, 
		float *u, float *v, float *w, float *bx, float *by, float *bz,
		float *dt, float *r, float *dx, FLOAT *t, 
		float *z, int *procnum, float *d1, float *x1, float *v1, 
		float *t1, int *nmax, FLOAT *xstart, FLOAT *ystart, 
		FLOAT *zstart, int *ibuff, int *imethod, int *idual,
		float *massthresh, int *level, int *np, FLOAT *xp, FLOAT *yp, 
		FLOAT *zp, float *up, float *vp, float *wp, float *mp, 
		float *tcp, float *tdp, float *dm, 
		int *type, int *npold, FLOAT *xpold, 
		FLOAT *ypold, FLOAT *zpold, float *upold, float *vpold, 
		float *wpold, float *mpold, float *tcpold, float *tdpold, float *dmold,
		float *nx_jet, float *ny_jet, float *nz_jet,
		int *typeold, PINT *idold, int *ctype, float *jlrefine, 
		float *temp, float *gamma, float *mu, int *nproc, int *nstar)
{

  int		i, j, k, index, ii, inew, n, bb, cc, nsinks, closest;
  int           xo, yo, zo;
#define MAX_SUPERCELL_NUMBER 1000
  int           n_cell, ind_cell[MAX_SUPERCELL_NUMBER];
  float		densthresh, maxdens, adddens, ugrid, vgrid, wgrid, m_cell;
  double	jeansthresh, jlsquared, dx2, dist2, total_density, nearestdx2;
  FLOAT		xpos, ypos, zpos, delx, dely, delz;
  double        DensityFloor;
  double        Pi = 3.1415926;
  float nx_cell[MAX_SUPERCELL_NUMBER], 
    ny_cell[MAX_SUPERCELL_NUMBER], nz_cell[MAX_SUPERCELL_NUMBER];
  double msun = 1.989e33;
  double umass = (*d1)*POW(*x1,3)/msun;

  //printf("Star Maker 9 running.......\n");

  /* Convert mass threshold to density */

  densthresh = *massthresh / POW(*dx,3);
  //densthresh = 1e-12/(*d1);
  dx2 = (*dx)*(*dx);
  //printf("Star Maker 8: densthresh = %g\n", densthresh);
  if (*jlrefine > 0) {
    jlsquared = ((double)((*gamma) * 3.14159 * 1.38e-16 / 6.673e-08) / 
		 ((double)(*d1) * 1.673e-24)) / POW(*x1,2) / (*mu) / POW((*jlrefine),2);
  }

  /* Set new particle index to number of created star particles */

  ii = *np;

  /* Look for existing sink particles */

  nsinks = 0;
  int *sink_index = new int[50000];
  for (n = 0; n < *npold; n++) {
    if (typeold[n] == *ctype) {
      sink_index[nsinks++] = n;
    }
  }

  /* For a 3D->1D index, put x+1, y+1, z+1 indices into nice variables */

  xo = 1;
  yo = *nx;
  zo = (*nx) * (*ny);

  //BigStarFormationDone = CommunicationMaxValue(BigStarFormationDone);

  /* Loop over grid looking for a cell with mass larger than massthres */
  //printf("BigStarFormationDone = %"ISYM" MyProcessorNumber = %"ISYM"\n", BigStarFormationDone,MyProcessorNumber);
  if(BigStarFormationDone == 0){
    if (*level == MaximumRefinementLevel) {
      float oldrho;
      float SinkCollapseDistance = SinkMergeDistance;
      for (k = *ibuff; k < *nz-*ibuff; k++) {
	for (j = *ibuff; j < *ny-*ibuff; j++) {
	  index = (k * (*ny) + j) * (*nx) + (*ibuff);
	  for (i = *ibuff; i < *nx-*ibuff; i++, index++) {

	    /* Finest level of refinement and density greater than threshold? */
	  
	    if (*jlrefine > 0)
	      jeansthresh = jlsquared * temp[index] / d[index];
	    //printf("jeansthresh = %g \n",jeansthresh);
	    //printf("jlsquared = %g \n",jlsquared);printf("temp[index] = %g \n",temp[index]);printf("d[index] = %g \n",d[index]);

	    if (r[index] == 0 && (d[index] > densthresh ||
				  (*jlrefine > 0 && dx2 > jeansthresh))) {
	      //printf("star_maker9: density above thresh-hold - will now make a new star?!\n");
	      xpos = *xstart + ((float) i - 0.5)*(*dx);
	      ypos = *ystart + ((float) j - 0.5)*(*dx);
	      zpos = *zstart + ((float) k - 0.5)*(*dx);

	      nearestdx2 = 1e20;
	      //float BigStarSeparation = (*x1)/4;
	      // printf("BigStarSeparation = %g = %g cgs \n", BigStarSeparation, BigStarSeparation*(*x1) );

	      for (cc = 0; cc < nsinks; cc++) {
	      
		n = sink_index[cc];
		if (mpold[n]*umass < 3.0) continue;
		delx = xpos - xpold[n];
		dely = ypos - ypold[n];
		delz = zpos - zpold[n];
		dist2 = delx*delx + dely*dely + delz*delz;

		/* If sink is within 5 cells of the closest one, then add to it */
		if (dist2 < POW(BigStarSeparation,2) && dist2 < nearestdx2) {
		  nearestdx2 = dist2;
		  closest = n;		  
		}

	      } // ENDFOR old particles
	      //printf("star_maker9: nearest old star = %"FSYM"\n",POW(nearestdx2,0.5) );

	      if (ii < *nmax) {

		// PUT BIG STAR FORMATION IF STATEMENT HERE
		if(BigStarFormation == 1){

		  /* Calculate change in density */

// 		  if (*jlrefine > 0)
// 		    maxdens = min(jlsquared * temp[index] / dx2, densthresh);
// 		  else
// 		    maxdens = densthresh;
// 		  oldrho = d[index];
// 		  adddens = d[index] - maxdens;
		  BigStarFormationDone = 1;
		  //StarParticleCreation = 0;
		  //StarParticleFeedback = 0;
		  CommunicationBroadcastValue(&BigStarFormationDone, MyProcessorNumber);
		  //CommunicationBroadcastValue(&StarParticleCreation, MyProcessorNumber);
		  //CommunicationBroadcastValue(&StarParticleFeedback, MyProcessorNumber);
	    
		  /* Remove mass from grid */
	    
// 		  d[index] = maxdens;
		  printf("BigStarFormation: Star made at %"FSYM", %"FSYM", %"FSYM" \n ", xpos, ypos, zpos);
		  //printf("now BigStarFormation = %"ISYM"\n", BigStarFormation);

		  if (*imethod == 2) {
		    ugrid = 0.5*(u[index] + u[index+xo]);
		    vgrid = 0.5*(v[index] + v[index+yo]);
		    wgrid = 0.5*(w[index] + w[index+zo]);
		  } else {
		    total_density = d[index] + d[index-xo] + d[index+xo] + d[index-yo] +
		      d[index+yo] + d[index-zo] + d[index+zo];

		    ugrid = (u[index]*d[index] + 
			     u[index-xo]*d[index-xo] + u[index+xo]*d[index+xo] +
			     u[index-yo]*d[index-yo] + u[index+yo]*d[index+yo] +
			     u[index-zo]*d[index-zo] + u[index+zo]*d[index+zo]) /
		      total_density;

		    vgrid = (v[index]*d[index] + 
			     v[index-xo]*d[index-xo] + v[index+xo]*d[index+xo] +
			     v[index-yo]*d[index-yo] + v[index+yo]*d[index+yo] +
			     v[index-zo]*d[index-zo] + v[index+zo]*d[index+zo]) /
		      total_density;

		    wgrid = (w[index]*d[index] + 
			     w[index-xo]*d[index-xo] + w[index+xo]*d[index+xo] +
			     w[index-yo]*d[index-yo] + w[index+yo]*d[index+yo] +
			     w[index-zo]*d[index-zo] + w[index+zo]*d[index+zo]) /
		      total_density;
		  }


		  printf("star_maker9: making new star, type = %"ISYM"\n",*ctype );
		  mp[ii] = 0.0; //adddens;
		  type[ii] = -*ctype;
	      
		  /* Set positions and velocities */
	    
		  xp[ii] = xpos;
		  yp[ii] = ypos;
		  zp[ii] = zpos;;
		  up[ii] = ugrid;
		  vp[ii] = vgrid;
		  wp[ii] = wgrid;

		  /* Set creation time */
	      
		  tcp[ii] = (float) *t;
		  tdp[ii] = 1.0e20;
		  dm[ii]  = 0.0; //adddens*POW(*dx,3);

		  ii++;



		}
		else {
		  printf("star_maker 9 called but not BigStarFormation \n");
		}

	      } // ENDIF create a new sink
	    
	    } // ENDIF make sink particle

	  } // ENDFOR i
	} // ENDFOR j
      } // ENDFOR k

    } // if (level == maxlevel)
  } // if BigStarForm

  if (ii > 0)
    printf("P(%"ISYM"): star_maker9[add]: %"ISYM" new sink particles\n", *nproc, ii);

  if (ii >= *nmax) {
    fprintf(stdout, "star_maker9: reached max new particle count");
    return FAIL;
  }

  delete sink_index;

  *np = ii;
  return SUCCESS;

}
