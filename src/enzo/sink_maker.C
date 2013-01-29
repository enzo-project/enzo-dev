/***********************************************************************
  CREATES AND ADDS TO SINK PARTICLES

  written by: Greg Bryan
  date:       November, 2002
  modified1:

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
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"

int sink_maker(int *nx, int *ny, int *nz, int *size, float *d, float *u, 
		float *v, float *w, float *dt, float *r, float *dx, FLOAT *t, 
		float *z, int *procnum, float *d1, float *x1, float *v1, 
		float *t1, int *nmax, FLOAT *xstart, FLOAT *ystart, 
		FLOAT *zstart, int *ibuff, int *imethod, 
		float *massthresh, int *level, int *np, FLOAT *xp, FLOAT *yp, 
		FLOAT *zp, float *up, float *vp, float *wp, float *mp, 
		float *tcp, float *tdp, int *type, int *npold, FLOAT *xpold, 
		FLOAT *ypold, FLOAT *zpold, float *upold, float *vpold, 
		float *wpold, float *mpold, float *tcpold, float *tdpold, 
		int *typeold, int *ctype, float *jlrefine, float *temp, float *JRCT)
{

  int		i, j, k, index, ii, inew, n, dim, bb, cc, nsinks, closest;
  int           xo, yo, zo;
  float         huge_number = 1e20;
  float		densthresh, maxdens, adddens, ugrid, vgrid, wgrid, my_temp;
  double	jeansthresh, jlsquared, dx2, dist2, total_density, nearestdx2;
  FLOAT		xpos, ypos, zpos, delx, dely, delz;

  // When refining by Jeans, make a sink particle when the density
  // exceeds the critical value times by this factor
  float         overflowFactor = 1.01;

  /* Convert mass threshold to density */

  densthresh = *massthresh / pow(*dx,3);

  dx2 = (*dx)*(*dx);
  if (*jlrefine > 0)
    jlsquared = ((double)(3.14159 * 1.38e-16 / 6.673e-08) / 
		 ((double)(*d1) * 1.673e-24)) / pow(*x1,2) / pow(*jlrefine,2);

  /* Set new particle index to number of created star particles */

  ii = *np;

  /* Look for existing sink particles */

  nsinks = 0;
  int *sink_index = new int[10000];
  for (n = 0; n < *npold; n++)
    if (typeold[n] == *ctype && tcpold[n] > 0 && tdpold[n] == 0)
      sink_index[nsinks++] = n;

  /* Merge any sink particles that are within 5 cells of each other */

  double mfrac_b, mfrac_c, total_mass;
  float index_b[3], index_c[3];

  for (i = 0; i < nsinks-1; i++) {
    
    bb = sink_index[i];
    nearestdx2 = huge_number;

    if (mpold[bb] < 0) continue;

    for (j = i+1; j < nsinks; j++) {

      cc = sink_index[j];
      if (mpold[cc] < 0) continue;

      delx = xpold[bb] - xpold[cc];
      dely = ypold[bb] - ypold[cc];
      delz = zpold[bb] - zpold[cc];
      dist2 = delx*delx + dely*dely + delz*delz;

      if (dist2 < 25.0*dx2 && dist2 < nearestdx2) {
//	printf("star_maker3[merge]: (part %"ISYM" + %"ISYM"): dx = %"GSYM", nearestdx = %"GSYM"\n",
//	       bb, cc, sqrt(dist2/dx2), sqrt(nearestdx2/dx2));
	nearestdx2 = dist2;
	closest = cc;
      }

    } // ENDFOR second old particle
    
    /* If there are particles to be merged, do it now (check if
       nearestdx2 has changed).  New particle has the center of mass
       of the two close particles */

    if (nearestdx2 < huge_number) {

      total_mass = mpold[bb] + mpold[closest];
      mfrac_b = mpold[bb] / total_mass;
      mfrac_c = mpold[closest] / total_mass;

      index_b[0] = (xpold[bb] - *xstart) / (*dx);
      index_b[1] = (ypold[bb] - *ystart) / (*dx);
      index_b[2] = (zpold[bb] - *zstart) / (*dx);
      index_c[0] = (xpold[closest] - *xstart) / (*dx);
      index_c[1] = (ypold[closest] - *ystart) / (*dx);
      index_c[2] = (zpold[closest] - *zstart) / (*dx);

      printf("star_maker3[merge1]: %"ISYM" %"FSYM" %"FSYM" %"FSYM" %"GSYM" %"FSYM"\n", bb, 
	     index_b[0], index_b[1], index_b[2], mpold[bb], mfrac_b);
      printf("star_maker3[merge2]: %"ISYM" %"FSYM" %"FSYM" %"FSYM" %"GSYM" %"FSYM"\n", closest, 
	     index_c[0], index_c[1], index_c[2], mpold[closest], mfrac_c);

      xpold[bb] = *xstart + (*dx) * 
	(FLOAT) (index_b[0]*mfrac_b + index_c[0]*mfrac_c);
      ypold[bb] = *ystart + (*dx) * 
	(FLOAT) (index_b[1]*mfrac_b + index_c[1]*mfrac_c);
      zpold[bb] = *zstart + (*dx) * 
	(FLOAT) (index_b[2]*mfrac_b + index_c[2]*mfrac_c);
      
      upold[bb] = upold[bb]*(FLOAT)mfrac_b + upold[closest]*(FLOAT)mfrac_c;
      vpold[bb] = vpold[bb]*(FLOAT)mfrac_b + vpold[closest]*(FLOAT)mfrac_c;
      wpold[bb] = wpold[bb]*(FLOAT)mfrac_b + wpold[closest]*(FLOAT)mfrac_c;
      mpold[bb] = total_mass;
    
      printf("star_maker3[merge3]: %"FSYM" %"FSYM" %"FSYM" %"GSYM"\n", 
	     (float) ((xpold[bb] - *xstart) / (*dx)), 
	     (float) ((ypold[bb] - *ystart) / (*dx)),
	     (float) ((zpold[bb] - *zstart) / (*dx)), 
	     mpold[bb]);

      // Set second particle to be ignored (no mass)
      tcpold[closest]                                  = 0.0;
      upold[closest] = vpold[closest] = wpold[closest] = 0.0;
      mpold[closest] = FLOAT_UNDEFINED;

    }  // ENDIF merge particle (nearestdx2 < 1e20)

  } // ENDFOR first old particle

  /* Remove deleted particle from sink particle list */

  int nRemoved = 0;
  for (n = 0; n < nsinks; n++) {
    if (mpold[sink_index[n]] <= 0) {
      for (bb = n+1; bb < nsinks; bb++) 
	sink_index[bb-1] = sink_index[bb];
      nRemoved++;
    }
  }
  nsinks -= nRemoved;
//  printf("star_maker3[remove]: Ignoring %"ISYM" sink particles.\n", nRemoved);

  /* For a 3D->1D index, put x+1, y+1, z+1 indices into nice variables */

  xo = 1;
  yo = *nx;
  zo = (*nx) * (*ny);

  /* Loop over grid looking for a cell with mass larger than massthres */

  for (k = *ibuff; k < *nz-*ibuff; k++) {
    for (j = *ibuff; j < *ny-*ibuff; j++) {
      index = (k * (*ny) + j) * (*nx) + (*ibuff);
      for (i = *ibuff; i < *nx-*ibuff; i++, index++) {

	/* Finest level of refinement and density greater than threshold? */

    my_temp = (*JRCT > 0) ? *JRCT : temp[index];

	if (*jlrefine > 0) 
	  jeansthresh = overflowFactor * jlsquared * my_temp / d[index];

//	printf("star_maker3[a]: %"ISYM" %"ISYM" %"ISYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", i, j, k,
//	       densthresh, d[index], jeansthresh, dx2, *jlrefine, jlsquared);

	if (r[index] == 0 && (d[index] > densthresh ||
			      (*jlrefine > 0 && dx2 > jeansthresh))) {

	  xpos = *xstart + ((float) i - 0.5)*(*dx);
	  ypos = *ystart + ((float) j - 0.5)*(*dx);
	  zpos = *zstart + ((float) k - 0.5)*(*dx);
	  
	  /* Calculate change in density */

	  if (*jlrefine > 0)
	    maxdens = min(jlsquared * my_temp / dx2, densthresh);
	  else
	    maxdens = densthresh;
	  
	  adddens = d[index] - maxdens;

	  /* Remove mass from grid */

	  d[index] = maxdens;

	  /* Get velocity at grid center (averaged over cells ... ONLY
	     for PPM) */

	  if (*imethod == 2) {
	    ugrid = 0.5*(u[index] + u[index+xo]);
	    vgrid = 0.5*(v[index] + v[index+yo]);
	    wgrid = 0.5*(w[index] + w[index+zo]);
	  } else {
	    ugrid = u[index];
	    vgrid = v[index];
	    wgrid = w[index];
	  }

	  /* Look for a nearby old sink particle to add the mass to */

	  inew = 1;
	  nearestdx2 = 1e20;
	  for (cc = 0; cc < nsinks; cc++) {
	    
	    n = sink_index[cc];

	    delx = xpos - xpold[n];
	    dely = ypos - ypold[n];
	    delz = zpos - zpold[n];
	    dist2 = delx*delx + dely*dely + delz*delz;

	    /* If sink is within 5 cells and closest one, then add to it */

	    if (dist2 < 25.0*dx2 && dist2 < nearestdx2) {
//	      printf("star_maker3[addold]: (part %"ISYM"): dx = %"GSYM", nearestdx = %"GSYM"\n",
//		     cc, sqrt(dist2/dx2), sqrt(nearestdx2/dx2));
	      nearestdx2 = dist2;
	      closest = n;
	    } // ENDIF check if closest	    

	  } // ENDFOR old particles

	  /* Add momentum and mass to nearest sink */

	  if (nearestdx2 < 1) {
	    
	    upold[closest] = (upold[closest] * mpold[closest] + ugrid*adddens) /
	      (mpold[closest] + adddens);
	    vpold[closest] = (vpold[closest] * mpold[closest] + vgrid*adddens) /
	      (mpold[closest] + adddens);
	    wpold[closest] = (wpold[closest] * mpold[closest] + wgrid*adddens) /
	      (mpold[closest] + adddens);
	    mpold[closest] = mpold[closest] + adddens;
	    
	    /* Record that a new particle is not needed */

	    inew = 0;

	  }  // ENDIF add to particle

	  /* Now look for nearby new sinks */

	  nearestdx2 = 1e20;
	  for (n = 0; n < ii; n++) {
	      
	    delx = xpos - xp[n];
	    dely = ypos - yp[n];
	    delz = zpos - zp[n];
	    dist2 = delx*delx + dely*dely + delz*delz;

	    /* If sink is within 5 cells, then add to it */

	    if (dist2 < 25.0*dx2 && dist2 < nearestdx2) {
//	      printf("star_maker3[addnew]: (part %"ISYM"): dx = %"GSYM", nearestdx = %"GSYM"\n",
//		     n, sqrt(dist2/dx2), sqrt(nearestdx2/dx2));
	      nearestdx2 = dist2;
	      closest = n;
	    } // ENDIF check if closest	    

	  } // ENDFOR new particles
	  
	  /* Add momentum and then mass */

	  if (nearestdx2 < 1) {

	    up[closest] = (up[closest] * mp[closest] + ugrid*adddens) /
	      (mp[closest] + adddens);
	    vp[closest] = (vp[closest] * mp[closest] + vgrid*adddens) /
	      (mp[closest] + adddens);
	    wp[closest] = (wp[closest] * mp[closest] + wgrid*adddens) /
	      (mp[closest] + adddens);
	    mp[closest] = mp[closest] + adddens;
	    
	    /* Record that a new particle is not needed */
	    
	    inew = 0;

	  } // ENDIF add to new particle

	  /* Create a new sink particle if necessary and if there's room */

	  if (inew == 1 && ii < *nmax) {

	    mp[ii] = adddens;
	    type[ii] = *ctype;
	    
	    /* Set positions and velocities */
	    
	    xp[ii] = xpos;
	    yp[ii] = ypos;
	    zp[ii] = zpos;
	    up[ii] = ugrid;
	    vp[ii] = vgrid;
	    wp[ii] = wgrid;

	    /* Set creation time */

	    tcp[ii] = (float) *t;
	    tdp[ii] = 0.0;

//	    printf("star_maker3[d]: Created new particle %"ISYM"\n", ii);
	    ii++;

	  } // ENDIF create a new sink

	} // ENDIF make sink particle

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  if (ii > 0)
    printf("sink_maker[add]: %"ISYM" new sink particles\n", ii);

//  if (nsinks > 0)
//    printf("sink_maker[sink]: %"ISYM" old sink particles\n", nsinks);

  if (ii >= *nmax) {
    fprintf(stdout, "sink_maker: reached max new particle count");
    return FAIL;
  }

  delete [] sink_index;

  *np = ii;
  return SUCCESS;

}
