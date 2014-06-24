/***********************************************************************
/
/  GRID CLASS (HANDLE THE EMISSIVITY FROM STAR PARTICLE FORMATION)
/
/  written by: Geoffrey So
/  date:       December, 2008
/  modified1:  June 10, 2009 Geoffrey So
/                fixed some comments when writing the Wiki page
/
/  PURPOSE: The emissivity is calculated for a region with star 
/           formation, then the value is fed into the radiation solver 
/           to calculate how well the gray radiation energy will 
/           propagate.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/


#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "StarParticleData.h"


/* Translating everything from star_maker2.src feedback section including paraphrasing the comments from that section */

int CalcEmiss(int *nx, int *ny, int *nz,
             float *d, float *dm, float *te, float *ge, float *u, float *v,
		       float *w, float *metal,
             int *idual, int *imetal, hydro_method *imethod, float *dt,
		       float *r, float *dx, FLOAT *t, float *z,
             float *d1, float *x1, float *v1, float *t1,
                       float *sn_param, float *m_eject, float *yield,
             int *nmax, FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
		       int *ibuff,
             FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up, float *vp, float *wp,
             float *mp, float *tdp, float *tcp, float *metalf,
	      float *justburn, float *EmissivityArray, float dtLevelAbove) {

  int i, j, k, n;
  float mform, tfactor, uv_energy, 
        minitial, xv1, xv2, dratio;
  float clight = 3e10;
  float msolar_e51 = 1800.0;

  /* uv_param taken from upper limits of
     Razoumov and Norman 2002 June 20 The Astrophysical Journal */
  //float uv_param = 1e-5;
  //printf("UV_PARAM IS %"FSYM" \n", uv_param);

  /* disabling clear of Emissivity field until a way to do it in AMR is found 
     in EvolveHierarchy and EvolveLevel.  Instead put the clear here manually
   */
  n = 0;
  int ArraySize = *nx * *ny * *nz;
  while (n < ArraySize) {
    EmissivityArray[n]=0;
    n++; 
  }

  // check initial EmissivityArray values 

  /*  
  int ArraySize = *nx * *ny * *nz;
  int ArrayCount;

  n = 0;

  while (n < ArraySize) {
    if(EmissivityArray[n]==0)
      n++;
    else {
      printf("AN ELEMENT IS %"FSYM" at %i \n", EmissivityArray[n], n);
      printf("dtFixed dtLevelAbove %"FSYM" %"FSYM" \n", *dt, dtLevelAbove);
      n = ArraySize+1;
    }
  }
  if(n == ArraySize)
    printf("ARRAY size %"ISYM" initially ZEROs, n at %i \n", ArraySize, n);
  */


  /* Loop over each star particle that was created to calculate 
     contribution to Emissivity */
  for (n = 0; n < *nmax; n++) {
    //printf("for loop over n = %"ISYM" \n",n);

    /* check if conditions are right for each particle */
    /* tcp is < 0 for non star particle, and mass has to be positive */
    if ( (tcp[n] > 0) && (mp[n] > 0) ) {
      /* next time if not seeing stars being made in fortran, check if 
	 feedback is turned on! */
      //printf("if loop 1\n");


      /* Comments:
	 The star particle creation algorithm partnered with this 
	 feedback algorithm creates a star particle instantaneously.
	 However, we do feedback as if the star particles are created 
	 over a long period of time (in code units), so the particle
	 actually loses mass over time in an exponentially decaying
	 way.

	 Determine how much of a given star particle would have been 
	 turned into stars during this timestep.  Then calculate the mass
	 which should have formed during this timestel dt using the integral
	 form of the Cen & Ostriker formula.
      */
      
      xv1 = (*t       - tcp[n])/tdp[n];

      /* if t-tcp >> tdp then ignore */
      if (xv1 <= 12.0) {
	//printf("if loop 2\n");

	xv2 = (*t + *dt - tcp[n])/tdp[n];

	/* Calculate the initial mass of the star particle being looped over */
	minitial = mp[n]/
	  (1.0 - *m_eject*(1.0 - (1.0 + xv1)*exp(-xv1)));

	/* Calculate the amount of mass that would have formed star in 
	   this timestep */
	mform = minitial * ((1.0 + xv1)*exp(-xv1) - (1.0 + xv2)*exp(-xv2));
	//mp is particle mass
	mform = max(min(mform, mp[n]), 0.0);

	/* Calculate the cell index of the star particle */
	/* 
	   This differs from Fortran counter part because the arrays 
	   in C start with index 0 instead of 1 in Fortran, so they 
	   are all missing the + 1 at the end 
	*/
	i = int((xp[n] - *xstart)/ *dx);
	j = int((yp[n] - *ystart)/ *dx);
	k = int((zp[n] - *zstart)/ *dx);

	/* 
	   Check if the star is within the bounds of the grid
	   If inside then calc star formation energy and update 
	   Emissivity or else give warning 
	*/
	if ( (i >= 0) && (i < *nx) &&
	     (j >= 0) && (j < *ny) &&
	     (k >= 0) && (k < *nz) ) {	  
	  //printf("if loop 3\n");
	  
	  /* 
	     Skip star formation and emissivity update if too 
	     little mass is formed 
	  */
	  if (mform/d[i,j,k] >= 1.0e-10) {
	    //printf("if loop 4\n");

	    /* Subtract ejected mass from particle (due to winds, supernovae) */
	    mp[n] = mp[n] - mform * *m_eject;

	    /* 
	       Calculate how much of the star formation energy in this
	       timestep would have gone into UV radiation energy 
	    */

	    // the following commented out is the old star_maker2 energy
	    /* 
	    uv_energy = uv_param * mform * pow(clight/ *v1,2) /
	      (d[i,j,k] + mform * *m_eject);
	    */

	    /* convert energy density in Enzo units to energy rate */
	    /* multiply by density*velocity^2/time conversion factors */
	    /*
	    uv_energy = uv_param * mform * (clight/ *v1) * (clight/ *v1) *
	      (*d1) * (*v1) * (*v1) / (dtLevelAbove * *t1);
	    */

	    /* noticed that the v1 cancels out from above, give the following */
            /* Units of uv_energy will be in erg/s/cm^3 */

	    uv_energy = uv_param * mform * (clight) * (clight) *
	      *d1 / (dtLevelAbove * *t1);


	    /* update the Emissivity_Array with <UV energy weighted by time> */
	    /*printf("before %i %i %i was %22.16e\n", EmissivityArray[i + *nx * (j + *ny * k)],i,j,k);*/
	    EmissivityArray[i + *nx * (j + *ny * k)] += uv_energy;
	    /*printf("after %i %i %i is %22.16e\n", EmissivityArray[i + *nx * (j + *ny * k)],i,j,k);*/
	  }
	}
	else {
	  printf("warning star particle out of grid in C %"ISYM", %"ISYM", %"ISYM" \n", i, j, k);
	}

      }
      /*
      else
	printf("NO STARS MADE \n");
      */
    }//10 in Fortran
    /*
    else
      printf("NO STARS MADE \n");
    */
  }//100 in Fortran

  return SUCCESS;
}
