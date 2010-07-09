/***********************************************************************
/
/  GRID CLASS (SOLVE THE GADGET COOLING/HEATING RATE EQUATIONS)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:  August 02, bwo - fixed units issue (thought it was energy dens
/                     but it was really energy per unit mass)
/
/  PURPOSE:
/      This routine loops over all of the cells in a grid, calculating
/      the new internal energy per unit mass of each cell using the GADGET
/      equilibrium cooling model described in Katz et. al. 1996 
/      (ApJS 105; 19-35)
/
/  RETURNS:
/    SUCCESS or FAIL
/
/  NOTE:  This routine is based on code in the GADGET code, written
/    by Volker Springel and Naoki Yoshida.  The subroutines are
/    essentially identical to the Gadget cooling routines, with an
/    Enzo-friendly wrapper around them.  As of the addition of this
/    comment to the code (Feb. 2004) the cooling code is private and
/    NOT FOR PUBLIC RELEASE until Volker makes the cooling code public
/    in Gadget.  See the Gadget web page for more details:
/
/      http://www.mpa-garching.mpg.de/~volker/gadget/
/
/      --B.W.O'Shea, Feb. 2004
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
#include "CosmologyParameters.h"
#include "Gadget.h"

void GadgetIonizeParamsTable(float redshift);

float GadgetCoolingRateFromU(float u, float rho, float *ne_guess, float redshift);

/*
  inputs:

  d                 comoving density
  e                 total energy per unit mass
  ge                gas energy per unit mass  <- only used if dualenergyformalism is true
  u,v,w             velocities in x,y,z directions
  cooltime          cooling time
  in,jn,kn          grid dimensions
  iexpand           Comoving coordinates flag
  imethod           hydro method flag
  idual             dual energy formalism flag
  idim              grid rank (we want 3-d!)
  is,js,ks          grid start
  ie,je,ke          grid end
  dt                timestep
  aye               a
  fh,               hydrogen fraction by mass
  utem,uxyz,uaye,   temp, length, a units
  urho,utim         density, time units
  gamma             gamma

 */


int grid::GadgetCoolingTime(float *d, float *e, float *ge, 
                 float *u, float *v, float *w,
		   float *cooltime,
                 int *in, int *jn, int *kn, 
                 int *iexpand, hydro_method *imethod, int *idual, int *idim,
                 int *is, int *js, int *ks, int *ie, int *je, 
                 int *ke, float *dt, float *aye,
                  float *fh, float *utem, float *uxyz, 
                 float *uaye, float *urho, float *utim,
                 float *gamma)
{
  
  if(debug) printf("In GadgetCoolingTime\n");

  /* local variables */
  int i,j,k;

  float u_new,u_old,  /* u_new, u_old are internal energy per unit mass*/
    totalenergypermass, gasenergypermass,
    density,
    vx,vy,vz,ne_guess,
    redshift,energy, ratefact,LambdaNet, coolingtime;

  float *fne_guess;  /* guess for electron fraction */

  redshift = (1.0 + InitialRedshift)/(*aye) - 1.0;

  /* we want to calculate the UV background ionization
     parameters once per redshift as well!  InitCool also calls 
     this, but simplicity makes it worth a single wasted call.  */

  GadgetIonizeParamsTable(redshift);

  if(debug) printf("GadgetCoolingTime: loop from %"ISYM" %"ISYM" %"ISYM" to %"ISYM" %"ISYM" %"ISYM"\n",
		   (*is),(*js),(*ks),(*ie),(*je),(*ke) );
  if(debug) printf("GadgetCoolingTime: grid dims %"ISYM" %"ISYM" %"ISYM"\n",
		   (*in),(*jn),(*kn));


  /* set first guess for electron fraction to be 1% and set a pointer to it. */
  ne_guess = 0.01;
  fne_guess = &ne_guess;
  

  /* loop over all cells in the grid, calculating new internal
     energy for each one. */
  for(k=(*ks); k<=(*ke);k++){
    for(j=(*js); j<=(*je);j++){
      for(i=(*is); i<=(*ie); i++){

	/* replace pointers with a variable for clarity -
	   note that these are still in CODE UNITS */

	vx=(*(u+k*(*in)*(*jn) + j*(*in) + i));
	vy=(*(v+k*(*in)*(*jn) + j*(*in) + i));
	vz=(*(w+k*(*in)*(*jn) + j*(*in) + i));
	totalenergypermass=(*(e+k*(*in)*(*jn) + j*(*in) + i));
	density=(*(d+k*(*in)*(*jn) + j*(*in) + i));
	
	/* calculate internal energy from energy, etc. - all in CODE UNITS */
	
	if((*imethod)==2){  /* hydromethod=2 is ZEUS hydro, 
				dual energy formalism is OFF */
	  
	  /* for whatever reason, 'total energy' in zeus hydro
	     is really 'gas energy'.  Much pain and suffering
	     will follow if you forget this. */

	  u_old = totalenergypermass;

	} else{  /* if it's not ZEUS, it must be PPM DE or PPM LR.  
		    Either way, we can use dual energy formalism */
	  
	  if((*idual)==1){  /* if dual energy formalism is on, we can
			    just use the gas energy */
	    
	    gasenergypermass=(*(ge+k*(*in)*(*jn) + j*(*in) + i));
	    u_old = gasenergypermass;
	    
	    // printf("new u_old from imethod=1 or 0, idual=1 is %e\n",u_old);
	  } else{  /* otherwise it's the same as for ZEUS hydro */
	    
	    u_old = (totalenergypermass 
		     - 0.5 * (vx*vx + vy*vy + vz*vz));
	    //printf("new u_old from imethod=1 or 0, idual=0 is %e\n",u_old);
	  }  /* end of if(idual  */
	}  /* end of if(hydro... */
	

	// everything in physical units!
	density *= (*urho);  // baryon density in proper CGS!
	nHcgs = XH * density / PROTONMASS;  // hydrogen number density in proper CGS
	ratefact = nHcgs * nHcgs / density;

	// internal energy per unit mass in proper CGS
	energy = u_old * POW( ( (*uxyz)*(1.0+redshift)/(*utim)/(1.0+InitialRedshift)),2.0);


	LambdaNet = GadgetCoolingRateFromU(energy, density, fne_guess, redshift);

	if(LambdaNet >= 0){
	  coolingtime = 0;
	} else {

	  coolingtime = u_old / (-ratefact * LambdaNet);

	}

	// cooling time back in code units
	(*(cooltime+k*(*in)*(*jn) + j*(*in) + i)) = coolingtime/(*utim);

      }  /* for(i=is; i<=ie; i++) */
    } /* for(j=js; j<=je;j++) */
  } /*  for(k=ks; k<=ke;k++) */
  
  return SUCCESS;
}

