/***********************************************************************
/
/  GRID CLASS (SOLVE THE GADGET COOLING/HEATING RATE EQUATIONS)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:  August 02, bwo - fixed units issue (thought it was energy dens
/                     but it was really energy per unit mass)
/  modified2:  Dec. 03, bwo - fixed yet more units issues.
/  modified3:  Feb. 04, bwo - changed the way that indices are
/                     calculated, neatened up the code a bit.
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

float Gadget_EquilibriumCooling(float u_old, float rho, float dt, 
		       float *ne_guess, 
		       float *utem, float *uxyz, float *uaye, float *urho, 
		       float *utim, float redshift);

void GadgetIonizeParamsTable(float redshift);

//void GadgetIonizeParams(float redshift);

//void GadgetInitCool(float redshift);

/*
  inputs:

  d                 comoving density
  e                 total energy per unit mass
  ge                gas energy per unit mass  <- only used if dualenergyformalism is true
  u,v,w             velocities in x,y,z directions
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


int grid::GadgetCalculateCooling(float *d, float *e, float *ge, 
                 float *u, float *v, float *w,
                 int *in, int *jn, int *kn, 
                 int *iexpand, hydro_method *imethod, int *idual, int *idim,
                 int *is, int *js, int *ks, int *ie, int *je, 
                 int *ke, float *dt, float *aye,
                  float *fh, float *utem, float *uxyz, 
                 float *uaye, float *urho, float *utim,
                 float *gamma)
{
  
  int GADGETDEBUG = 0;

  if(GADGETDEBUG) printf("In GadgetCalculateCooling\n");

  /* local variables */
  int i,j,k,index;

  float u_new,u_old,  /* u_new, u_old are internal energy per unit mass*/
    totalenergypermass, gasenergypermass,
    density,
    vx,vy,vz,ne_guess,
    redshift;

  float *fne_guess;  /* guess for electron fraction */

  redshift = (1.0 + InitialRedshift)/(*aye) - 1.0;

  /* have cooling tables been generated?  if not, generate them and
     set the values of various global variables. */

  //if(HasGadgetCoolingBeenInitialized == FALSE) {
  //  GadgetInitCool(redshift);
  //  printf("GadgetCalculateCooling:  Set Cooling tables\n");
  //  HasGadgetCoolingBeenInitialized = TRUE;
  // }

  /* we want to calculate the UV background ionization
     parameters once per redshift as well!  InitCool also calls 
     this, but simplicity makes it worth a single wasted call.  */
  //GadgetIonizeParams(redshift);

  GadgetIonizeParamsTable(redshift);

  if(GADGETDEBUG) printf("GadgetCalculateCooling: loop from %"ISYM" %"ISYM" %"ISYM" to %"ISYM" %"ISYM" %"ISYM"\n",
		   (*is),(*js),(*ks),(*ie),(*je),(*ke) );
  if(GADGETDEBUG) printf("GadgetCalculateCooling: grid dims %"ISYM" %"ISYM" %"ISYM"\n",
		   (*in),(*jn),(*kn));

  fflush(stdout);

  /* set first guess for electron fraction to be 1%
     and set a pointer to it. This pointer will be
     handed to Gadget_EquilibriumCooling and will thereafter be 
     set inside the code. */

  ne_guess = 0.01;
  fne_guess = &ne_guess;

  /* loop over all cells in the grid, calculating new internal
     energy for each one. */
  for(k=(*ks); k<=(*ke);k++){
    for(j=(*js); j<=(*je);j++){
      for(i=(*is); i<=(*ie); i++){
	if(i==j && j==k && GADGETDEBUG){
          printf("Gadget: %"ISYM" %"ISYM" %"ISYM" starting loop!\n",i,j,k);
	}

	index = i + (*in)*(j + (*jn)*k);

	/* replace pointers with a variable for clarity -
	   note that these are still in CODE UNITS */
	vx=(*(u+index));
	vy=(*(v+index));
	vz=(*(w+index));
	totalenergypermass=(*(e+index));

	if( (*idual) == 1){
	  gasenergypermass=(*(ge+index));
	} else {
	  gasenergypermass = 0.0;
	}

	density=(*(d+index));
	
	if(i==j && j==k && GADGETDEBUG){
	  printf("Gadget:  i,j,k = %"ISYM" %"ISYM" %"ISYM"\n",i,j,k);
	  printf("density:  %e  totalenergypermass: %e  gasenergypermass: %e\n",
		 density,totalenergypermass,gasenergypermass);
	  printf("vx %e vy %e vz %e\n",vx,vy,vz);
	  printf("imethod: %"ISYM"  idual:  %"ISYM"\n",(*imethod),(*idual));
	  fflush(stdout);
	}

	/* calculate internal energy from energy, etc. - all in CODE UNITS */
	
	if((*imethod)==2){  /* hydromethod=2 is ZEUS hydro, 
				dual energy formalism is OFF */
	  
	  /* for whatever reason, 'total energy' in zeus hydro
	     is really 'gas energy'.  Much pain and suffering
	     will follow if you forget this. */

	  u_old = totalenergypermass;

	  //printf("new u_old from imethod=2 is:  %e\n",u_old);
	} else{  /* if it's not ZEUS, it must be PPM DE or PPM LR.  
		    Either way, we can use dual energy formalism */
	  
	  if((*idual)==1){  /* if dual energy formalism is on, we can
			    just use the gas energy */
	    u_old = gasenergypermass;
	    
	    // printf("new u_old from imethod=1 or 0, idual=1 is %e\n",u_old);
	  } else{  /* otherwise it's the same as for ZEUS hydro */
	    
	    u_old = (totalenergypermass 
		     - 0.5 * (vx*vx + vy*vy + vz*vz));
	    //printf("new u_old from imethod=1 or 0, idual=0 is %e\n",u_old);
	  }  /* end of if(idual  */
	}  /* end of if(hydro... */
	
	/* get new internal energy per unit mass -- make sure to have debug
	   flag.  Note that everything passed in is in code units, and everything
	   that comes out should ALSO be in code units.  */

	u_new=Gadget_EquilibriumCooling(u_old,density,(*dt),fne_guess,
				       utem,uxyz,uaye,urho,utim,redshift);
	/* error checking:  If something's screwed up, Gadget_EquilibriumCooling
	   gives us a negative energy density, so return an error */
	
	if(u_new < 0.0){
	  ENZO_VFAIL("GadgetCalculateCooling:  incorrect internal energy calculated:  %e  Exiting.\n",u_new)
	}	
	
	/* store the new energy value - basically reverse the process above.
	   All of this is in CODE UNITS */
	
	
	if((*imethod)==2){  /* hydromethod=2 is ZEUS hydro, 
				dual energy formalism is OFF */

	  /* see above comment - 'total energy' is really 'gas
	     energy' in ZEUS hydro */

	  (*(e+index)) = u_new;
	  
	} else{  /* if it's not ZEUS, it must be PPM DE or PPM LR.  
		    Either way, we can use dual energy formalism */
	  
	  if((*idual)==1){  /* if dual energy formalism is on, we can
			    just use the gas energy - but we must also 
			    set total energy, otherwise energy isn't 
			    conserved.  Which is bad, and does funny
			    things.  */
	    
	    (*(ge+index)) = u_new;

	    (*(e+index)) = (u_new + 0.5 * 
						  (vx*vx + vy*vy + vz*vz));
	    
	  } else{  /* otherwise it's the same as for ZEUS hydro */
	    
	    (*(e+index)) = (u_new + 0.5 * 
						  (vx*vx + vy*vy + vz*vz));
	    
	  }  /* end of if(idual  */
	}  /* end of if(hydro... */

	if(i==j && j==k && GADGETDEBUG){

	  printf("Gadget: %"ISYM" %"ISYM" %"ISYM" finishing loop!\n",i,j,k);
	}
	
      }  /* for(i=is; i<=ie; i++) */
    } /* for(j=js; j<=je;j++) */
  } /*  for(k=ks; k<=ke;k++) */
  
  return SUCCESS;
}

