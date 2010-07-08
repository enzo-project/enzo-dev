/***********************************************************************
/
/  GRID CLASS (ACTUALLY SOLVE THE GADGET EQUILIBRIUM COOLING)
/
/  written by: Brian O'Shea
/  date:  June 2002
/  modified1: August 2002, bwo.  Fixed units problem.  (formerly energy/volume,
/                                now energy/unit mass)
/
/  PURPOSE:  This routine uses the GADGET equilibrium cooling model 
/      described in Katz et. al. 1996 (ApJS 105; 19-35) to calcualte
/      the new internal energy density of a cell, given the old
/      internal energy density, density, timestep, and a guess at
/      the electron number density.  All of these are given to the
/      routine in CODE UNITS
/            
/  RETURNS:
/    the new internal energy per unit mass, in code units.
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

float GadgetCoolingRateFromU(float u, float rho, float *ne_guess, float redshift);

float grid::Gadget_EquilibriumCooling(float u_old, float rho, float dt, 
		       float *ne_guess, 
		       float *utem, float *uxyz, float *uaye, float *urho, 
		       float *utim, float redshift)
{

  //if(debug) printf("in Gadget_EquilibriumCooling\n");

  float u, du;
  float u_lower, u_upper;
  float ratefact;
  float LambdaNet;
  int iter = 0;


 //if(debug) printf("u_old is: %e %e\n",u_old,u_old);

  DoCool_u_old_input = u_old;
  DoCool_rho_input = rho;
  DoCool_dt_input = dt;
  DoCool_ne_guess_input = *ne_guess;

  /* convert to physical cgs units */
  rho *= (*urho);

  /* energy units are energy/unit mass, which correspons to velocity_units^2 
     see:  http://zeus.ncsa.uiuc.edu/~gbryan/amr_guide/output.html
     for units.
     This is now physical CGS units!  */
  u_old *= POW( ( (*uxyz)*(1.0+redshift)/(*utim)/(1.0+InitialRedshift)),2.0);

  //if(debug) printf("in Gadget_EquilibriumCooling:2\n");

  // time in CGS
  dt *= (*utim);

  /* in gadget:
     rho *= All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
     u_old *= All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
     dt *= All.UnitTime_in_s / All.HubbleParam; 
     don't erase this until we KNOW it's useless */

  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
  ratefact = nHcgs * nHcgs / rho;

  u = u_old;
  u_lower = u;
  u_upper = u;

  //if(debug) printf("in Gadget_EquilibriumCooling:3\n");

  LambdaNet = GadgetCoolingRateFromU(u, rho, ne_guess,redshift);

  /* bracketing */
  //if(debug) printf("in Gadget_EquilibriumCooling:4\n");

  if(u - u_old - ratefact * LambdaNet * dt < 0)	/* heating */
    {
      u_upper *= sqrt(1.1);
      u_lower /= sqrt(1.1);
      while(u_upper - u_old - ratefact * GadgetCoolingRateFromU(u_upper, rho, ne_guess,redshift) * dt < 0)
	{
	  u_upper *= 1.1;
	  u_lower *= 1.1;
	}
    }

  if(u - u_old - ratefact * LambdaNet * dt > 0)
    {
      u_lower /= sqrt(1.1);
      u_upper *= sqrt(1.1);
      while(u_lower - u_old - ratefact * GadgetCoolingRateFromU(u_lower, rho, ne_guess,redshift) * dt > 0)
	{
	  u_upper /= 1.1;
	  u_lower /= 1.1;
	}
    }
  //if(debug) printf("in Gadget_EquilibriumCooling:5\n");

  do
    {
      u = 0.5 * (u_lower + u_upper);

      LambdaNet = GadgetCoolingRateFromU(u, rho, ne_guess,redshift);

      if(u - u_old - ratefact * LambdaNet * dt > 0)
	{
	  u_upper = u;
	}
      else
	{
	  u_lower = u;
	}

      du = u_upper - u_lower;

      iter++;

      if(iter >= (MAXITER - 10))
	printf("u= %e\n", u);
    }
  while(fabs(du / u) > 1.0e-6 && iter < MAXITER);

  //if(debug) printf("in Gadget_EquilibriumCooling:6\n");

  if(iter >= MAXITER)
    {
      printf("failed to converge in DoCooling()\n");
      printf("DoCool_u_old_input=%e\nDoCool_rho_input= %e\nDoCool_dt_input= %e\nDoCool_ne_guess_input= %e\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      u = -1.0;  /* if we screw up return negative energy density so the code knows to shut down */
      return u;
    }

  /* u *= All.UnitDensity_in_cgs / All.UnitPressure_in_cgs; */	/* to internal units */

  /* back into code units */
  u /= POW( ( (*uxyz)*(1.0+redshift)/(*utim)/(1.0+InitialRedshift)),2.0);

  //if(debug) printf("exiting Gadget_EquilibriumCooling\n");

  return u; /* make sure this is in CODE UNITS! */

}
