/***********************************************************************
/
/  GRID CLASS (SOLVE THE GADGET EQUILIBRIUM COOLING)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:  calculate the cooling rate de/dt using gadget equilibrium
/        cooling and return de/dt
/
/  RETURNS:
/    de/dt  (the cooling rate in CGS units)
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

float Gadgetconvert_u_to_temp(float u, float rho, float *ne_guess);
float GadgetCoolingRate(float logT, float rho, float *nelec, float redshift);

float grid::GadgetCoolingRateFromU(float u, float rho, float *ne_guess,float redshift)
{
  float temp;

  //if(debug) printf("In GadgetCoolingRateFromU:1\n");

  temp = Gadgetconvert_u_to_temp(u, rho, ne_guess);

  //if(debug) printf("In GadgetCoolingRateFromU:2\n");

  return GadgetCoolingRate(log10(temp), rho, ne_guess,redshift);
}
