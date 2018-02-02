/***********************************************************************
/
/  COMPUTES THE EXPANSION FACTORS (A & DADT) AT THE REQUESTED TIME
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "CosmologyParameters.h"
 
// function prototypes
 
int CosmologyTableComputeExpansionFactor(FLOAT time, FLOAT *a);

 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt)
{
 
  /* Error check. */
 
  if (InitialTimeInCodeUnits == 0) {
    ENZO_FAIL("The cosmology parameters seem to be improperly set.\n");
  }
 
  /* Find Omega due to curvature. */
 
  float OmegaCurvatureNow = 1 - OmegaMatterNow -
    OmegaLambdaNow - OmegaRadiationNow;
 
  /* Convert the time from code units to Time * H0 (c.f. CosmologyGetUnits). */
 
  float TimeUnits = 2.52e17/sqrt(OmegaMatterNow)/HubbleConstantNow/
                    POW(1 + InitialRedshift,FLOAT(1.5));
 
  FLOAT TimeHubble0 = time * TimeUnits * (HubbleConstantNow*3.24e-18);

  /* Interpolate from a(t) table. */
 
   if (CosmologyTableComputeExpansionFactor(TimeHubble0, a) == FAIL) {
      ENZO_FAIL("Error in CosmologyTableComputeExpansionFactor.\n");
  }
 
  /* Compute the derivative of the expansion factor (Peebles93, eq. 13.3). */
 
  FLOAT TempVal = (*a)/(1 + InitialRedshift);
  *dadt = sqrt( 2.0/(3.0*OmegaMatterNow*(*a)) *
	       (OmegaMatterNow + OmegaCurvatureNow*TempVal +
		OmegaLambdaNow*TempVal*TempVal*TempVal +
                OmegaRadiationNow/TempVal));
 
  return SUCCESS;
}
