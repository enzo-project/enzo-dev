/***********************************************************************
/
/  COSMOLOGY: COMPUTES THE TIME (IN CODE UNITS) FROM THE GIVEN REDSHIFT
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

int CosmologyTableComputeTimeFromRedshift(FLOAT z, FLOAT *time);


int CosmologyComputeTimeFromRedshift(FLOAT Redshift, FLOAT *TimeCodeUnits)
{
 
  FLOAT TimeHubble0;
 
  /* Find Omega due to curvature. */
 
  float OmegaCurvatureNow = 1 - OmegaMatterNow -
    OmegaLambdaNow - OmegaRadiationNow;

  /* Interpolate from a(t) table. */
 
  if (CosmologyTableComputeTimeFromRedshift(Redshift, &TimeHubble0) == FAIL) {
    ENZO_FAIL("Error in CosmologyTableComputeTime.\n");
  }

  /* Now convert from Time * H0 to code units (see also CosmologyGetUnits). */
 
  float TimeUnits = 2.52e17/sqrt(OmegaMatterNow)/HubbleConstantNow/
                    POW(1 + InitialRedshift, FLOAT(1.5));
 
  *TimeCodeUnits = TimeHubble0 / (HubbleConstantNow*3.24e-18) / TimeUnits;
 
 
  return SUCCESS;
}
