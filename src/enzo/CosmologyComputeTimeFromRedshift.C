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
 
double arccosh(double x);
double arcsinh(double x);
 
int CosmologyComputeTimeFromRedshift(FLOAT Redshift, FLOAT *TimeCodeUnits)
{
 
  /* Error check [INCOMPLETE]. */
 
  FLOAT eta, TimeHubble0 = FLOAT_UNDEFINED;
 
  /* Find Omega due to curvature. */
 
  float OmegaCurvatureNow = 1 - OmegaMatterNow - OmegaLambdaNow;
 
  /* 1) For a flat universe with OmegaMatterNow = 1, things are easy. */
 
  if (OmegaMatterNow == 1 && OmegaLambdaNow == 0)
    TimeHubble0 = 2.0/3.0/POW(1+Redshift, FLOAT(1.5));
 
#define INVERSE_HYPERBOLIC_EXISTS
 
#ifdef INVERSE_HYPERBOLIC_EXISTS
 
  /* 2) For OmegaMatterNow < 1 and OmegaLambdaNow == 0 see
        Peebles 1993, eq. 13-3, 13-10. */
 
  if (OmegaMatterNow < 1 && OmegaLambdaNow == 0) {
    eta = arccosh(1 + 2*(1-OmegaMatterNow)/OmegaMatterNow/(1+Redshift));
    TimeHubble0 = OmegaMatterNow/(2*POW(1.0-OmegaMatterNow, 1.5))*
		  (sinh(eta) - eta);
  }
 
  /* 3) For OmegaMatterNow > 1 && OmegaLambdaNow == 0, use sin/cos. */
 
  if (OmegaMatterNow > 1 && OmegaLambdaNow == 0) {
    eta = acos(1 - 2*(1-OmegaMatterNow)/OmegaMatterNow/(1+Redshift));
    TimeHubble0 = OmegaMatterNow/(2*POW(1.0-OmegaMatterNow, 1.5))*
                  (eta - sin(eta));
  }
 
  /* 4) For flat universe, with non-zero OmegaLambdaNow, see eq. 13-20. */
 
  if (fabs(OmegaCurvatureNow) < 1.0e-3 && OmegaLambdaNow != 0)
    TimeHubble0 = 2.0/3.0/sqrt(1-OmegaMatterNow)*
		    arcsinh(sqrt((1-OmegaMatterNow)/OmegaMatterNow)/
		           POW(1+Redshift,FLOAT(1.5))             );
 
#endif /* INVERSE_HYPERBOLIC_EXISTS */
 
  /* 5) Someday, we'll implement the general case... */
 
  if (TimeHubble0 == FLOAT_UNDEFINED) {
    ENZO_FAIL("Cosmology selected is not implemented.\n");

  }
 
  /* Now convert from Time * H0 to code units (see also CosmologyGetUnits). */
 
  float TimeUnits = 2.52e17/sqrt(OmegaMatterNow)/HubbleConstantNow/
                    POW(1 + InitialRedshift, FLOAT(1.5));
 
  *TimeCodeUnits = TimeHubble0 / (HubbleConstantNow*3.24e-18) / TimeUnits;
 
 
  return SUCCESS;
}
