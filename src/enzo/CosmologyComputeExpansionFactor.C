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
#include "macros_and_parameters.h"
#include "CosmologyParameters.h"
 
#define OMEGA_TOLERANCE 1.0e-5
 
#ifdef p4
#define ETA_TOLERANCE 1.0e-5
#endif
#ifdef p8
#define ETA_TOLERANCE 1.0e-10
#endif
#ifdef p16
#define ETA_TOLERANCE 1.0e-20
#endif
 
// function prototypes
 
double arccosh(double x);
double arcsinh(double x);
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt)
{
 
  /* Error check. */
 
  if (InitialTimeInCodeUnits == 0) {
    fprintf(stderr, "The cosmology parameters seem to be improperly set.\n");
    return FAIL;
  }
 
  *a = FLOAT_UNDEFINED;
 
  /* Find Omega due to curvature. */
 
  float OmegaCurvatureNow = 1 - OmegaMatterNow - OmegaLambdaNow;
 
  /* Convert the time from code units to Time * H0 (c.f. CosmologyGetUnits). */
 
  float TimeUnits = 2.52e17/sqrt(OmegaMatterNow)/HubbleConstantNow/
                    POW(1 + InitialRedshift,FLOAT(1.5));
 
  FLOAT TimeHubble0 = time * TimeUnits * (HubbleConstantNow*3.24e-18);
 
  /* 1) For a flat universe with OmegaMatterNow = 1, it's easy. */
 
  if (fabs(OmegaMatterNow-1) < OMEGA_TOLERANCE &&
      OmegaLambdaNow < OMEGA_TOLERANCE)
    *a      = POW(time/InitialTimeInCodeUnits, FLOAT(2.0/3.0));
 
#define INVERSE_HYPERBOLIC_EXISTS
 
#ifdef INVERSE_HYPERBOLIC_EXISTS
 
  FLOAT eta, eta_old, x;
  int i;
 
  /* 2) For OmegaMatterNow < 1 and OmegaLambdaNow == 0 see
        Peebles 1993, eq. 13-3, 13-10.
	Actually, this is a little tricky since we must solve an equation
	of the form eta - sinh(eta) + x = 0..*/
 
  if (OmegaMatterNow < 1 && OmegaLambdaNow < OMEGA_TOLERANCE) {
 
    x = 2*TimeHubble0*POW(1.0 - OmegaMatterNow, 1.5) / OmegaMatterNow;
 
    /* Compute eta in a three step process, first from a third-order
       Taylor expansion of the formula above, then use that in a fifth-order
       approximation.  Then finally, iterate on the formula itself, solving for
       eta.  This works well because parts 1 & 2 are an excellent approximation
       when x is small and part 3 converges quickly when x is large. */
 
    eta = POW(6*x, FLOAT(1.0/3.0));                     // part 1
    eta = POW(120*x/(20+eta*eta), FLOAT(1.0/3.0));      // part 2
    for (i = 0; i < 40; i++) {                          // part 3
      eta_old = eta;
      eta = arcsinh(eta + x);
      if (fabs(eta-eta_old) < ETA_TOLERANCE) break;
    }
    if (i == 40) {
      fprintf(stderr, "Case 2 -- no convergence after %"ISYM" iterations.\n", i);
      return FAIL;
    }
 
    /* Now use eta to compute the expansion factor (eq. 13-10, part 2). */
 
    *a = OmegaMatterNow/(2*(1 - OmegaMatterNow))*(cosh(eta) - 1);
    *a *= (1 + InitialRedshift);    // to convert to code units, divide by [a]
  }
 
  /* 3) For OmegaMatterNow > 1 && OmegaLambdaNow == 0, use sin/cos.
        Easy, but skip it for now. */
 
  if (OmegaMatterNow > 1 && OmegaLambdaNow < OMEGA_TOLERANCE) {
  }
 
  /* 4) For flat universe, with non-zero OmegaLambdaNow, see eq. 13-20. */
 
  if (fabs(OmegaCurvatureNow) < OMEGA_TOLERANCE &&
      OmegaLambdaNow > OMEGA_TOLERANCE) {
    *a = POW(FLOAT(OmegaMatterNow/(1 - OmegaMatterNow)), FLOAT(1.0/3.0)) *
         POW(FLOAT(sinh(1.5 * sqrt(1.0 - OmegaMatterNow)*TimeHubble0)),
	     FLOAT(2.0/3.0));
    *a *= (1 + InitialRedshift);    // to convert to code units, divide by [a]
  }
 
#endif /* INVERSE_HYPERBOLIC_EXISTS */
 
  /* Compute the derivative of the expansion factor (Peebles93, eq. 13.3). */
 
  FLOAT TempVal = (*a)/(1 + InitialRedshift);
  *dadt = sqrt( 2.0/(3.0*OmegaMatterNow*(*a)) *
	       (OmegaMatterNow + OmegaCurvatureNow*TempVal +
		OmegaLambdaNow*TempVal*TempVal*TempVal));
 
  /* Someday, we'll implement the general case... */
 
  if ((*a) == FLOAT_UNDEFINED) {
    fprintf(stderr, "Cosmology selected is not implemented.\n");
    return FAIL;
  }
 
  return SUCCESS;
}
