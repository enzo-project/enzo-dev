/***********************************************************************
/
/  GRID CLASS (COMPUTE THE RATIO OF PRESSURE TO GRAVITY FOR ONE-ZONE
/              COLLAPSE)
/
/  written by: Britton Smith
/  date:       December, 2012
/  modified1:
/
/  PURPOSE: Computes an approximate ratio of pressure to gravitational 
/           force as modifier to the freefall collapse solution.
/           This follows equations 6-9 of Omukai et al (2005).
/
/  RETURNS:
/
************************************************************************/
 
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
 
int grid::ComputeOneZoneCollapseFactor(float *force_factor)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;  

  /* Compute the size of the fields. */
 
  int i, j, k, t, der, index, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];


  /* Initialize values. */

  for (i = 0; i < size; i++) {
    force_factor[i] = 0.0;
  }

  if (!TestProblemData.OneZoneFreefallUseEffectiveGamma) {
    return SUCCESS;
  }
 
  /* Check for density and pressure history. */
  der = 1;
  if (freefall_density[1] == NULL)
    return SUCCESS;
  if (freefall_density[2] != NULL)
    der = 2;

  /* Find Density, if possible. */
 
  int DensNum;  
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0)
    ENZO_FAIL("Cannot find density.");

  float gamma_eff;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) { // nothing
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) { // metallicity
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) { // energy

	index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];

	/* Calculate the effective adiabatic index, dlog(p)/dlog(rho). */
        gamma_eff = log10(freefall_pressure[0][index] /
                          freefall_pressure[1][index]) /
          log10(freefall_density[0][index] /
                freefall_density[1][index]);
        if (der > 1) {
          gamma_eff += 0.5 * ((log10(freefall_pressure[1][index] /
                                     freefall_pressure[2][index]) /
                               log10(freefall_density[1][index] /
                                     freefall_density[2][index])) - gamma_eff);
        }
        gamma_eff = min(gamma_eff, (4./3.));

	if (gamma_eff < 0.83) {
	  force_factor[index] = 0.0;
	}
	else if (gamma_eff < 1.0) {
	  force_factor[index] = 0.6 + 2.5 * (gamma_eff - 1) -
	    6.0 * POW((gamma_eff - 1.0), 2.);
	}
	else {
	  force_factor[index] = 1.0 + 0.2 * (gamma_eff - (4./3.)) -
	    2.9 * POW((gamma_eff - (4./3.)), 2.);
	}
	force_factor[index] = max(force_factor[index], 0.0);
	force_factor[index] = min(force_factor[index], 0.95);

      }
    }
  }
 
  return SUCCESS;
}
