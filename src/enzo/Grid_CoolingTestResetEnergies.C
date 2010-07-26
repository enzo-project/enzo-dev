/***********************************************************************
/
/  GRID CLASS (Adjust internal energies for cooling test problem.)
/
/  written by: Britton Smith
/  date:       February, 2008
/  modified1:
/
/  PURPOSE: Adjust internal energies to keep temperature constant as mu changes.
/
/  RETURNS: FAIL or SUCCESS
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

#define MH 1.67e-24
#define DEFAULT_MU 0.6

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::CoolingTestResetEnergies()
{

  /* Return if this doesn't concern us. */

  if ((ProblemType != 62) || !(TestProblemData.ResetEnergies)) {
    return SUCCESS;
  }

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  // Mu never varies for MultiSpecies = 0, so nothing to do.
  if (MultiSpecies < 1) {
    return SUCCESS;
  }

  /* declarations */

  int i, j, k, index;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum;

  float temperatureSlope;
  float mu;

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  /* Find fields: density, total energy, velocity1-3. */
 
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  /* Find Multi-species fields. */

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      return FAIL;
    }

  temperatureSlope = log10(TestProblemData.MaximumTemperature / TestProblemData.MinimumTemperature) /
    (GridEndIndex[2] - GridStartIndex[2]);

  // Reset internal energies from temperature.

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) { // Temperature
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) { // Metallicity
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) { // H NumberDensity

	index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];

	// calculate mu

	mu = BaryonField[DeNum][index] + BaryonField[HINum][index] + BaryonField[HIINum][index] + 
	  (BaryonField[HeINum][index] + BaryonField[HeIINum][index] + BaryonField[HeIIINum][index])/4.0;
	if (MultiSpecies > 1) {
	  mu += BaryonField[HMNum][index] + (BaryonField[H2INum][index] + BaryonField[H2IINum][index])/2.0;
	}
	if (MultiSpecies > 2) {
	  mu += (BaryonField[DINum][index] + BaryonField[DIINum][index])/2.0 + (BaryonField[HDINum][index]/3.0);
	}
	mu = BaryonField[DensNum][index] / mu;

	BaryonField[TENum][index] = pow(10,((temperatureSlope * (k-GridStartIndex[2])) + log10(TestProblemData.MinimumTemperature))) /
	  TemperatureUnits / mu / (Gamma-1.0);

	if (DualEnergyFormalism)
	  BaryonField[GENum][index] = BaryonField[TENum][index];

      }
    }
  }

  return SUCCESS;
}
