/***********************************************************************
/
/  GRID CLASS (SOLVE THE ANALYTICAL SOLUTION FOR FREE-FALL COLLAPSE)
/
/  written by: Britton Smith
/  date:       October, 2010
/  modified1:  
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
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
#include "fortran.def"
#include "CosmologyParameters.h"
#include "phys_constants.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

int grid::SolveOneZoneFreefall()
{

  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  this->DebugCheck("SolveRadiativeCooling");

  /* Declarations */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  FLOAT a = 1.0, dadt;
    
  /* Find fields: density, total energy, velocity1-3. */

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Find Multi-species fields. */

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  /* Compute size of the current grid. */

  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  /* If using cosmology, compute the expansion factor and get units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  if (ComovingCoordinates) {

    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	== FAIL) {
            ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
    }

    aUnits = 1.0/(1.0 + InitialRedshift);

  }

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }

  /* Metal cooling codes. */

  int MetalNum = 0, SNColourNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

  /* Calculate new density and energy. */

  float FreefallTimeConstant = POW(((32 * GravitationalConstant) / (3 * pi)), 0.5);
  float NewDensity = POW((TestProblemData.OneZoneFreefallConstant - 
			  (0.5 * FreefallTimeConstant * Time)), -2.);
  float DensityRatio = NewDensity / BaryonField[DensNum][0];

  /* Update all cells. */

  for (i = 0;i < size;i++) {

    /* Update enegy. */

    BaryonField[TENum][i] += (Gamma - 1) * BaryonField[TENum][i] * 
      FreefallTimeConstant * POW(BaryonField[DensNum][i], 0.5) * dtFixed;
    if (DualEnergyFormalism) {
      BaryonField[GENum][i] = BaryonField[TENum][i];
    }

    /* Update density. */

    BaryonField[DensNum][i] = NewDensity;

    /* Update species fields. */

    if (MultiSpecies) {
      BaryonField[DeNum][i] *= DensityRatio;
      BaryonField[HINum][i] *= DensityRatio;
      BaryonField[HIINum][i] *= DensityRatio;
      BaryonField[HeINum][i] *= DensityRatio;
      BaryonField[HeIINum][i] *= DensityRatio;
      BaryonField[HeIIINum][i] *= DensityRatio;
      if (MultiSpecies > 1) {
	BaryonField[HMNum][i] *= DensityRatio;
	BaryonField[H2INum][i] *= DensityRatio;
	BaryonField[H2IINum][i] *= DensityRatio;
      }
      if (MultiSpecies > 2) {
	BaryonField[DINum][i] *= DensityRatio;
	BaryonField[DIINum][i] *= DensityRatio;
	BaryonField[HDINum][i] *= DensityRatio;
      }
    }

    if (MetalFieldPresent) {
      BaryonField[MetalNum][i] *= DensityRatio;
      if (MultiMetals) {
	BaryonField[MetalNum+1][i] *= DensityRatio;
	BaryonField[MetalNum+2][i] *= DensityRatio;
      }
    }
  }

  return SUCCESS;

}
