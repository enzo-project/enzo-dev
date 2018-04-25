/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID TO A UNIFORM POOL OF GAS)
/
/  written by: Britton Smith
/  date:       February, 2008
/  modified1:
/
/  PURPOSE: Only different from InitializeUniformGrid in how 
/           species fractions are initialized.
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

int grid::CoolingTestInitializeGrid()
{
  /* declarations */

  int dim, i, j, k, size, field, GCM, index;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, MetalNum;

  int ExtraField[2];

  float HNumberDensity, HNumberDensitySlope, metallicitySlope, temperatureSlope;
  float mu;

  /* create fields */
 
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1)
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  int colorfields = NumberOfBaryonFields;

  // Enzo's standard multispecies (primordial chemistry - H, D, He)
  if (TestProblemData.MultiSpecies) {
    FieldType[DeNum     = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum     = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum    = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum    = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum   = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum  = NumberOfBaryonFields++] = HeIIIDensity;
    if (TestProblemData.MultiSpecies > 1) {
      FieldType[HMNum   = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum  = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (TestProblemData.MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }

  //  Metal fields, including the standard 'metallicity' as well 
  // as two extra fields
  if (TestProblemData.UseMetallicityField) {
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity;

    if(TestProblemData.MultiMetals){
      FieldType[ExtraField[0] = NumberOfBaryonFields++] = ExtraType0;
      FieldType[ExtraField[1] = NumberOfBaryonFields++] = ExtraType1;
    }
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* compute size of fields */

  if (GridRank < 3) {
    fprintf(stderr, "GridRank must be 3 for Cooling Test.\n");
    return FAIL;
  }

  // Get the units so we can convert temperature later
 
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits,
	       Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* allocate fields */

  this->AllocateGrids();

  // Set up a 3d grid that varies over density, metallicity, and temperature.
  // Do density and metallicity first, and temperature later when species fractions are set.

  HNumberDensitySlope = log10(TestProblemData.MaximumHNumberDensity / TestProblemData.MinimumHNumberDensity) /
    (GridEndIndex[0] - GridStartIndex[0]);
  metallicitySlope = log10(TestProblemData.MaximumMetallicity / TestProblemData.MinimumMetallicity) /
    (GridEndIndex[1] - GridStartIndex[1]);
  temperatureSlope = log10(TestProblemData.MaximumTemperature / TestProblemData.MinimumTemperature) /
    (GridEndIndex[2] - GridStartIndex[2]);

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) { // Temperature
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) { // Metallicity
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) { // H NumberDensity

	index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];

	BaryonField[0][index] = MH * pow(10,((HNumberDensitySlope * (i-GridStartIndex[0])) + log10(TestProblemData.MinimumHNumberDensity))) /
	  TestProblemData.HydrogenFractionByMass / DensityUnits;

	BaryonField[MetalNum][index] = pow(10,((metallicitySlope * (j-GridStartIndex[1])) + log10(TestProblemData.MinimumMetallicity))) *
	  CoolData.SolarMetalFractionByMass * BaryonField[0][index];

      }
    }
  }
 
  /* set velocities */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      BaryonField[vel+dim][i] = 0.0;
 
  /* set density of color fields to user-specified values (if user doesn't specify, 
     the defaults are set in SetDefaultGlobalValues.  Do some minimal amount of error
     checking to try to ensure charge conservation when appropriate */
  for (i = 0; i < size; i++){

    // Set multispecies fields!
    // this attempts to set them such that species conservation is maintained,
    // using the method in CosmologySimulationInitializeGrid.C
    if(TestProblemData.MultiSpecies){

      BaryonField[HINum][i] = TestProblemData.HI_Fraction * 
	TestProblemData.HydrogenFractionByMass * BaryonField[0][i];
 
      BaryonField[HeINum][i] =  TestProblemData.HeI_Fraction *
	BaryonField[0][i] * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeIINum][i] = TestProblemData.HeII_Fraction *
	BaryonField[0][i] * (1.0-TestProblemData.HydrogenFractionByMass);

      BaryonField[HeIIINum][i] =
	(1.0 - TestProblemData.HydrogenFractionByMass) * BaryonField[0][i] -
	BaryonField[HeINum][i] - BaryonField[HeIINum][i];

      if(TestProblemData.MultiSpecies > 1){
	BaryonField[HMNum][i] = TestProblemData.HM_Fraction *
	  TestProblemData.HydrogenFractionByMass * BaryonField[0][i];

	BaryonField[H2INum][i] = 2 * TestProblemData.H2I_Fraction *
	  TestProblemData.HydrogenFractionByMass * BaryonField[0][i];

	BaryonField[H2IINum][i] = 2 * TestProblemData.H2II_Fraction * 
	  TestProblemData.HydrogenFractionByMass * BaryonField[0][i];
      }

      // HII density is calculated by subtracting off the various ionized fractions
      // from the total
      BaryonField[HIINum][i] = TestProblemData.HydrogenFractionByMass * BaryonField[0][i]
	- BaryonField[HINum][i];
      if (MultiSpecies > 1)
	BaryonField[HIINum][i] -= (BaryonField[HMNum][i] + BaryonField[H2IINum][i]
				  + BaryonField[H2INum][i]);

      // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
      // density for convenience) is calculated by summing up all of the ionized species.
      // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
      // calculating mass density, not number density (because the BaryonField values are 4x as
      // heavy for helium for a single electron)
      BaryonField[DeNum][i] = BaryonField[HIINum][i] +
	0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];
      if (MultiSpecies > 1)
	BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] -
	  BaryonField[HMNum][i];

      BaryonField[DeNum][i] = max(BaryonField[DeNum][i], tiny_number);

      // Set deuterium species (assumed to be a negligible fraction of the total, so not
      // counted in the conservation)
      if(TestProblemData.MultiSpecies > 2){
	BaryonField[DINum ][i]  = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HINum][i];
	BaryonField[DIINum][i] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HIINum][i];
	BaryonField[HDINum][i] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * BaryonField[H2INum][i];
      }

    } // if(TestProblemData.MultiSpecies)

  } // for (i = 0; i < size; i++)

  // Now set internal energy from temperature.

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) { // Temperature
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) { // Metallicity
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) { // H NumberDensity

	index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];

	// calculate mu

	if (TestProblemData.MultiSpecies == 0) {
	  mu = DEFAULT_MU;
	}

	else {
	  mu = BaryonField[DeNum][index] + BaryonField[HINum][index] + BaryonField[HIINum][index] + 
	    (BaryonField[HeINum][index] + BaryonField[HeIINum][index] + BaryonField[HeIIINum][index])/4.0;
	  if (MultiSpecies > 1) {
	    mu += BaryonField[HMNum][index] + (BaryonField[H2INum][index] + BaryonField[H2IINum][index])/2.0;
	  }
	  if (MultiSpecies > 2) {
	    mu += (BaryonField[DINum][index] + BaryonField[DIINum][index])/2.0 + (BaryonField[HDINum][index]/3.0);
	  }
	  mu = BaryonField[0][index] / mu;
	}

	BaryonField[1][index] = pow(10,((temperatureSlope * (k-GridStartIndex[2])) + log10(TestProblemData.MinimumTemperature))) /
	  TemperatureUnits / mu / (Gamma-1.0);

	if (DualEnergyFormalism)
	  BaryonField[2][index] = BaryonField[1][index];

      }
    }
  }

  return SUCCESS;
}
