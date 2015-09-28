/***********************************************************************
/
/  GRID CLASS (COMPUTE THE TEMPERATURE FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the pressure at the requested time.  The pressure here is
//   just the ideal-gas equation-of-state.
 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
 
/* Set the mean molecular mass. */
 
#define MU_METAL 16.0
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
 
int grid::ComputeTemperatureField(float *temperature,int IncludeCRs)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  int DensNum, result;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
 
  /* If Gadget equilibrium cooling is on, call the appropriate routine,
     then exit - don't use the rest of the routine. */

  if(GadgetEquilibriumCooling){
    if(DualEnergyFormalism)
      result = this->GadgetComputeTemperatureDEF(Time, temperature);
    else
      result = this->GadgetComputeTemperature(Time,temperature);

    if(result == FAIL) {
      ENZO_FAIL("Error in grid->ComputePressure: Gadget.");
    }
    return SUCCESS;
  }

  /* Compute the pressure first. */
 
  this->ComputePressure(Time, temperature,0,IncludeCRs);
 
  /* Compute the size of the fields. */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0)
    ENZO_FAIL("Cannot find density.");

  /* Find metal field */

  int MetalNum = 0;
  bool MetalFieldPresent = false;
  float inv_metal_mol = 1.0 / MU_METAL;

  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) == -1)
    MetalNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1);
   
  if (ProblemType == 60 || ProblemType == 61) { //AK
    for (i = 0; i < size; i++) {
      if (BaryonField[DensNum][i] <= 0.0)
	temperature[i] = 1.0;
      else
	temperature[i] /=  BaryonField[DensNum][i];
    }
    return SUCCESS;
  }
 
  float TemperatureUnits = 1, number_density;
  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1;
 
  /* Find the temperature units. */
 
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  /* For Sedov Explosion compute temperature without floor */

  float mol_weight = Mu, min_temperature = 1.0;   // Mu is now read from parameter file

  // this : would not be needed as we can specify it in parameter file:
  if (ProblemType == 7) {//AK for Sedov explosion test
    mol_weight = 1.0;
    min_temperature = tiny_number;
  }

  if (MultiSpecies == FALSE)
 
    /* If the multi-species flag is not set,
       Compute temperature T = p/d and assume mu = Mu (global data). */
 
    for (i = 0; i < size; i++)
      temperature[i] = max((TemperatureUnits*temperature[i]*mol_weight
		         /max(BaryonField[DensNum][i], tiny_number)),
			 min_temperature);
  else {
 
    /* Find Multi-species fields. */
 
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
			  HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
 
    /* Compute temperature with mu calculated directly. */
 
    for (i = 0; i < size; i++) {
 
      number_density =
	0.25*(BaryonField[HeINum][i]  + BaryonField[HeIINum][i] +
	      BaryonField[HeIIINum][i]                        ) +
              BaryonField[HINum][i]   + BaryonField[HIINum][i]  +
              BaryonField[DeNum][i];
 
      /* Add in H2. */
 
      if (MultiSpecies > 1)
	number_density += BaryonField[HMNum][i]   +
	  0.5*(BaryonField[H2INum][i]  + BaryonField[H2IINum][i]);

      if (MetalFieldPresent)
	number_density += BaryonField[MetalNum][i] * inv_metal_mol;
 
      /* Ignore deuterium. */
 
      temperature[i] *= TemperatureUnits/max(number_density, tiny_number);
      temperature[i] = max(temperature[i], MINIMUM_TEMPERATURE);
    }
  }
 
  return SUCCESS;
}
