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
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
 
/* Set the mean molecular mass. */
 
#define DEFAULT_MU 0.6
 
/* This is minimum returned temperature. (K) */
 
#define MINIMUM_TEMPERATURE 1.0
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
 
 
int grid::ComputeTemperatureField(float *temperature)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  int DensNum, result;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
 
  /* Compute the pressure first. */
 
  if (DualEnergyFormalism)
    result = this->ComputePressureDualEnergyFormalism(Time, temperature);
  else
    result = this->ComputePressure(Time, temperature);
 
  if (result == FAIL) {
    fprintf(stderr, "Error in grid->ComputePressure.\n");
    return FAIL;
  }
 
  /* Compute the size of the fields. */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    fprintf(stderr, "Cannot find density.\n");
    return FAIL;
  }
 
 
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
  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1, MassUnits=1;
 
  /* Find the temperature units. */
 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }
 
  if (MultiSpecies == FALSE)
 
    /* If the multi-species flag is not set,
       Compute temperature T = p/d and assume mu = DEFAULT_MU. */
 
    for (i = 0; i < size; i++)
      temperature[i] = max((TemperatureUnits*temperature[i]*DEFAULT_MU
		         /max(BaryonField[DensNum][i], tiny_number)),
			 MINIMUM_TEMPERATURE);
  else {
 
    /* Find Multi-species fields. */
 
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      return FAIL;
    }
 
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
 
      /* Ignore deuterium. */
 
      temperature[i] *= TemperatureUnits/max(number_density, tiny_number);
      temperature[i] = max(temperature[i], MINIMUM_TEMPERATURE);
    }
  }
 
  return SUCCESS;
}
