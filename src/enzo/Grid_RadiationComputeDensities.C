/***********************************************************************
/
/  GRID CLASS (COMPUTE VARIOUS DENSITY COMBINATIONS FOR THE RADIATION FIELD)
/
/  written by: Greg Bryan
/  date:       October, 1999
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
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int grid::RadiationComputeDensities(int level)
{
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0 || MultiSpecies < 1)
    return SUCCESS;
 
  /* Initialize. */
 
  int dim, i, j, k, index, size = 1, itemp;
  float ne, nZni, nHI, nHII, nHeI, nHeII, nHeIII, Volume = 1;
  for (dim = 0; dim < GridRank; dim++) {
    Volume *= CellWidth[dim][0];
    size *= GridDimension[dim];
  }
 
  /* Compute temperature field. */
 
  float *temperature = new float[size];
  this->ComputeTemperatureField(temperature);
 
  /* Find Multi-species fields. */
 
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }
 
  /* Get units. */
 
  float DensityUnits = 1, LengthUnits = 1,
    TemperatureUnits = 1, TimeUnits = 1, VelocityUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
 
  float factor = double(LengthUnits)*CellWidth[0][0]*
                 (double(DensityUnits)/double(1.67e-24));
 
  /* Loop over grid and compute densities. */
 
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
 
	/* If this cell is not further refined, the add it in. */
 
	if (BaryonField[NumberOfBaryonFields][index] == 0) {
 
	  /* Compute number densities, but missing factor of rho_units/mh. */
 
	  ne   = BaryonField[DeNum][index];
	  nZni = (     BaryonField[HIINum][index] +
		  0.25*BaryonField[HeIINum][index] +
	               BaryonField[HeIIINum][index]);
	  nHI    =      BaryonField[HINum][index];
	  nHII   =      BaryonField[HIINum][index];
	  nHeI   = 0.25*BaryonField[HeINum][index];
	  nHeII  = 0.25*BaryonField[HeIINum][index];
	  nHeIII = 0.25*BaryonField[HeIIINum][index];
 
	  /* Compute temperature bin for this cell. */
 
	  itemp = int((log10(temperature[index])-
		       RadiationData.TemperatureBinMinimum)/
		      RadiationData.TemperatureBinWidth);
	  itemp = max(min(itemp, RadiationData.NumberOfTemperatureBins-1), 0);
 
	  /* Add density into appropriate bin. The densities are in number
	     densities but missing factor of LengthUnits/mh.  The Volume factor
	     is for average (sum over Volumes = 1). */
 
	  RadiationData.FreeFreeDensity[level][itemp] += ne*nZni*Volume;
	  RadiationData.HIIFreeBoundDensity[level][itemp] += ne*nHII*Volume;
	  RadiationData.HeIIFreeBoundDensity[level][itemp] += ne*nHeII*Volume;
	  RadiationData.HeIIIFreeBoundDensity[level][itemp] +=
	    ne*nHeIII*Volume;
 
	  /* mean HI, etc. densities.  If RadiationFieldType == 11 then
	     use an approximate self-shielding. */
 
	  if (RadiationData.RadiationShield == TRUE) {
	  //	  if (RadiationFieldType == 11) {

	    RadiationData.HIMeanDensity[level] += nHI*Volume*
	      exp(-RadiationData.HIAveragePhotoionizationCrossSection*
		   nHI*factor);
	    RadiationData.HeIMeanDensity[level] += nHeI*Volume;
	      exp(-RadiationData.HeIAveragePhotoionizationCrossSection*
		   nHeI*factor);
	    RadiationData.HeIIMeanDensity[level] += nHeII*Volume;
	      exp(-RadiationData.HeIIAveragePhotoionizationCrossSection*
		   nHeII*factor);
	  } else {
	    RadiationData.HIMeanDensity[level] += nHI*Volume;
	    RadiationData.HeIMeanDensity[level] += nHeI*Volume;
	    RadiationData.HeIIMeanDensity[level] += nHeII*Volume;
	  }
 
	}
      }
    }
 
  /* Clean up */
 
  delete [] temperature;
  delete [] BaryonField[NumberOfBaryonFields];
  BaryonField[NumberOfBaryonFields] = NULL;
 
  return SUCCESS;
}
