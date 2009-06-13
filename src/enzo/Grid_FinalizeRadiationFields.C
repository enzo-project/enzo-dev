#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (FINALIZE THE PHOTO-RATES BY DIVIDING BY NUMBER OF PARTICLES)
/
/  written by: John Wise
/  date:       June, 2009
/  modified1:
/
/  PURPOSE: 
/
************************************************************************/


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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);

int grid::FinalizeRadiationFields(void)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

#ifdef TRANSFER

  int i, j, k, index, dim;
  float CellVolume = 1;

  for (dim = 0; dim < GridRank; dim++)
    CellVolume *= CellWidth[dim][0];

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifySpeciesFields.\n");
    ENZO_FAIL("");
  }

  /* Find radiative transfer fields. */

  int kphHINum, gammaHINum, kphHeINum, gammaHeINum, kphHeIINum, gammaHeIINum,
    kdissH2INum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				      gammaHeINum, kphHeIINum, gammaHeIINum, 
				      kdissH2INum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifyRadiativeTransferFields.\n");
    ENZO_FAIL("");
  }

  /* Get units. */

  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    MassUnits, DensityUnits; 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, PhotonTime) == FAIL) {
    fprintf(stdout, "Error in GetUnits.\n");
    ENZO_FAIL("");
  }

  float DensityConversion = DensityUnits / 1.673e-24;
  float factor = DensityConversion * CellVolume;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	BaryonField[kphHINum][index] /= (factor * BaryonField[HINum][index]);
	BaryonField[gammaHINum][index] /= (factor * BaryonField[HINum][index]);
	BaryonField[kphHeINum][index] /= (factor * BaryonField[HeINum][index]);
	BaryonField[gammaHeINum][index] /= (factor * BaryonField[HeINum][index]);
	//BaryonField[kphHeIINum][index] /= (factor * BaryonField[HeIINum][index]);
	BaryonField[gammaHeIINum][index] /= (factor * BaryonField[HeIINum][index]);
      } // ENDFOR i
    } // ENDFOR j
  

#endif /* TRANSFER */  
  
  return SUCCESS;
}
