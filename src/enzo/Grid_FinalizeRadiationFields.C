#define DEBUG 0
#define MYPROC MyProcessorNumber == ProcessorNumber

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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

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
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }

  int HeHIINum, DMNum   , HDIINum
    , CINum   , CIINum  , CONum     , CO2Num   , OINum   , OHNum
    , H2ONum  , O2Num   , SiINum    , SiOINum  , SiO2INum
    , CHNum   , CH2Num  , COIINum   , OIINum   , OHIINum , H2OIINum, H3OIINum, O2IINum
    , MgNum   , AlNum   , SNum      , FeNum
    , SiMNum  , FeMNum  , Mg2SiO4Num, MgSiO3Num, Fe3O4Num
    , ACNum   , SiO2DNum, MgONum    , FeSNum   , Al2O3Num
    , DustNum ;
  if (IdentifySpeciesFieldsMD( HeHIINum, DMNum   , HDIINum
                             , CINum   , CIINum  , CONum     , CO2Num   , OINum   , OHNum
                             , H2ONum  , O2Num   , SiINum    , SiOINum  , SiO2INum
                             , CHNum   , CH2Num  , COIINum   , OIINum   , OHIINum , H2OIINum,  H3OIINum,  O2IINum
                             , MgNum   , AlNum   , SNum      , FeNum
                             , SiMNum  , FeMNum  , Mg2SiO4Num, MgSiO3Num, Fe3O4Num
                             , ACNum   , SiO2DNum, MgONum    , FeSNum   , Al2O3Num
                             , DustNum ) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFieldsMD.\n");
  }

  /* Find radiative transfer fields. */

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  int kdissHDINum, kphCINum, kphOINum, kdissCONum, kdissOHNum, kdissH2ONum;
  IdentifyRadiativeTransferFieldsMD(kdissHDINum, kphCINum, kphOINum, kdissCONum, kdissOHNum, kdissH2ONum);

  /* Get units. */

  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits; 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, PhotonTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  float DensityConversion = DensityUnits / mh;
  float factor = DensityConversion * CellVolume;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	BaryonField[kphHINum][index] /= factor * BaryonField[HINum][index];
	BaryonField[gammaNum][index] /= factor * BaryonField[HINum][index]; //divide by N_HI = n_HI*(dx)^3
      } // ENDFOR i
    } // ENDFOR j

  if (RadiativeTransferHydrogenOnly == FALSE)
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  BaryonField[kphHeINum][index] /= 
	    0.25 * factor * BaryonField[HeINum][index];
	  BaryonField[kphHeIINum][index] /= 
	    0.25 * factor * BaryonField[HeIINum][index];
	} // ENDFOR i
      } // ENDFOR j
  
   if (MultiSpecies > 1 && !RadiativeTransferFLD)
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  if(!RadiativeTransferUseH2Shielding) {
	    BaryonField[kdissH2INum][index] /= 
	     0.5 * factor * BaryonField[H2INum][index];
	  }
	  BaryonField[kphHMNum][index] /= 
	    1.0 * factor * BaryonField[HMNum][index];
	  BaryonField[kdissH2IINum][index] /= 
	    0.5 * factor * BaryonField[H2IINum][index];
	} // ENDFOR i
      } // ENDFOR j

   if (MultiSpecies > 2)
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  if(!RadiativeTransferUseH2Shielding) {
	    BaryonField[kdissHDINum][index] /= 
	      factor/3.0 * BaryonField[HDINum][index];
          }
	} // ENDFOR i
      } // ENDFOR j

   if (MetalChemistry)
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  BaryonField[kphCINum][index] /= 
	    factor/12.0 * BaryonField[CINum][index];
	  BaryonField[kphOINum][index] /= 
	    factor/16.0 * BaryonField[OINum][index];
	  if(!RadiativeTransferUseH2Shielding) {
	    BaryonField[kdissCONum][index] /= 
	      factor/28.0 * BaryonField[CONum][index];
	    BaryonField[kdissOHNum][index] /= 
	      factor/17.0 * BaryonField[OHNum][index];
	    BaryonField[kdissH2ONum][index] /= 
	      factor/18.0 * BaryonField[H2ONum][index];
          }
	} // ENDFOR i
      } // ENDFOR j

   if(this->IndexOfMaximumkph >= 0)
     this->MaximumkphIfront /= (factor * BaryonField[HINum][IndexOfMaximumkph]);

#endif /* TRANSFER */  
  
  return SUCCESS;
}
