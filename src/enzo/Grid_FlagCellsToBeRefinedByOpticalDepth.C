/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY OPTICAL DEPTH)
/
/  written by: John H. Wise
/  date:       December, 2005
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
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
#include "Grid.h"

/* function prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

#define MAX_TAU 1

int grid::FlagCellsToBeRefinedByOpticalDepth()
{
  /* declarations */

  int i, j, k, dim;

  /* error check */

  if (FlaggingField == NULL) {
    ENZO_FAIL("Flagging Field is undefined.\n");
  }

  /* compute size */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }

  /* Find radiative transfer fields. */

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  /* Get density units. */

  float DensityUnits, LengthUnits, VelocityUnits, TimeUnits,
        TemperatureUnits;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Calculate conversion factor to optical depth */

  float sigmaHI = 6.0e-18 * LengthUnits;
  float ConvertToProperNumberDensity =
    (float) (double(DensityUnits)/double(1.67e-24));
  float OpticalDepthConversion = 
    ConvertToProperNumberDensity * CellWidth[0][0] * sigmaHI;

  float tau0, tau1, tau2;
  float inv_dt_sec = 1.0 / (dtFixed * TimeUnits);

  /* Loop over grid. */

  int NumberOfFlaggedCells_TAU = 0;
  int index, offset = 1;
  float avgTau = 0, avg_kph = 0;
  float minTau = 1e20, maxTau = -1e20;
  float minkph = 1e20, maxkph = -1e20;

  //  for (dim = 0; dim < GridRank; dim++) {

    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = (j + k*GridDimension[1])*GridDimension[0];
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	  tau0 = OpticalDepthConversion * BaryonField[HINum][index+i-offset];
	  tau1 = OpticalDepthConversion * BaryonField[HINum][index+i];
	  tau2 = OpticalDepthConversion * BaryonField[HINum][index+i+offset];
	  if (BaryonField[kphHINum][index+i] > tiny_number && // inv_dt_sec &&
	      //	      max(tau0, MAX(tau1,tau2)) > MAX_TAU) {
	      tau1 > MAX_TAU) {

	    FlaggingField[index+i]++;
	    NumberOfFlaggedCells_TAU++;

	    avgTau  += tau1;//tau0 + tau1 + tau2;
	    avg_kph += BaryonField[kphHINum][index+i];
	    minTau = min(tau1, minTau);
	    maxTau = max(tau1, maxTau);
	    minkph = min(BaryonField[kphHINum][index+i], minkph);
	    maxkph = max(BaryonField[kphHINum][index+i], maxkph);

//	    if (dim == 0)
//	      printf("FlagTau: kph = %"FSYM", idx = %"ISYM"/%"ISYM" %"ISYM"/%"ISYM" %"ISYM"/%"ISYM" (%"ISYM")\n",
//		     BaryonField[kphHINum][index+i], i, GridDimension[0], 
//		     j, GridDimension[1], k, GridDimension[2], index+i);
	    
	    /* Flag boundary cells, but don't flag ghost cells */

//	    if (k > GridStartIndex[2] && j > GridStartIndex[1] &&
//		i > GridStartIndex[0])
//	      FlaggingField[index-offset]++;
//
//	    if (k < GridEndIndex[2] && j < GridEndIndex[1] &&
//		i < GridEndIndex[0])
//	      FlaggingField[index+offset]++;

	  }  /* ENDIF flag */
	} /* ENDFOR i */
      } /* ENDFOR j */

    //    offset *= GridDimension[dim];

    //  } /* ENDFOR dimesnions */

  /* Count number of flagged Cells. */

  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i] > 0)
      NumberOfFlaggedCells++;

  if (NumberOfFlaggedCells_TAU) {

    fprintf(stdout, "FlagCellsOpticalDepth: %"ISYM"\n", NumberOfFlaggedCells_TAU);
    fprintf(stdout, "FlagCellsOpticalDepth: avg(kph) = %"GSYM", avg(tau) = %"GSYM"\n",
	    avg_kph/(NumberOfFlaggedCells_TAU), avgTau/NumberOfFlaggedCells_TAU);
    fprintf(stdout, "FlagCellsOpticalDepth: MIN(kph) = %"GSYM", MAX(kph) = %"GSYM"\n",
	    minkph, maxkph);
    fprintf(stdout, "FlagCellsOpticalDepth: MIN(tau) = %"GSYM", MAX(tau) = %"GSYM"\n",
	    minTau, maxTau);
  }

  return NumberOfFlaggedCells;

}
