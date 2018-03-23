/***********************************************************************
/
/  GRID CLASS (APPROXIMATE DETECTION METHOD OF IONIZATION FRONTS)
/
/  written by: John Wise
/  date:       August, 2008
/  modified1:
/
/  PURPOSE: 
/
/  RETURNS: TRUE or FALSE
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

int grid::DetectIonizationFrontApprox(float TemperatureUnits)
{

  if (MyProcessorNumber != ProcessorNumber)
    return FALSE;

  if (!this->RadiationPresent())
    return FALSE;

  int i, j, k, n, index, ret, off[3];
  //  float max_kph = 0, min_kph = 1e20;
  float maxDT, dT;
  float TOLERANCE = 5e3;  // if abs[dT(i+1) - dT(i)] > TOLERANCE, then we return TRUE

  int eNum, kphHINum, gammaNum, kphHeINum, kphHeIINum, 
    kdissH2INum;

  if (DualEnergyFormalism)
    eNum = FindField(InternalEnergy, FieldType, NumberOfBaryonFields);
  else
    eNum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields);

  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum);

  maxDT = -1e20;
  TOLERANCE /= TemperatureUnits;

  // For convenience, offsets in 1D array corresponding to a shift in the grid,
  // i.e. i+1, j+1, k+1

  off[0] = 1;
  off[1] = GridDimension[0];
  off[2] = GridDimension[0] * GridDimension[1];

  // Search to the second-to-rightmost cell since we don't want to
  // include ghost zones because they may not have been updated yet,
  // especially when using the coupled RT/energy solver.
  for (k = GridStartIndex[2]; k < GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
      for (i = GridStartIndex[0]; i < GridEndIndex[0]; i++, index++) {
	if (BaryonField[kphHINum][index] > 0)
	  for (n = 0; n < 3; n++) {
	    dT = fabs(BaryonField[eNum][index+off[n]] - BaryonField[eNum][index]);
	    maxDT = max(maxDT, dT);
	  } // ENDFOR n
      } // ENDFOR i
    } // ENDFOR j

  ret = (maxDT > TOLERANCE) ? TRUE : FALSE;
  //printf("DetectIFront[P%"ISYM"]: max (ret) = %"GSYM" (%"ISYM")\n", MyProcessorNumber, maxDT, ret);

  return ret;

}
