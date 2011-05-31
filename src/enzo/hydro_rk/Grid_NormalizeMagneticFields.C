/***********************************************************************
/
/  GRID CLASS (PREPARE NORMALIZATION FOR MAGNETIC FIELDS)
/            
/
/  written by: Tom Abel
/  date:       June, 2011
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
 
int grid::NormalizeMagneticFields(Eflt fac)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

  //  if (fac == 0.) 
  //    return SUCCESS; // wouldn't call this a normalization

  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
  fprintf(stderr, "Grid_NormalizeMagneticFields: \n");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    fprintf(stderr, "GPRFN: Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }

  int i, j, k, index, dim;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
	BaryonField[TENum][index] -= 0.5*(pow(BaryonField[B1Num][index], 2) + 
					  pow(BaryonField[B2Num][index], 2) + 
					  pow(BaryonField[B3Num][index], 2) )
	                                  /BaryonField[DensNum][index];

	for (dim = 0; dim < GridRank; dim++) {
	  int bf = B1Num + dim;
	  BaryonField[bf][index] *= fac; 
	}
	BaryonField[TENum][index] += 0.5*(pow(BaryonField[B1Num][index], 2) + 
					  pow(BaryonField[B2Num][index], 2) + 
					  pow(BaryonField[B3Num][index], 2) )
	                                  /BaryonField[DensNum][index];

      }

  return SUCCESS;
}
