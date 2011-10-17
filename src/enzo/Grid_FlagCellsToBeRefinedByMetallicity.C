/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINE BY Metallicity)
/
/  written by: Brian O'Shea
/  date:       August 2006
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
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
 
int FindField(int field, int farray[], int numfields);
 
int grid::FlagCellsToBeRefinedByMetallicity(int level)
{
  /* declarations */
 
  int i, j, k, index, dim, NumberOfFlaggedCells = 0;
 
  /* Return if this grid is not on this processor. */
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* error check */
   if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 

  /* if we're already at the max level that we want to refine to,
     return!  Otherwise we want to refine here. */
  if(MetallicityRefinementMinLevel <= level) return NumberOfFlaggedCells;


  /* Find metallicity field.  If no metallicity field exists, quit and
     yell at user. */
  int MetalNum, SNColourNum;
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);

  if (MetalNum < 0 && SNColourNum < 0) {
    fprintf(stderr,"FlagCellsToBeRefinedByMetallicity: no metallicity field!\n");
    return -1;
  }

  /* Find fields: density, total energy, velocity1-3. */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return -1;
  }


  /* compute size */ 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  if (MetalNum > 0)
    for(i=0; i<size; i++)
      if( (BaryonField[MetalNum][i]/BaryonField[DensNum][i])/0.022 
	  >= MetallicityRefinementMinMetallicity &&
	  (BaryonField[DensNum][i] > MetallicityRefinementMinDensity))
	FlaggingField[i]++;
  if (SNColourNum > 0)
    for(i=0; i<size; i++)
      if( (BaryonField[SNColourNum][i]/BaryonField[DensNum][i])/0.022 
	  >= MetallicityRefinementMinMetallicity &&
	  (BaryonField[DensNum][i] > MetallicityRefinementMinDensity))
	FlaggingField[i]++;


  /* Count number of flagged Cells. */ 
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }

  return NumberOfFlaggedCells;
 
}
