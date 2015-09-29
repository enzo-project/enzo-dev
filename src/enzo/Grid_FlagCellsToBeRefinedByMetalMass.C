/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY METAL MASS)
/
/  written by: John Wise
/  date:       March, 2012
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
 
int grid::FlagCellsToBeRefinedByMetalMass(int level)
{

  int i, dim, size, method;
  const float Zsun = CoolData.SolarMetalFractionByMass;
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* error check */

  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  /* Search for mass flagging field */

  int MassFlaggingMethod = INT_UNDEFINED;
  for (method = 0; method < MAX_FLAGGING_METHODS; method++)
    if (CellFlaggingMethod[method] == 2) {
      MassFlaggingMethod = method;
      break;
    }

  if (MassFlaggingMethod == INT_UNDEFINED)
    ENZO_FAIL("Metal mass refinement requires baryon refinement.");

  /* compute size */

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Compute cell volume */
 
  float CellVolume = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    CellVolume *= CellWidth[dim][0];
 
  /* Compute the ModifiedMinimumMass */
 
  float ModifiedMinimumMassForRefinement =
    (MetallicityForRefinement*Zsun) * 
    MinimumMassForRefinement[MassFlaggingMethod] * 
    POW(RefineBy, level*MinimumMassForRefinementLevelExponent[MassFlaggingMethod]);

  /* Metal cooling field numbers. */

  int MetalNum = 0, SNColourNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

  // Double check if there's a metal field when we have metal cooling
  if (MetalCooling && MetalFieldPresent == FALSE)
    ENZO_FAIL("ERROR: No metal field found. "
	      "Restart and turn OFF MetalMass refinement.");

  /* If both metal fields (Pop I/II and III) exist, create a field
     that contains their sum */

  float *MetalPointer;
  float *TotalMetals = NULL;

  if (MetalNum != -1 && SNColourNum != -1) {
    TotalMetals = new float[size];
    for (i = 0; i < size; i++)
      TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
    MetalPointer = TotalMetals;
  } // ENDIF both metal types
  else {
    if (MetalNum != -1)
      MetalPointer = BaryonField[MetalNum];
    else if (SNColourNum != -1)
      MetalPointer = BaryonField[SNColourNum];
  } // ENDELSE both metal types

  /* Flag points */

  for (i = 0; i < size; i++)
    FlaggingField[i] += 
      (CellVolume * MetalPointer[i] > ModifiedMinimumMassForRefinement) ? 1 : 0;
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i] > 0)
      NumberOfFlaggedCells++;
 
  delete [] TotalMetals;
 
  return NumberOfFlaggedCells;
 
}
