/***********************************************************************
/
/  GRID CLASS (ZERO THE SOLUTION VECTOR UNDER THE SUBGRID SPECIFIED)
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:
/
/  PURPOSE:
/
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
 
 
int grid::ZeroSolutionUnderSubgrid(grid *Subgrid, int FieldsToZero,
				   float Value, int AllProcessors,
				   int IncludeGhostZones)
{
 
  /* Return if this grid is not on this processor. */
  //  printf("Valeu:%f\n", Value); 
  if (AllProcessors == FALSE)
    if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;
 
  /* declarations */
 
  int i, j, k, dim, field, index;
  int SubgridStart[MAX_DIMENSION], SubgridEnd[MAX_DIMENSION];
 
  /* If FieldsToZero is ZERO_UNDER_SUBGRID_FIELD, and Subgrid is NULL,
     initialize the under subgrid fields (which is the next unused
     BaryonField even though this pretty much a dumb thing to do). */
 
  if (FieldsToZero == ZERO_UNDER_SUBGRID_FIELD && Subgrid == NULL) {
 
//    printf("ZeroSUS - ZERO_UNDER_SUBGRID_FIELD && Subgrid == NULL\n");
 
    int size = 1;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      size *= GridDimension[dim];
 
 
    delete [] BaryonField[NumberOfBaryonFields];
 
    if ( (Unigrid == 1) && ((ProblemType == 30) || (ProblemType == 60) || (ProblemType >= 400)) ) //AK
    {
       printf("ZeroSUS - ZERO_UNDER_SUBGRID_FIELD && Subgrid == NULL\n");
       printf("ZeroSUS - Unigrid: %"ISYM"\n", Unigrid);
       printf("ZeroSUS - ProblemType: %"ISYM"\n", ProblemType);
       printf("ZeroSUS - NumberOfBaryonFields: %"ISYM"\n", NumberOfBaryonFields);
       printf("ZeroSUS - Zero field size: %"ISYM"\n", (int) (size*sizeof(float)));
       printf("ZeroSUS - CAUTION - allocating bogus BaryonField size = 1\n");
       size = 1;
    }
 
//    else
//    {
//       printf("ZeroSUS - NumberOfBaryonFields: %"ISYM"\n", NumberOfBaryonFields);
//       printf("ZeroSUS - Zero field size: %"ISYM"\n", (int) (size*sizeof(float)));
//    }
 
    BaryonField[NumberOfBaryonFields] = new float[size];
    for (i = 0; i < size; i++)
      BaryonField[NumberOfBaryonFields][i] = 0.0;
    return SUCCESS;
  }
 
 
//  printf("ZeroSUS - Subgrid not NULL\n");
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    SubgridStart[dim] = 0;
    SubgridEnd[dim] = 0;
  }
 
  /* Compute start and stop indices of the active region of the subgrid
     within this grid (and check to make sure subgrid is within this grid). */

  FLOAT Left, Right, SubLeft, SubRight;
  int NumberOfGhostZones;
 
  for (dim = 0; dim < GridRank; dim++) {

    if (IncludeGhostZones == FALSE) {

      if (Subgrid->GridRightEdge[dim] <= GridLeftEdge[dim] ||
	  Subgrid->GridLeftEdge[dim]  >= GridRightEdge[dim])
	return SUCCESS;
      
      SubgridStart[dim] = 
	nint((Subgrid->GridLeftEdge[dim] - GridLeftEdge[dim]) / 
	     CellWidth[dim][0]) + GridStartIndex[dim];
      SubgridEnd[dim] = 
	nint((Subgrid->GridRightEdge[dim] - GridLeftEdge[dim]) / 
	     CellWidth[dim][0]) + GridStartIndex[dim] - 1;

      SubgridStart[dim] = max(SubgridStart[dim], GridStartIndex[dim]);
      SubgridEnd[dim]   = min(SubgridEnd[dim], GridEndIndex[dim]);

    } 
    
    // Include ghost zones
    else {

      NumberOfGhostZones = GridStartIndex[dim];
      Left = GridLeftEdge[dim] - CellWidth[dim][0] * NumberOfGhostZones;
      Right = GridRightEdge[dim] + CellWidth[dim][GridDimension[dim]-1]
	* NumberOfGhostZones;

      // We only want to mark where the other grids (no ghost zones)
      // overlap with this grid (including the ghost zones).
      
      if (Subgrid->GridRightEdge[dim] <= Left || 
	  Subgrid->GridLeftEdge[dim] >= Right)
	return SUCCESS;

      SubgridStart[dim] = nint((Subgrid->GridLeftEdge[dim] - Left) / 
			       CellWidth[dim][0]);
      SubgridEnd[dim] = nint((Subgrid->GridRightEdge[dim] - Left) / 
			     CellWidth[dim][0]) - 1;

      SubgridStart[dim] = max(SubgridStart[dim], 0);
      SubgridEnd[dim]   = min(SubgridEnd[dim], GridDimension[dim]-1);
    }
 
//    printf("  ZeroSUS: %"ISYM", %"ISYM", %"ISYM"\n", dim, SubgridStart[dim], SubgridEnd[dim]);
  }
 
  /* Now that there is overlap, take the appropriate action. */
 
  if (FieldsToZero == ZERO_ALL_FIELDS){
 
//    printf("  ZeroSUS == ZERO_ALL_FIELDS\n");
 
    if (RandomForcing) {
      printf("Forcing Fields will be modified (ZeroSUS == ZERO_ALL_FIELDS)\n");
      WARNING_MESSAGE;
    }
 
    /* Zero all fields in this region. */
 
    for (field = 0; field < NumberOfBaryonFields; field++)
      for (k = SubgridStart[2]; k <= SubgridEnd[2]; k++)
	for (j = SubgridStart[1]; j <= SubgridEnd[1]; j++) {
	  index = (k*GridDimension[1] + j)*GridDimension[0];
	  for (i = SubgridStart[0]; i <= SubgridEnd[0]; i++)
	    BaryonField[field][index + i] = 0;
	}
  }
  else if (FieldsToZero == ZERO_UNDER_SUBGRID_FIELD) {
 
//    printf("  ZeroSUS = ZERO_UNDER_SUBGRID_FIELD\n");
 
    if (BaryonField[NumberOfBaryonFields] == NULL)
      ENZO_FAIL("UNDER_SUBGRID_FIELD not allocated.");
 
    /* Set points under this subgrid to Value (typically 1). */
 
//    printf("    ZeroSUS Value = %10.4e\n", Value);
 
    if ( (Unigrid == 1) && ((ProblemType == 30) || (ProblemType == 60) || (ProblemType >= 400)) ) //AK
    {
      printf("    Ignoring BaryonField Assignment!\n");
    }
    else
    {
    for (k = SubgridStart[2]; k <= SubgridEnd[2]; k++)
      for (j = SubgridStart[1]; j <= SubgridEnd[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0];
	for (i = SubgridStart[0]; i <= SubgridEnd[0]; i++)
	  BaryonField[NumberOfBaryonFields][index + i] = Value;
      }
    }
 
  }
 
  else {
    ENZO_VFAIL("FieldsToZero = %"ISYM" not recognized.\n", FieldsToZero)
  }
 
  return SUCCESS;
 
}
