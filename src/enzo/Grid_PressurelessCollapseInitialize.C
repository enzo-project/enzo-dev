/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A PRESSURELESS COLLAPSE)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
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
 
int grid::PressurelessCollapseInitializeGrid(int PressurelessCollapseDirection,
				      float PressurelessCollapseInitialDensity,
				        int PressurelessCollapseNumberOfCells)
{
  /* declarations */
 
  int Divisor, index, size, dim, i;
 
  /* error check */
 
  if (PressurelessCollapseDirection > GridRank-1 ||
      GridDimension[PressurelessCollapseDirection] == 1) {
    ENZO_FAIL("PressurelessCollapseDirection not properly defined.\n");
  }
 
  /* create fields */
 
  NumberOfBaryonFields = 2 + GridRank;
  FieldType[0] = Density;
  FieldType[1] = TotalEnergy;
  FieldType[2] = Velocity1;
  FieldType[3] = Velocity2;
  FieldType[4] = Velocity3;
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Set NumberOfCells if left unspecified. */
 
  if (PressurelessCollapseNumberOfCells == INT_UNDEFINED)
    PressurelessCollapseNumberOfCells =
      GridDimension[PressurelessCollapseDirection] - 2;
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  this->AllocateGrids();
 
  /* set fields */
 
  Divisor = 1;
  for (dim = 0; dim < PressurelessCollapseDirection; dim++)
    Divisor *= GridDimension[dim];
 
  /* set density, total energy and velocity in problem dimension */
 
  for (i = 0; i < size; i++) {
    index = i/Divisor % GridDimension[PressurelessCollapseDirection];
     if (index > GridStartIndex[PressurelessCollapseDirection] &&
	 index < GridEndIndex[PressurelessCollapseDirection])
       *(BaryonField[0] + i) = PressurelessCollapseInitialDensity;
    else
       *(BaryonField[0] + i) = tiny_number;
    *(BaryonField[1] + i) = tiny_number;
    *(BaryonField[2] + i) = 0;
    if (GridRank > 1)
      *(BaryonField[3] + i) = 0;
    if (GridRank > 2)

      *(BaryonField[4] + i) = 0;
  }
 
  return SUCCESS;
}
