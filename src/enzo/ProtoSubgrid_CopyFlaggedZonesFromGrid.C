/***********************************************************************
/
/  PROTOSUBGRID CLASS (SETUP A PROTOSUBGRID FROM A GRID OBJECT)
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
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
 
extern "C" void FORTRAN_NAME(copy3dint)(int *source, int *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
 
int ProtoSubgrid::CopyFlaggedZonesFromGrid(grid *Grid)
{
  /* Initialize. */
 
  int dim, Zero[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Zero[dim] = 0;
 
  /* Error check */
 
  if (Grid->FlaggingField == NULL) {
    ENZO_FAIL("FlaggingField absent in grid!");
  }
 
  /* Set scalars. */
 
  GridRank = Grid->GridRank;
 
  int size = 1;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    GridLeftEdge[dim]  = Grid->GridLeftEdge[dim];
    GridRightEdge[dim] = Grid->GridRightEdge[dim];
    StartIndex[dim]    = Grid->GridStartIndex[dim];
    EndIndex[dim]      = Grid->GridEndIndex[dim];
    GridDimension[dim] = EndIndex[dim] - StartIndex[dim] + 1;
    size *= GridDimension[dim];
  }
 
  /* Allocate and copy GridFlaggingField. */
 
  GridFlaggingField = new int[size];
 
  FORTRAN_NAME(copy3dint)(Grid->FlaggingField, GridFlaggingField,
			   Grid->GridDimension, Grid->GridDimension+1,
                             Grid->GridDimension+2,
			   GridDimension, GridDimension+1, GridDimension+2,
			   Zero, Zero+1, Zero+2,
			   Grid->GridStartIndex, Grid->GridStartIndex+1,
                             Grid->GridStartIndex+2);
 
  return SUCCESS;
}
