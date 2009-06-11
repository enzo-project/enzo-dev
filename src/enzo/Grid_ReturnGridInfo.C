/***********************************************************************
/
/  GRID CLASS (RETURNS INFORMATION ABOUT THIS GRID)
/
/  written by: Greg Bryan
/  date:       May, 1995
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
 
/* Returns some grid info. */
 
int grid::ReturnGridInfo(int *Rank, int Dims[], FLOAT Left[], FLOAT Right[])
{
 
  *Rank = GridRank;
 
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    Dims[dim] = GridDimension[dim];
    Left[dim] = GridLeftEdge[dim];
    Right[dim] = GridRightEdge[dim];
  }
 
  return SUCCESS;
}
 
