/***********************************************************************
/
/  GRID CLASS (COMPUTES GRID DERIVED QUANTITES SUCH AS BOUNDARY FLUXES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
//  Compute derived quantites
//
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
void grid::PrepareGridDerivedQuantities()
{
 
  /*  Baryons only: set up field quantities and allocate fields
      (we assume here the grid is uniform in each dimension) */
 
  FLOAT AllCellWidth, GridLeftIncludingBoundary;
  int dim, j;
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* create cell position descriptors */
 
    CellLeftEdge[dim] = new FLOAT[GridDimension[dim]];
    CellWidth[dim]    = new FLOAT[GridDimension[dim]];
 
    /* assuming uniform grid, compute the cell width and the left edge
       including the ghostzones (for this dimension) */
 
    AllCellWidth    = (GridRightEdge[dim] - GridLeftEdge[dim])/
                      FLOAT(GridEndIndex[dim] - GridStartIndex[dim] + 1);
    GridLeftIncludingBoundary = GridLeftEdge[dim] -
                                AllCellWidth*FLOAT(GridStartIndex[dim]);
 
    /* based on these quantities, set the cell position descriptors */
 
    for (j = 0; j < GridDimension[dim]; j++) {
      CellLeftEdge[dim][j] = GridLeftIncludingBoundary + AllCellWidth*FLOAT(j);
      CellWidth[dim][j]    = AllCellWidth;
    }
 
  } // end: loop over dims
 
}
 
