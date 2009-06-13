/***********************************************************************
/
/  PROTOSUBGRID CLASS (COPY A NEW SUBGRID TO A NEW PROTOSUBGRID OBJECT)
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
 
int ProtoSubgrid::CopyToNewSubgrid(int GridDim, int GridStart, int GridEnd,
				   ProtoSubgrid *NewSubgrid)
{
 
  /* First copy scalar information. */
 
  NewSubgrid->GridRank = GridRank;
 
  /* Copy all triplet information (although one dim will be changed later). */
 
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    NewSubgrid->GridDimension[dim] = GridDimension[dim];
    NewSubgrid->GridLeftEdge[dim]  = GridLeftEdge[dim];
    NewSubgrid->GridRightEdge[dim] = GridRightEdge[dim];
    NewSubgrid->StartIndex[dim]    = StartIndex[dim];
    NewSubgrid->EndIndex[dim]      = EndIndex[dim];
  }
 
  /* Now reset the GridDim values. */
 
  FLOAT CellWidth = (GridRightEdge[GridDim] - GridLeftEdge[GridDim])/
                    FLOAT(GridDimension[GridDim]);
  NewSubgrid->GridDimension[GridDim] = GridEnd - GridStart + 1;
  NewSubgrid->GridLeftEdge[GridDim]  = GridLeftEdge[GridDim] +
    FLOAT(GridStart - StartIndex[GridDim])*CellWidth;
  NewSubgrid->GridRightEdge[GridDim] = NewSubgrid->GridLeftEdge[GridDim] +
    FLOAT(NewSubgrid->GridDimension[GridDim])*CellWidth;
  NewSubgrid->StartIndex[GridDim]    = GridStart;
  NewSubgrid->EndIndex[GridDim]      = GridEnd;
 
  /* Allocate new FlaggingField. */
 
  NewSubgrid->GridFlaggingField = new int[NewSubgrid->GridDimension[0]*
					  NewSubgrid->GridDimension[1]*
					  NewSubgrid->GridDimension[2]];
 
  /* Now copy the Portion of the FlaggingField. */
 
  FORTRAN_NAME(copy3dint)(GridFlaggingField, NewSubgrid->GridFlaggingField,
                       GridDimension, GridDimension+1, GridDimension+2,
                       NewSubgrid->GridDimension, NewSubgrid->GridDimension+1,
                          NewSubgrid->GridDimension+2,
                       StartIndex, StartIndex+1, StartIndex+2,
                       NewSubgrid->StartIndex, NewSubgrid->StartIndex+1,
                          NewSubgrid->StartIndex+2);
 
  return SUCCESS;
}
