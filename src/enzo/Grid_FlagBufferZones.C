/***********************************************************************
/
/  GRID CLASS (FLAG BUFFER CELLS AND REMOVE BOUNDARY (+1) CELLS)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
// This routine flags all cells which are adjacent (diagonally counts) to
//  an already flagged Cell.  It also removes all flagged cells which are
//  in the boundary.
 
#include <stdio.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::FlagBufferZones()
{
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  int i, j, k, dim;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    exit(EXIT_FAILURE);
  }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Flag buffer points (i.e. all cells that surround any existing cell.
     At the same time, copy the TempBuffer into the FlaggingField. */
 
  int Index, Offset = 1;
 
  /* allocate a temporary Flagging field buffer */
 
  int *TempBuffer = new int[size];
  for (i = 0; i < size; i++)
    TempBuffer[i] = 0;
 
  /* loop over each active face */
 
  for (int ibuffer = 0; ibuffer<NumberOfBufferZones; ibuffer++){
    Offset = 1;
    //Do this for all cells except last ghost zone.
    int kstart = 1;
    int kend = GridDimension[2]-1;
    int jstart = 1;
    int jend = GridDimension[1]-1;
    if(GridRank < 3){
      kstart=0;
      kend=1;
    }
    if(GridRank < 2){
      jstart=0;
      jend=1;
    }

    for (dim = 0; dim < GridRank; dim++){
      if (GridDimension[dim] > 1) {

        /* flag buffer zones to the left & right for this dimension */

        for(k=kstart;k<kend;k++){
          for(j=jstart; j<jend;j++){
            Index = (k*GridDimension[1] + j)*GridDimension[0] + 1;
            for(i=1; i<GridDimension[0]-1;i++, Index++){
              TempBuffer[Index] += FlaggingField[Index         ] + 
                FlaggingField[Index - Offset] +
                FlaggingField[Index + Offset];
            }
          }
        }

        /* Update offset. */

        Offset *= GridDimension[dim];

        /* Copy TempBuffer back to FlaggingField. */

        for (i = 0; i < size; i++)
          FlaggingField[i] = TempBuffer[i];
      }
    } 
  }
  

  /* Remove points in ghost zones. */
 
  if (GridDimension[0] > 1)
    for (k = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++) {
        Index = (k*GridDimension[1] + j)*GridDimension[0];
        for (i = 0; i < GridStartIndex[0]; i++){
          FlaggingField[Index+i] = 0;
        }
        for (i = GridEndIndex[0]+1; i < GridDimension[0]; i++){
          FlaggingField[Index+i] = 0;
        }
      }
 
  if (GridDimension[1] > 1)
    for (k = 0; k < GridDimension[2]; k++) {
 
      for (j = 0; j < GridStartIndex[1]; j++) {
        Index = (k*GridDimension[1] + j)*GridDimension[0];
        for (i = 0; i < GridDimension[0]; i++, Index++){
          FlaggingField[Index] = 0;
        }
      }

      for (j = GridEndIndex[1]+1; j < GridDimension[1]; j++) {
        Index = (k*GridDimension[1] + j)*GridDimension[0];
        for (i = 0; i < GridDimension[0]; i++, Index++){
          FlaggingField[Index] = 0;
        }
      }
    }
 
  if (GridDimension[2] > 1)
    for (j = 0; j < GridDimension[1]; j++) {
 
      for (k = 0; k < GridStartIndex[2]; k++) {
        Index = (k*GridDimension[1] + j)*GridDimension[0];
        for (i = 0; i < GridDimension[0]; i++, Index++){
          FlaggingField[Index] = 0;
        }
      }
 
      for (k = GridEndIndex[2]+1; k < GridDimension[2]; k++) {
        Index = (k*GridDimension[1] + j)*GridDimension[0];
        for (i = 0; i < GridDimension[0]; i++, Index++){
          FlaggingField[Index] = 0;
        }
      }
    }
 
  /* Set flag to zero if outside allowed refine region. */
 
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
 
      Index = (k*GridDimension[1] + j)*GridDimension[0] +
	GridStartIndex[0];
 
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, Index++) {
 
	if (CellLeftEdge[0][i] + 0.5*CellWidth[0][i] <
	    RefineRegionLeftEdge[0]                  ||
	    CellLeftEdge[0][i] + 0.5*CellWidth[0][i] >
	    RefineRegionRightEdge[0])
	  FlaggingField[Index] = 0;
 
	if (GridRank > 1)
	  if (CellLeftEdge[1][j] + 0.5*CellWidth[1][j] <
	      RefineRegionLeftEdge[1]                  ||
	      CellLeftEdge[1][j] + 0.5*CellWidth[1][j] >
	      RefineRegionRightEdge[1])
	    FlaggingField[Index] = 0;
 
	if (GridRank > 2)
	  if (CellLeftEdge[2][k] + 0.5*CellWidth[2][k] <
	      RefineRegionLeftEdge[2]                  ||
	      CellLeftEdge[2][k] + 0.5*CellWidth[2][k] >
	      RefineRegionRightEdge[2])
	    FlaggingField[Index] = 0;
	
	FlaggingField[Index] = min(FlaggingField[Index], 1);
      }
    }
 
  /* clean up */
 
  delete [] TempBuffer;
 
  /* If debuging, count up the number of flagged cells & report. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i] > 0)
      NumberOfFlaggedCells++;
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridEndIndex[dim] - GridStartIndex[dim] + 1;
 
  if (debug1)
    printf("FlagBufferZones: NumberOfFlaggedCells = %"ISYM" (%.1"FSYM"%%)\n",
	   NumberOfFlaggedCells, float(NumberOfFlaggedCells)*100.0/
	   float(size));
 
  return NumberOfFlaggedCells;
 
}
