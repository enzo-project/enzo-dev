#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (COMPUTE A CELL CENTERED FIELD FROM VERTEX-CENTERED)
/
/  written by: John Wise
/  date:       February, 2008
/  modified1:
/
/  PURPOSE: Takes vertex-centered field in InterpolatedField[] and 
/           computes the cell-centered counterpart are puts it in 
/           BaryonField[].
/
************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"

int grid::ComputeCellCenteredField(int Num)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Error Check */

  if (InterpolatedField[Num] == NULL) {
    ENZO_VFAIL("Interpolated field #%"ISYM" does not exist.\n", Num)
  }

  /* declarations */

  int i,j,k, index, vc_index, dim, size = 1;
  int vi[8];
  float weight;

  /* Compute the size of the new grid adds one cell to the active region */

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (i = 0; i < size; i++)
    BaryonField[Num][i] = 0.0f;

  /* For radiation, take into account child grids and ghost zones,
     both which have no radiation ever. */

#ifdef TRANSFER
  if (SubgridMarker != NULL && 
      InterpolatedField[NumberOfBaryonFields] == NULL) {
    BaryonField[NumberOfBaryonFields] = new float[size];
    for (i = 0; i < size; i++)
      BaryonField[NumberOfBaryonFields][i] = 0.0f;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
	  BaryonField[NumberOfBaryonFields][index] =
	    (this == SubgridMarker[index]) ? 1 : 0;
      } // ENDFOR j
    this->ComputeVertexCenteredField(NumberOfBaryonFields);
  } // ENDIF SubGridMarker
#endif

  if (GridRank == 1) 
    for (i = 0; i <= GridEndIndex[0]-GridStartIndex[0]; i++) 
      BaryonField[Num][GRIDINDEX(i,0,0)] = 0.5f * 
	(InterpolatedField[Num][i] + InterpolatedField[Num][i+1]);

  if (GridRank == 2) 
    for (j = 0; j <= GridEndIndex[1]-GridStartIndex[1]; j++) 
      for (i = 0; i <= GridEndIndex[0]-GridStartIndex[0]; i++) 
	BaryonField[Num][GRIDINDEX(i,j,0)] = 0.25f * 
	  (InterpolatedField[Num][VCGRIDINDEX(i,j,0)] + 
	   InterpolatedField[Num][VCGRIDINDEX(i+1,j,0)] +
	   InterpolatedField[Num][VCGRIDINDEX(i,j+1,0)] + 
	   InterpolatedField[Num][VCGRIDINDEX(i+1,j+1,0)]);

  if (GridRank == 3) 
#ifdef TRANSFER
    if (SubgridMarker == NULL) {
#endif
      for (k = 0; k <= GridEndIndex[2]-GridStartIndex[2]; k++) 
	for (j = 0; j <= GridEndIndex[1]-GridStartIndex[1]; j++) {
	  index = GRIDINDEX(0,j,k);
	  for (i = 0; i <= GridEndIndex[0]-GridStartIndex[0]; i++, index++)
	    BaryonField[Num][index] = 0.125f * 
	      (InterpolatedField[Num][VCGRIDINDEX(i,j,k)] + 
	       InterpolatedField[Num][VCGRIDINDEX(i+1,j,k)] +
	       InterpolatedField[Num][VCGRIDINDEX(i,j+1,k)] + 
	       InterpolatedField[Num][VCGRIDINDEX(i+1,j+1,k)] +
	       InterpolatedField[Num][VCGRIDINDEX(i,j,k+1)] + 
	       InterpolatedField[Num][VCGRIDINDEX(i+1,j,k+1)] +
	       InterpolatedField[Num][VCGRIDINDEX(i,j+1,k+1)] + 
	       InterpolatedField[Num][VCGRIDINDEX(i+1,j+1,k+1)]);
	} // ENDFOR j
#ifdef TRANSFER
    } // ENDIF no subgrid marker
    else {

      for (k = 0; k <= GridEndIndex[2]-GridStartIndex[2]; k++) 
	for (j = 0; j <= GridEndIndex[1]-GridStartIndex[1]; j++) {
	  index = GRIDINDEX(0,j,k);
	  vc_index = VCGRIDINDEX(0,j,k);
	  this->NeighborVertexIndices(vc_index, vi);
	  for (i = 0; i <= GridEndIndex[0]-GridStartIndex[0]; i++, index++) {
	    weight = InterpolatedField[NumberOfBaryonFields][vi[0]+i] +
	      InterpolatedField[NumberOfBaryonFields][vi[1]+i] +
	      InterpolatedField[NumberOfBaryonFields][vi[2]+i] +
	      InterpolatedField[NumberOfBaryonFields][vi[3]+i] +
	      InterpolatedField[NumberOfBaryonFields][vi[4]+i] +
	      InterpolatedField[NumberOfBaryonFields][vi[5]+i] +
	      InterpolatedField[NumberOfBaryonFields][vi[6]+i] +
	      InterpolatedField[NumberOfBaryonFields][vi[7]+i];

	    BaryonField[Num][index] = InterpolatedField[Num][vi[0]+i] +
	      InterpolatedField[Num][vi[1]+i] + InterpolatedField[Num][vi[2]+i] + 
	      InterpolatedField[Num][vi[3]+i] + InterpolatedField[Num][vi[4]+i] + 
	      InterpolatedField[Num][vi[5]+i] + InterpolatedField[Num][vi[6]+i] + 
	      InterpolatedField[Num][vi[7]+i];

	    if (weight > 0)

	      BaryonField[Num][index] /= weight;
	    
	  } // ENDFOR i
	} // ENDFOR j
      
    } // ENDELSE subgrid marker
#endif /* TRANSFER */
  
  return SUCCESS;
}
