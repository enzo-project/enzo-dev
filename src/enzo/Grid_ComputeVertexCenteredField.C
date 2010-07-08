/***********************************************************************
/
/  GRID CLASS (COMPUTE A VERTEX CENTERED FIELD)
/
/  written by: Tom Abel
/  date:       July, 2007
/  modified1:
/
/  PURPOSE: Compute a vertex centered field and stores in 
/           InterpolatedField[].
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

int grid::ComputeVertexCenteredField(int Num)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Return if the interpolated field is already calculated */

  if (InterpolatedField[Num] != NULL)
    return SUCCESS;

  /* declarations */

  int i,j,k, index, dim, size = 1;

  /* Error Check */

  if (BaryonField[Num] == NULL) {
    ENZO_FAIL("grid::ComputeVertexCenteredField called with inconsistent "
	    "field number");
  }

  /* Compute the size of the new grid adds one cell to the active region */

  for (dim = 0; dim < GridRank; dim++)
    size *= GridEndIndex[dim]-GridStartIndex[dim]+2 ;

  InterpolatedField[Num] = new float[size];

  int ci = 0;
  int vi[8], v;

  if (GridRank == 1) 
    for (i = 0; i <= (GridEndIndex[0]-GridStartIndex[0]+1); i++) 
      InterpolatedField[Num][i] = 0.5 * 
	(BaryonField[Num][GRIDINDEX(i,0,0)] + 
	 BaryonField[Num][GRIDINDEX(i-1,0,0)]);

  if (GridRank == 2) 
    for (j = 0; j <=(GridEndIndex[1]-GridStartIndex[1]+1); j++) 
      for (i = 0; i <=(GridEndIndex[0]-GridStartIndex[0]+1); i++) 
	InterpolatedField[Num][ci++] = 0.25 *
	  (BaryonField[Num][GRIDINDEX(i,j,0)] + 
	   BaryonField[Num][GRIDINDEX(i-1,j,0)] +
	   BaryonField[Num][GRIDINDEX(i,j-1,0)] + 
	   BaryonField[Num][GRIDINDEX(i-1,j-1,0)]);

  if (GridRank == 3) 
    for (k = 0; k <=(GridEndIndex[2]-GridStartIndex[2]+1); k++) 
      for (j = 0; j <=(GridEndIndex[1]-GridStartIndex[1]+1); j++) {
	index = GRIDINDEX(-1,j-1,k-1);
	this->NeighborIndices(index, vi);
	for (i = 0; i <=(GridEndIndex[0]-GridStartIndex[0]+1); i++) {
	  InterpolatedField[Num][ci++] = 0.125f *
	    (BaryonField[Num][vi[0]+i] + BaryonField[Num][vi[1]+i] +
	     BaryonField[Num][vi[2]+i] + BaryonField[Num][vi[3]+i] +
	     BaryonField[Num][vi[4]+i] + BaryonField[Num][vi[5]+i] +
	     BaryonField[Num][vi[6]+i] + BaryonField[Num][vi[7]+i]);
	} // ENDFOR i
      } // ENDFOR j

  return SUCCESS;
}

/***********************************************************************/

float grid::ComputeInterpolatedValue(int Num, int vci, int vcj, int vck, 
				     float mx, float my, float mz)
{

  if (InterpolatedField[Num] == NULL)

    return FLOAT_UNDEFINED;

  float *vc = InterpolatedField[Num];
  int index, vi[8];
  index = VCGRIDINDEX(vci,vcj,vck);
  this->NeighborVertexIndices(index, vi);

  return
    vc[vi[0]] * (1-mx)*(1-my)*(1-mz) +
    vc[vi[1]] * (mx)  *(1-my)*(1-mz) +
    vc[vi[2]] * (1-mx)*(my)  *(1-mz) +
    vc[vi[3]] * (mx)  *(my)  *(1-mz) +
    vc[vi[4]] * (1-mx)*(1-my)*(mz)   +
    vc[vi[5]] * (mx)  *(1-my)*(mz)   +
    vc[vi[6]] * (1-mx)*(my)  *(mz)   +
    vc[vi[7]] * (mx)  *(my)  *(mz);


}

/************************************************************************/

int grid::NeighborIndices(int index, int vi[])
{

  vi[0] = index;					// BOTTOM-SW (ORIGIN)
  vi[1] = index+1;					// BOTTOM-SE (i+1,j,k)
  vi[2] = index+GridDimension[0];			// BOTTOM-NW (i,j+1,k)
  vi[3] = vi[2]+1;					// BOTTOM-NE (i+1,j+1,k)
  vi[4] = index+GridDimension[0]*GridDimension[1];	// TOP-SW (i,j,k+1)
  vi[5] = vi[4]+1;					// TOP-SE (i+1,j,k+1)
  vi[6] = vi[4]+GridDimension[0];			// TOP-NW (i,j+1,k+1)
  vi[7] = vi[6]+1;					// TOP-NE (i+1,j+1,k+1)

  return SUCCESS;

}

/************************************************************************/

int grid::NeighborVertexIndices(int index, int vi[])
{

  int dim, vcdim[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    vcdim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 2;

  vi[0] = index;			// BOTTOM-SW (ORIGIN)
  vi[1] = index+1;			// BOTTOM-SE (i+1,j,k)
  vi[2] = index+vcdim[0];		// BOTTOM-NW (i,j+1,k)
  vi[3] = vi[2]+1;			// BOTTOM-NE (i+1,j+1,k)
  vi[4] = index+vcdim[0]*vcdim[1];	// TOP-SW (i,j,k+1)
  vi[5] = vi[4]+1;			// TOP-SE (i+1,j,k+1)
  vi[6] = vi[4]+vcdim[0];		// TOP-NW (i,j+1,k+1)
  vi[7] = vi[6]+1;			// TOP-NE (i+1,j+1,k+1)

  return SUCCESS;

}
