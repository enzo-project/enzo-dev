/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE PRESENCE OF SHEAR)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:  Alexei Kritsuk, Aug.2004
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
 
int grid::FlagCellsToBeRefinedByShear()
{
  /* declarations */
 
  int i, j, k, index, dim;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Allocate temp array and fill with zeroes */
 
  float *DelVelocity = new float[size];
  for (i = 0; i < size; i++)
    DelVelocity[i] = 0.0;
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* loop over active dimensions */
 
  int Offset = 1;
  for (dim = 0; dim < GridRank; dim++)
    if (GridDimension[dim] > 1) {
 
      /* For shear: [(du/dy)^2 + (dv/dx)^2 + (dw/dy)^2 +
                     (du/dz)^2 + (dv/dz)^2 + (dw/dx)^2  ] > parameter,
	 where:  du/dy = [u(j-1) - u(j+1)]/[2dy],
	 assume: dx=dy=dz=1;
	 parameter ~ (sound/dx)^2
         assume: sound = 1
      */
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
 
	    index = i + j*GridDimension[0] +
	                k*GridDimension[1]*GridDimension[0];
 
	    float DelVel1 = 0.0, DelVel2 = 0.0;
	    switch (dim) {
	    case 0:
	      if (GridRank > 1)
		DelVel1 = BaryonField[Vel1Num+1][index - Offset] -
		          BaryonField[Vel1Num+1][index + Offset];
	      if (GridRank > 2)
		DelVel2 = BaryonField[Vel1Num+2][index - Offset] -
		          BaryonField[Vel1Num+2][index + Offset];
	      break;
 
	    case 1:
	      DelVel1 = BaryonField[Vel1Num][index - Offset] -
		        BaryonField[Vel1Num][index + Offset];
	      if (GridRank > 2)
		DelVel2 = BaryonField[Vel1Num+2][index - Offset] -
		          BaryonField[Vel1Num+2][index + Offset];
	      break;
 
	    case 2:
	      if (GridRank > 1)
		DelVel1 = BaryonField[Vel1Num+1][index - Offset] -
		          BaryonField[Vel1Num+1][index + Offset];
	      DelVel2 = BaryonField[Vel1Num][index - Offset] -
		        BaryonField[Vel1Num][index + Offset];
	      break;
 
	    default:
	      break;
	    }
 
	    DelVel1 *= DelVel1;
	    DelVel2 *= DelVel2;
	    DelVelocity[index] += DelVel1 + DelVel2;
	    if (dim == GridRank-1)
	      FlaggingField[index] +=
		(DelVelocity[index] > MinimumShearForRefinement) ? 1 : 0;
	  }
 
      Offset *= GridDimension[dim];
 
    }  // end loop over dimension
 
  /* clean up */
 
  delete DelVelocity;
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i])

      NumberOfFlaggedCells++;
 
  return NumberOfFlaggedCells;
 
}
