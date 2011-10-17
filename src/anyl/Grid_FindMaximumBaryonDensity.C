/***********************************************************************
/
/  GRID CLASS (FIND MAXIMUM BARYON DENSITY)
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../enzo/macros_and_parameters.h"
#include "../enzo/typedefs.h"
#include "../enzo/global_data.h"
#include "../enzo/Fluxes.h"
#include "../enzo/GridList.h"
#include "../enzo/ExternalBoundary.h"
#include "../enzo/Grid.h"

/* function prototypes */



int grid::FindMaximumBaryonDensity(FLOAT Position[MAX_DIMENSION],
				   float *MaxDensity)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  
  int i, j, k, index;

  /*
  if (ComovingCoordinates != 1) {  
    InitialRedshift = 0; 
    //FinalRedshift = 0;
    HubbleConstantNow = 0.7; 
    OmegaMatterNow = 0.3;
    OmegaLambdaNow = 0.7;
    //    float ComovingBoxSize = 1;
    //    float MaxExpansionRate = 1;
  } 
  */

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);

  /* Loop over grid. */

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	if (BaryonField[DensNum][index] > *MaxDensity) {
	  *MaxDensity = BaryonField[DensNum][index];
	  Position[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  if (GridRank > 0)
	    Position[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  if (GridRank > 1)
	    Position[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
    }
	}
    }

  return SUCCESS;
}
