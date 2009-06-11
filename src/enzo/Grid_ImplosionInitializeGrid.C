/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR IMPLOSION PROBLEM)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, December 2004.
/
/  PURPOSE: Sets density and total energy in the lower left corner of
/           the domain.
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
 
int grid::ImplosionInitializeGrid(float ImplosionDiamondDensity,
				  float ImplosionDiamondTotalEnergy)
{
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* set fields in a "diamond" region: y=0.15-x. */
 
  int index, jndex, i;
  for (i = 0; i < size; i++) {
    index = i % GridDimension[0];
    jndex = (i-index)/GridDimension[0];
    if (*(CellLeftEdge[0] + index) + 0.5*(*(CellWidth[0] + index))
	     + *(CellLeftEdge[1] + jndex) + 0.5*(*(CellWidth[1] + jndex))
	< 0.1517) { // must be 0.15 but this creates an imperfect initial front
                   // 0.151 is good for 400^2 simulation to fix the
      // imperfection. 0.1517 is good for L2x2 start-up to match the unigrid.
      BaryonField[0][i] = ImplosionDiamondDensity;
      BaryonField[1][i] = ImplosionDiamondTotalEnergy;
    }
  }
 
  return SUCCESS;
}
