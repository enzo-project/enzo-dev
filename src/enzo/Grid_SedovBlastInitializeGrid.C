/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR 2D SEDOV BLAST WAVE TEST) 
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, January 2005.
/
/  PURPOSE: Sets the total energy in the initial explosion region
/           Assumes 2D.
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

int grid::SedovBlastInitializeGrid(float dr,
				   float SedovBlastInnerTotalEnergy)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* set fields in the initial explosion region: x^2+y^2 < dr^2. */

  int index, jndex, i;
  float zonex, zoney, dr2 = dr*dr;
  for (i = 0; i < size; i++) {
    index = i % GridDimension[0];
    jndex = (i-index)/GridDimension[0];
    zonex = *(CellLeftEdge[0] + index) + 0.5*(*(CellWidth[0] + index));
    zoney = *(CellLeftEdge[1] + jndex) + 0.5*(*(CellWidth[1] + jndex));
    zonex -= 0.5*CellWidth[0][0]*(GridEndIndex[0] - GridStartIndex[0] + 1) +
      GridLeftEdge[0];
    zoney -= 0.5*CellWidth[1][0]*(GridEndIndex[1] - GridStartIndex[1] + 1) +
      GridLeftEdge[1];
    if (zonex*zonex + zoney*zoney < dr2)
      BaryonField[1][i] = SedovBlastInnerTotalEnergy;
  }

  return SUCCESS;
}
