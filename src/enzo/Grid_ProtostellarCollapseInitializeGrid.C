/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE PROTOSTELLAR CORE COLLAPSE TEST) 
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, June 2005.
/
/  PURPOSE: Sets the initial core region.
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

int grid::ProtostellarCollapseInitializeGrid(
	    float ProtostellarCollapseCoreDensity,
            float ProtostellarCollapseCoreEnergy,
            float ProtostellarCollapseCoreRadius,
	    float ProtostellarCollapseAnguarVelocity)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* set fields in the rotating core region: x^2+y^2+z^2 < dr^2. */

  int index, jndex, kndex, i;
  float zonex, zoney, zonez;
  float r2 = ProtostellarCollapseCoreRadius*ProtostellarCollapseCoreRadius;
  float xcenter = (DomainRightEdge[0] - DomainLeftEdge[0])/2.0;
  float ycenter = (DomainRightEdge[1] - DomainLeftEdge[1])/2.0;
  float zcenter = (DomainRightEdge[2] - DomainLeftEdge[2])/2.0;
  for (i = 0; i < size; i++) {

    index = i % GridDimension[0];
    jndex = ((i-index) % (GridDimension[0]*GridDimension[1]))/GridDimension[0];
    kndex = (i-index-jndex*GridDimension[0])/GridDimension[0]/GridDimension[1];

    zonex = *(CellLeftEdge[0] + index) + 0.5*(*(CellWidth[0] + index)) - xcenter;
    zoney = *(CellLeftEdge[1] + jndex) + 0.5*(*(CellWidth[1] + jndex)) - ycenter;
    zonez = *(CellLeftEdge[2] + kndex) + 0.5*(*(CellWidth[2] + kndex)) - zcenter;;

    if (zonex*zonex + zoney*zoney + zonez*zonez < r2) {
      BaryonField[0][i] = ProtostellarCollapseCoreDensity;
      // thermal energy
      BaryonField[1][i] = ProtostellarCollapseCoreEnergy;
      // solid body clockwise rotation
      BaryonField[2][i] =   ProtostellarCollapseAnguarVelocity*zoney;
      BaryonField[3][i] = - ProtostellarCollapseAnguarVelocity*zonex;
      BaryonField[4][i] = 0.0; // the rotation axis is || to z.
      // total energy
      BaryonField[1][i] += (BaryonField[2][i]*BaryonField[2][i] +
			    BaryonField[3][i]*BaryonField[3][i])/2.0;
    }
  }

  return SUCCESS;
}
