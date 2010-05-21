/***********************************************************************
/
/  SET THE POTENTIAL FIELD BOUNDARY CONDITIONS OF AN ISOLATED GRID
/
/  written by: Greg Bryan
/  date:       April, 2005
/  modified1:
/
/  PURPOSE:  THIS ROUTINE SETS THE POTENTIAL GHOST ZONES BY ASSUMING ZERO
/            GRADIENT.  This is not perfect for the isolated topgrid potential
/            but if the grid is isolated, nothing should be outside so it
/            should be ok.
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

#define GINDEX(i1,i2,i3) (((i3)*GravitatingMassFieldDimension[1]+(i2))*GravitatingMassFieldDimension[0]+(i1))

int grid::SetIsolatedPotentialBoundary()
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (PotentialField == NULL || GravitatingMassFieldCellSize == FLOAT_UNDEFINED) {
    ENZO_FAIL("Potential NULL or gravity unitialized.\n");

  }

  /* Set start index and dimension of active part of potential field. */

  int dim, i, j, k, GravStart[] = {0,0,0}, GravEnd[] = {0,0,0};
  for (dim = 0; dim < GridRank; dim++) {
    GravStart[dim] = nint((GridLeftEdge[dim] -
	       GravitatingMassFieldLeftEdge[dim])/GravitatingMassFieldCellSize);
    GravEnd[dim] = nint((GridRightEdge[dim] -
	       GravitatingMassFieldLeftEdge[dim])/GravitatingMassFieldCellSize)-1;
  }

  /* First, copy potential values along boundaries into i-direction. */

  for (k = GravStart[2]; k <= GravEnd[2]; k++)
    for (j = GravStart[1]; j <= GravEnd[1]; j++) {
      for (i = 0; i < GravStart[0]; i++)
	PotentialField[GINDEX(i,j,k)] = PotentialField[GINDEX(GravStart[0],j,k)];
      for (i = GravEnd[0]+1; i < GravitatingMassFieldDimension[0]; i++)
	PotentialField[GINDEX(i,j,k)] = PotentialField[GINDEX(GravEnd[0],j,k)];
    }


  /* Next copy along the j-direction. */

  for (k = GravStart[2]; k <= GravEnd[2]; k++)
    for (i = 0; i < GravitatingMassFieldDimension[0]; i++) {
      for (j = 0; j < GravStart[1]; j++)
	PotentialField[GINDEX(i,j,k)] = PotentialField[GINDEX(i,GravStart[1],k)];
      for (j = GravEnd[1]+1; j < GravitatingMassFieldDimension[1]; j++)
	PotentialField[GINDEX(i,j,k)] = PotentialField[GINDEX(i,GravEnd[1],k)];
    }

  /* Finally copy along the k-direction. */

  for (j = 0; j < GravitatingMassFieldDimension[1]; j++)
    for (i = 0; i < GravitatingMassFieldDimension[0]; i++) {
      for (k = 0; k < GravStart[2]; k++)
	PotentialField[GINDEX(i,j,k)] = PotentialField[GINDEX(i,j,GravStart[2])];
      for (k = GravEnd[2]+1; k < GravitatingMassFieldDimension[2]; k++)
	PotentialField[GINDEX(i,j,k)] = PotentialField[GINDEX(i,j,GravEnd[2])];
    }

  return SUCCESS;
}

