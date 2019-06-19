/***********************************************************************
/
/  FIND THE TOTAL MASS IN THE CONTROL REGION
/
/  written by: John Regan
/  date:       February 2017
/  modified1:
/
/  PURPOSE: get the total mass
/
************************************************************************/
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"

float grid::FindMassinRegion(FLOAT *cellpos, FLOAT radius, int DensNum)
{

  int i = 0, j = 0, k = 0;
  FLOAT pos[3];
  FLOAT dx = CellWidth[0][0];
  float TotalMass = 0.0;
  int index = 0;
  float *density = BaryonField[DensNum];
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i, j, k);
	pos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	pos[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	/* Use only values within the control volume */
	if(sqrt((cellpos[0] - pos[0])*(cellpos[0] - pos[0]) +
		(cellpos[1] - pos[1])*(cellpos[1] - pos[1]) +
		(cellpos[2] - pos[2])*(cellpos[2] - pos[2])) < radius)
	{
	  TotalMass += density[index]*dx*dx*dx;
	}
      }
    }
  }

  return TotalMass;
}

float grid::FindMassinGrid(int DensNum)
{
  int i = 0, j = 0, k = 0;
  FLOAT dx = CellWidth[0][0];
  float TotalMass = 0.0;
  int index = 0;
  float *density = BaryonField[DensNum];
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i, j, k);
	TotalMass += density[index]*dx*dx*dx;
      }
    }
  }

  return TotalMass;
}

