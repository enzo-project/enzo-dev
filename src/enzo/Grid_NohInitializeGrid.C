/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR NOH PROBLEM)
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:  Alexei Kritsuk, May 2005
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
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

int grid::NohInitializeGrid(float d0, float p0, float u0)
{
  /* declarations */

  int index, size, dim, i, j, k;

  /* error check */

  if (GridRank < 2 || GridRank > 3) {
    ENZO_FAIL("GridRank must be 2 or 3\n");
  }

  /* create fields */

  NumberOfBaryonFields = 2 + GridRank;
  FieldType[0] = Density;
  FieldType[1] = TotalEnergy;
  FieldType[2] = Velocity1;
  FieldType[3] = Velocity2;
  if (GridRank > 2)
    FieldType[4] = Velocity3;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  this->AllocateGrids();

  /* set fields */

  float radius, xx, yy, zz = 0.0;
  int ishift = GridStartIndex[0], jshift = GridStartIndex[1];
  int kshift = GridStartIndex[2];

  if (NohProblemFullBox == 1) {
    if (GridDimension[0]%2 != 0 || GridDimension[1]%2 != 0)
      ERROR_MESSAGE;
    if (GridRank > 2)
      if (GridDimension[2]%2 != 0)
	ERROR_MESSAGE;
    ishift = GridDimension[0]/2;
    jshift = GridDimension[1]/2;
    if (GridRank > 2)
      kshift = GridDimension[2]/2;
  }

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++) {
      index = k*GridDimension[0]*GridDimension[1] + j*GridDimension[0];
      for (i = 0; i < GridDimension[0]; i++) {
	BaryonField[0][index+i]  = d0;
	BaryonField[1][index+i]  = p0/(Gamma-1.0)/d0;
	xx = i + 0.5 - ishift;
	yy = j + 0.5 - jshift;
	if (GridRank > 2)
	  zz = k + 0.5 - kshift;
	radius = sqrt(xx*xx + yy*yy + zz*zz);
	BaryonField[2][index+i]  = u0*xx/radius;
	BaryonField[3][index+i]  = u0*yy/radius;
	if (GridRank > 2)
	  BaryonField[4][index+i]  = u0*zz/radius;
	if (HydroMethod != Zeus_Hydro)

	  BaryonField[1][index+i] += 0.5*u0*u0;
      }
    }
  
  return SUCCESS;
}
