/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID TO A DOUBLE MACH REFLECTION PROBLEM)
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:
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
 
int grid::DoubleMachInitializeGrid(float d0, float e0, float u0, float v0,
				   float w0)
{
  /* declarations */
 
  int index, size, dim, field, i, j, k;
  float xx;
 
  /* error check */
 
  if (GridRank < 2) {
    ENZO_FAIL("GridRank must be > 1\n");
  }
 
  /* create fields */
 
  NumberOfBaryonFields = 2 + GridRank;
  FieldType[0] = Density;
  FieldType[1] = TotalEnergy;
  FieldType[2] = Velocity1;
  FieldType[3] = Velocity2;
  FieldType[4] = Velocity3;
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  this->AllocateGrids();
 
  /* set fields */
 
  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0];
      xx = 1.0/6.0 + (CellLeftEdge[1][j] + 0.5*CellWidth[1][j])/sqrt(3.0);
      for (i = 0; i < GridDimension[0]; i++) {
	if (CellLeftEdge[0][i] + 0.5*CellWidth[0][i] < xx) {
	  BaryonField[0][index+i] = d0;
	  BaryonField[1][index+i] = e0/d0 + 0.5*(u0*u0 + v0*v0 + w0*w0);
	  BaryonField[2][index+i] = u0;
	  BaryonField[3][index+i] = v0;
	  if (GridRank > 2) BaryonField[4][index+i] = w0;
	} else {
	  BaryonField[0][index+i] = 1.4;
	  BaryonField[1][index+i] = 2.5/1.4;
	  BaryonField[2][index+i] = 0.0;
	  BaryonField[3][index+i] = 0.0;
	  if (GridRank > 2) BaryonField[4][index+i] = 0.0;

	}
      }
    }
 
  return SUCCESS;
}
