/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (SETS OUTFLOW BOUNDARY CONDITIONS FOR DOUBLE MACH)
/
/  written by: Greg Bryan
/  date:       Marc, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
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
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
int ExternalBoundary::SetDoubleMachBoundary(FLOAT time, FLOAT CellLeftEdge[],
					    FLOAT CellWidth[])
{
  /* declarations */
 
  int i, k, field, index;
  float x, xx = 1.0/6.0 + (1.0 + 20.0*time)/sqrt(3.0);
  boundary_type tmp;
 
  for (k = 0; k < BoundaryDimension[2]; k++) {
    index = k*BoundaryDimension[0];
    for (i = 0; i < BoundaryDimension[0]; i++) {
 
      /* Set right side of Dimension 1. */
 
      x = CellLeftEdge[i] + 0.5*CellWidth[i];
      if (x < xx)
	tmp = inflow;
      else
	tmp = outflow;
      for (field = 0; field < NumberOfBaryonFields; field++)
	BoundaryType[field][1][1][index+i] = tmp;
 
      /* Set left side of Dimension 1. */
 
      if (x < 1.0/6.0)
	tmp = inflow;
      else
	tmp = reflecting;
      for (field = 0; field < NumberOfBaryonFields; field++)
	BoundaryType[field][1][0][index+i] = tmp;
 
    }
  }
 
  return SUCCESS;
 
}
