/***********************************************************************
/
/  GRID CLASS (INTERPOLATE POSITIONS FROM THE ACCELERATION FIELD)
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
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
 
extern "C" void PFORTRAN_NAME(cic_interp)(FLOAT *posx, FLOAT *posy,
			FLOAT *posz, int *ndim, int *npositions,
                        float *sumfield, float *field, FLOAT *leftedge,
                        int *dim1, int *dim2, int *dim3, FLOAT *cellsize);
 
int grid::InterpolatePositions(FLOAT *Position[], int dim, float *Field,
			       int Number)
{
  if (Number == 0 || MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* Set the pointer to the AccelerationField or the PotentialField. */
 
  float *InterpolationField = AccelerationField[dim];
  if (dim == GridRank)
    InterpolationField = PotentialField;
 
  /* Error check. */
 
  if (InterpolationField == NULL) {
    ENZO_VFAIL("AccelerationField[%"ISYM"] absent.\n", dim)
  }
 
  if (GravitatingMassFieldCellSize <= 0) {
    ENZO_FAIL("GravitatingMassFieldCellSize undefined.\n");

  }
 
  /* Set the left edge of the field. */
 
  FLOAT LeftEdge[MAX_DIMENSION];
  for (int i = 0; i < GridRank; i++)
    LeftEdge[i] = CellLeftEdge[i][0];
//    LeftEdge[i] = CellLeftEdge[i][0] - ((dim == i)? (0.5*CellWidth[i][0]) : 0);
 
  /* Interpolate from field. */
 
  PFORTRAN_NAME(cic_interp)(Position[0], Position[1], Position[2], &GridRank,
			   &Number, Field, InterpolationField, LeftEdge,
			   GridDimension, GridDimension+1, GridDimension+2,
			   &GravitatingMassFieldCellSize);
 
  return SUCCESS;
}
