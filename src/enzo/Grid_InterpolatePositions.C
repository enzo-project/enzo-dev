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

void InterpolatePositionsPileUpTSC1D(FLOAT *Position[], int Number, float *SumField,
                                     float *Field, FLOAT LeftEdge[],
                                     int EffectiveDim[], FLOAT CellSize);
void InterpolatePositionsPileUpTSC2D(FLOAT *Position[], int Number, float *SumField,
                                     float *Field, FLOAT LeftEdge[],
                                     int EffectiveDim[], int Dimension[],
                                     FLOAT CellSize);
void InterpolatePositionsPileUpTSC3D(FLOAT *Position[], int Number, float *SumField,
                                     float *Field, FLOAT LeftEdge[],
                                     int EffectiveDim[], int Dimension[],
                                     FLOAT CellSize);

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

  if (GravitySolverType == GRAVITY_SOLVER_FAST) {

    PFORTRAN_NAME(cic_interp)(Position[0], Position[1], Position[2], &GridRank,
                              &Number, Field, InterpolationField, LeftEdge,
                              GridDimension, GridDimension+1, GridDimension+2,
                              &GravitatingMassFieldCellSize);

  } else if (GravitySolverType == GRAVITY_SOLVER_APM) {

    /* Use TSC with the APM solver. */

    int EffectiveDimension[MAX_DIMENSION], ActualDimension[MAX_DIMENSION];
    FLOAT CellSize = GravitatingMassFieldCellSize;
    int i;

    for (i = 0; i < GridRank; i++) {
      EffectiveDimension[i] = GridDimension[i];
      ActualDimension[i] = GridDimension[i];
    }

    /* 1D case. */
    if (GridRank == 1)
      InterpolatePositionsPileUpTSC1D(Position, Number, Field, InterpolationField,
                                      LeftEdge, EffectiveDimension, CellSize);

    /* 2D Isolated case. */
    if (GridRank == 2)
      InterpolatePositionsPileUpTSC2D(Position, Number, Field, InterpolationField,
                                      LeftEdge, EffectiveDimension,
                                      ActualDimension, CellSize);

    /* 3D Isolated case. */
    if (GridRank == 3)
      InterpolatePositionsPileUpTSC3D(Position, Number, Field, InterpolationField,
                                      LeftEdge, EffectiveDimension,
                                      ActualDimension, CellSize);
  } // end: if (GravitySolverType == GRAVITY_SOLVER_APM)

  return SUCCESS;
}
