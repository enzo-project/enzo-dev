/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID TO A SHOCK TUBE PROBLEM)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::ShockTubeInitializeGrid(int ShockTubeDirection,
				  float ShockTubeBoundary,
				  float ShockTubeDensity[],
				  float ShockTubePressure[],
				  float ShockTubeVelocity[])
{
 
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
 
  /* declarations */
 
  int Divisor, index, side, size, dim, field, i;
 
  /* error check */
 
  if (ShockTubeDirection > GridRank-1 ||
      GridDimension[ShockTubeDirection] == 1) {
    fprintf(stderr, "ShockTubeDirection is not properly defined.\n");
    return FAIL;
  }
 
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  for (field = 0; field < NumberOfBaryonFields; field++)
    BaryonField[field] = new float[size];
 
  /* set fields */
 
  Divisor = 1;
  for (dim = 0; dim < ShockTubeDirection; dim++)
    Divisor *= GridDimension[dim];
 
  /* set density, total energy and velocity in problem dimension */
 
  for (i = 0; i < size; i++) {
    index = i/Divisor % GridDimension[ShockTubeDirection];
    side = (*(CellLeftEdge[ShockTubeDirection] + index) +
	    0.5*(*(CellWidth[ShockTubeDirection] + index)) > ShockTubeBoundary)
      ? 1 : 0;
    BaryonField[0][i] = ShockTubeDensity[side];
    BaryonField[1][i] = ShockTubePressure[side] /
                  ((Gamma - 1.0)*ShockTubeDensity[side]  ) +
                  0.5*ShockTubeVelocity[side]*ShockTubeVelocity[side];
    BaryonField[2+ShockTubeDirection][i] = ShockTubeVelocity[side];
  }
 
  /* set transverse velocities */
 
  for (field = 2; field < 2+GridRank; field++)
    if (field != ShockTubeDirection+2)
      for (i = 0; i < size; i++)
	BaryonField[field][i] = 0.0;
 
  return SUCCESS;
}
