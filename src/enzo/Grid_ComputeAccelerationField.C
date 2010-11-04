/***********************************************************************
/
/  GRID CLASS (DIFFERENCES THE POTENTIAL TO GET THE ACCELERATION FIELD)
/
/  written by: Greg Bryan
/  date:       January, 1998
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
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
extern "C" void FORTRAN_NAME(comp_accel)(float *source, float *dest1,
			    float *dest2, float *dest3, int *ndim, int *iflag,
                           int *sdim1, int *sdim2, int *sdim3, int *ddim1,
			    int *ddim2, int *ddim3,
                           int *start1, int *start2, int *start3,
			    float *delx, float *dely, float *delz);
extern "C" void FORTRAN_NAME(smooth)(float *source1, float *source2,
				     float *source3, int *ndim, int *sdim1,
				     int *sdim2, int *sdim3, int *nsmooth);
 
 
int grid::ComputeAccelerationField(int DifferenceType, int level)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (SelfGravity == 0) 
    return SUCCESS;

  /* declarations */
 
  int dim, size = 1, Offset[MAX_DIMENSION] = {0,0,0};
  float CellSize[MAX_DIMENSION] = {1,1,1};
 
  /* Compute adot/a at time = t+1/2dt (time-centered). */
 
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) == FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.\n");
    }
 
  /* Set cell size. */
 
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    CellSize[dim] = a * GravitatingMassFieldCellSize;
    Offset[dim] = nint((CellLeftEdge[dim][0] -
			GravitatingMassFieldLeftEdge[dim])/ CellWidth[dim][0]);
  }
 
  /* Loop over dimensions and difference acceleration. */
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* Allocate acceleration field. */
 
    if (AccelerationField[dim] != NULL) {
      delete [] AccelerationField[dim];
    }
 
    AccelerationField[dim] = new float[size];
 
  }
 
  /* Difference potential. */
 
  FORTRAN_NAME(comp_accel)(PotentialField, AccelerationField[0],
      AccelerationField[1], AccelerationField[2], &GridRank, &DifferenceType,
	    GravitatingMassFieldDimension, GravitatingMassFieldDimension+1,
	      GravitatingMassFieldDimension+2,
	    GridDimension, GridDimension+1, GridDimension+2,
            Offset, Offset+1, Offset+2, CellSize, CellSize+1, CellSize+2);
 
  /* Smooth if necessary. */
 
#define NO_SMOOTH_ACCEL
#ifdef SMOOTH_ACCEL
  int nsmooth = max(level - MaximumGravityRefinementLevel, 0);
  if (nsmooth > 0) {

 
    nsmooth = nint(0.5*POW(RefineBy, nsmooth-1));
    FORTRAN_NAME(smooth)(AccelerationField[0], AccelerationField[1],
			 AccelerationField[2], &GridRank,
			 GridDimension, GridDimension+1, GridDimension+2,
			 &nsmooth);
  }
#endif /* SMOOTH_ACCEL */
 
  return SUCCESS;
}
 
