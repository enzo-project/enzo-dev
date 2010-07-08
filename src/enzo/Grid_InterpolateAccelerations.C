/***********************************************************************
/
/  GRID CLASS (INTERPOLATE PARTICLE AND GRID ACCELERATIONS)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "communication.h"
 
/* function prototypes */
 
extern "C" void FORTRAN_NAME(int_grid_cic)(float *source,
				   float *dest, int *ndim,
				   int *sdim1, int *sdim2, int *sdim3,
				   int *ddim1, int *ddim2, int *ddim3,
				   float *start1, float *start2, float *start3,
				   int *refine1, int *refine2, int *refine3);
 
 
/* EvolveHierarchy function */
 
int grid::InterpolateAccelerations(grid *FromGrid)
{
 
  /* Return if this grid is not on this processor. */
 
  if (this->CommunicationMethodShouldExit(FromGrid))
    return SUCCESS;
 
  /* Declarations. */
 
  float GridOffset[MAX_DIMENSION];
  int GridStart[MAX_DIMENSION], GridEnd[MAX_DIMENSION], GridDim[MAX_DIMENSION],
      GridActiveDim[MAX_DIMENSION], Refinement[MAX_DIMENSION], dim, size = 1;
 
  /* Compute refinement factors. */
 
  FromGrid->ComputeRefinementFactors(this, Refinement);
 
  /* Set unused dims to zero. */
 
  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    GridOffset[dim] = 0;
    GridStart[dim] = GridEnd[dim] = 0;
    GridDim[dim] = GridActiveDim[dim] = 1;
  }
 
  /* Compute the GridOffset (in grid units) and GridStartIndex and
     the region dim (in grid units). */
 
  for (dim = 0; dim < GridRank; dim++) {
    GridOffset[dim] = (CellLeftEdge[dim][0] -
		       FromGrid->CellLeftEdge[dim][0])/
		       CellWidth[dim][0];
    GridDim[dim]    = GridDimension[dim];
#ifdef UNUSED
    GridStart[dim]  = nint(GridOffset[dim]/Refinement[dim]) - 1;
    GridEnd[dim]    = nint((GridStart[dim] + GridDim[dim])/Refinement[dim])+2;
#endif /* UNUSED */
    GridStart[dim]  = 0;
    GridEnd[dim]    = FromGrid->GridDimension[dim] - 1;
    GridActiveDim[dim] = GridEnd[dim] - GridStart[dim] + 1;
    size *= GridDimension[dim];
 
    if (GridOffset[dim] < 0) {
      ENZO_VFAIL("GridOffset[%"ISYM"] = %"GSYM" < 0.\n", dim, GridOffset[dim])
    }
  }
 
  /* If posting a receive, then record details of call. */

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE &&
      ProcessorNumber != FromGrid->ProcessorNumber) {
    CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
    CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = FromGrid;
    CommunicationReceiveCallType[CommunicationReceiveIndex] = 10;
  }
#endif /* USE_MPI */

  /* Copy data from other processor if needed (modify GridDim and
     GridStartIndex to reflect the fact that we are only coping part
     of the grid. */
 
  if (ProcessorNumber != FromGrid->ProcessorNumber) {
    FromGrid->CommunicationSendRegion(FromGrid, ProcessorNumber,
	    ACCELERATION_FIELDS, NEW_ONLY, GridStart, GridActiveDim);
    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	CommunicationDirection == COMMUNICATION_SEND)
      return SUCCESS;
#ifdef UNUSED
    for (dim = 0; dim < GridRank; dim++) {
      GridOffset[dim] = (CellLeftEdge[dim][GridStart[dim]] -
		       FromGrid->GravitatingMassFieldLeftEdge[dim])/
		       CellWidth[dim][0];
      GridDim[dim] = GridEnd[dim] - GridStart[dim] + 1;
      GridStart[dim] = 0;
      GridEnd[dim] = GridDim[dim] - 1;
    }
#endif /* UNUSED */
  }
 
  /* Return if this is not our concern. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (FromGrid->AccelerationField[0] == NULL) {
        ENZO_FAIL("FromGrid->AccelerationField is NULL.");
  }
 
  /* Allocate acceleration fields. */
 
  for (dim = 0; dim < GridRank; dim++) {
    delete [] AccelerationField[dim];
    AccelerationField[dim] = new float[size];
  }
 
  /* --------------------------------------------------- */
  /* Interpolate accelerations to particles.             */
 
  if (NumberOfParticles > 0) {
 
    /* Clear particle accelerations. */
 
    this->ClearParticleAccelerations();
 
    /* Move particles 1/2 step forward in preparation for interpolation. */
 
    this->UpdateParticlePosition(0.5*dtFixed);
 
    /* Interpolate the accelerations back to the grid and particles.
       (this assumes that FromGrid's acdeleration fields have been copied
        from the other processor in their entirety).  Also, note hack
       to prevent error. */
 
    int hold = FromGrid->ProcessorNumber;
    FromGrid->ProcessorNumber = ProcessorNumber;
    int DiffType = DIFFERENCE_TYPE_NORMAL;
    this->InterpolateParticlePositions(FromGrid, DiffType);
    FromGrid->ProcessorNumber = hold;
 
    /* Move particles 1/2 step backwards to return to original positions. */
 
    this->UpdateParticlePosition(-0.5*dtFixed);
 
  } // end: if (NumberOfParticles > 0)
 
  /* --------------------------------------------------- */
  /* Interpolate acceleration field for cells. */
 
  if (NumberOfBaryonFields > 0)
    for (dim = 0; dim < GridRank; dim++) {
 
      /* If using zeus, then shift the offset since accelerations are
	 face centered. */
 
      if (HydroMethod == Zeus_Hydro)
	GridOffset[dim] += 0.5*(FromGrid->CellWidth[dim][0] -
				CellWidth[dim][0])/CellWidth[dim][0];
 
      /* Call a fortran routine to do the work. */
 
      FORTRAN_NAME(int_grid_cic)(FromGrid->AccelerationField[dim],
				 AccelerationField[dim], &GridRank,
				 FromGrid->GridDimension,
				 FromGrid->GridDimension+1,
				 FromGrid->GridDimension+2,
				 GridDimension, GridDimension+1,
				     GridDimension+2,
				 GridOffset, GridOffset+1, GridOffset+2,
				 Refinement, Refinement+1, Refinement+2);
 
      /* Shift back, if necessary. */
 
      if (HydroMethod == Zeus_Hydro)
	GridOffset[dim] -= 0.5*(FromGrid->CellWidth[dim][0] -
				CellWidth[dim][0])/CellWidth[dim][0];
 
    } // ENDFOR dims

  /* Clean up if we have transfered data. */

  if (MyProcessorNumber != FromGrid->ProcessorNumber)

    FromGrid->DeleteAllFields();
 
  return SUCCESS;
}
