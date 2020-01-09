/***********************************************************************
/
/  GRID CLASS (ADDS PARENTAL POTENTIAL FIELD)
/  based on Grid_AddParentAccelerationField.C
/
/  written by: JC Passy
/  date:       March, 2012
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "communication.h"

/* function prototypes */

extern "C" void FORTRAN_NAME(prolong_add)(float *source, float *dest, int *ndim,
					  int *sdim1, int *sdim2, int *sdim3,
					  int *ddim1, int *ddim2, int *ddim3,
					  int *start1, int *start2, int *start3,
					  int *refine1, int *refine2,int *refine3);

int grid::AddParentPotentialFieldAPM(grid *ParentGrid)
{

  if (ParentGrid == NULL)
    return FAIL;

  /* Return if this doesn't involve us. */

  if (this->CommunicationMethodShouldExit(ParentGrid))
    return SUCCESS;

  /* Allocate space for potential field, if necessary. */

  int dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (MyProcessorNumber == ProcessorNumber &&
      CommunicationDirection != COMMUNICATION_POST_RECEIVE)
    for (dim = 0; dim < GridRank; dim++)
      if (PotentialField == NULL)
	ENZO_FAIL("PotentialField must already exist.\n")

  /* Declarations. */

  int ParentOffset[MAX_DIMENSION], ParentStartIndex[MAX_DIMENSION],
      ParentTempDim[MAX_DIMENSION], Refinement[MAX_DIMENSION],
      ParentDim[MAX_DIMENSION];

  /* Compute refinement factors. */

  ParentGrid->ComputeRefinementFactors(this, Refinement);

  /* Set unused dims to zero. */

  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    ParentOffset[dim] = ParentStartIndex[dim] = 0;
    ParentTempDim[dim] = ParentDim[dim] = 1;
  }

  /* Compute the ParentOffset (in grid units) and ParentStartIndex and
     the region dim (in parent units). */

  /*
  for (dim = 0; dim < GridRank; dim++) {
    ParentOffset[dim] = nint((CellLeftEdge[dim][0] -
			      ParentGrid->CellLeftEdge[dim][0])/CellWidth[dim][0]);
    ParentStartIndex[dim] = ParentOffset[dim]/Refinement[dim] - 1;
    ParentTempDim[dim] = (ParentOffset[dim] + GridDimension[dim] - 1)/Refinement[dim] -
			 ParentStartIndex[dim] + 2;
    ParentDim[dim] = ParentGrid->GridDimension[dim];
    if (ParentStartIndex[dim] < 0 ||
	ParentStartIndex[dim]+ParentTempDim[dim] > ParentDim[dim]) {
      fprintf(stderr, "ParentStartIndex[%d] = %d ParentTempDim = %d(%d).\n",
	      dim, ParentStartIndex[dim], ParentTempDim[dim], ParentDim[dim]);
      fprintf(stderr, "GridDimension = %d  ParentGridDimension = %d  CLE = %g PCLE = %g\n",
	      GridDimension[dim], ParentDim[dim], CellLeftEdge[dim][0], ParentGrid->CellLeftEdge[dim][0]);
      return FAIL;
    }
  }
  */

  // from Grid_PreparePotentialField.C

  for (dim = 0; dim < GridRank; dim++) {
    ParentOffset[dim] = nint((GravitatingMassFieldLeftEdge[dim] -
			      ParentGrid->GravitatingMassFieldLeftEdge[dim])/
                             GravitatingMassFieldCellSize);
    ParentStartIndex[dim] = ParentOffset[dim]/Refinement[dim] - 1;
    ParentTempDim[dim] = (ParentOffset[dim] +
			  GravitatingMassFieldDimension[dim]-1)/Refinement[dim] -
      ParentStartIndex[dim] + 3;
    ParentDim[dim] = ParentGrid->GravitatingMassFieldDimension[dim];
    if (ParentStartIndex[dim] < 0 ||
        ParentStartIndex[dim]+ParentTempDim[dim] > ParentDim[dim]) {
      ENZO_VFAIL("ParentStartIndex[%"ISYM"] = %"ISYM" ParentTempDim = %"ISYM"(%"ISYM").\n",
		 dim, ParentStartIndex[dim], ParentTempDim[dim], ParentDim[dim])
	}
  }

  /* Copy data from other processor if needed (modify ParentDim and
     ParentStartIndex to reflect the fact that we are only coping part of
     the grid. */

  /* If posting a receive, then record details of call. */

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE &&
      ProcessorNumber != ParentGrid->ProcessorNumber) {
    CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
    CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = ParentGrid;
    CommunicationReceiveCallType[CommunicationReceiveIndex] = 24;
  }
#endif /* USE_MPI */

  if (ProcessorNumber != ParentGrid->ProcessorNumber) {
    ParentGrid->CommunicationSendRegion(ParentGrid, ProcessorNumber,
					POTENTIAL_FIELD, NEW_ONLY, ParentStartIndex, ParentTempDim);
    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
        CommunicationDirection == COMMUNICATION_SEND) {

      return SUCCESS;
    }
    for (dim = 0; dim < GridRank; dim++) {
      ParentOffset[dim] -= Refinement[dim]*ParentStartIndex[dim];
      ParentDim[dim] = ParentTempDim[dim];
    }
  }

  /* Return if this is not our concern. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Interpolate */
  FORTRAN_NAME(prolong_add)(ParentGrid->PotentialField,
			    PotentialField, &GridRank,
			    ParentDim, ParentDim+1, ParentDim+2,
			    GravitatingMassFieldDimension,
			    GravitatingMassFieldDimension+1,
			    GravitatingMassFieldDimension+2,
			    ParentOffset, ParentOffset+1, ParentOffset+2,
			    Refinement, Refinement+1, Refinement+2);

  /* Clean up parent. */

  if (MyProcessorNumber != ParentGrid->ProcessorNumber) {
    delete [] ParentGrid->PotentialField;
    ParentGrid->PotentialField = NULL;
  }

  return SUCCESS;

}
