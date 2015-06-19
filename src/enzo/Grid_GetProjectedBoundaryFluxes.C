/***********************************************************************
/
/  GRID CLASS: PROJECT (DOWNSAMPLE) THIS GRIDS FLUXES BY GIVEN REFINEMENT
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1: David Collins, 2005
/              Updated algebra so Cosmological Expansion is also
/              conservative.  This fix also came with fixes to euler.src and
/              Grid_GetProjectedBoundaryFluxes.C, so make sure you get those.
/
/  PURPOSE:
/
/  NOTE: This routine assumes that the fluxes structure and current grid
/        have the same baryon fields.
/
************************************************************************/
 
// Use the refinement factors in the arguement to down sample the fluxes
//   and put them into the fluxes structure in the argument.  Also,
//   fill out the rest of the fluxes structure.

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
#include "communication.h"
 
int CommunicationSendFluxes(fluxes *Fluxes, int ToProc, int NumberOfFields,
			    int Rank);
int CommunicationReceiveFluxes(fluxes *Fluxes, int FromProc,
			       int NumberOfFields, int Rank);
void DeleteFluxes(fluxes *Fluxes);
 
 
int grid::GetProjectedBoundaryFluxes(grid *ParentGrid, fluxes &ProjectedFluxes)
{
 
  /* Return if this doesn't involve us. */
  
  if (ParentGrid->CommunicationMethodShouldExit(this) ||
      NumberOfBaryonFields == 0)
    return SUCCESS;

  /* Compute the subgrid's refinement factors (for step #19) */

  int RefinementFactors[MAX_DIMENSION];
  ParentGrid->ComputeRefinementFactors(this, RefinementFactors);
 
  int i, j, k, i1, j1, k1, dim, Dims[3], field, size;
  int index1, index2;
  int ProjectedDims[MAX_DIMENSION];
 
  /* If the BoundaryFluxes structure doesn't exist yet, then create it. */
 
  if (BoundaryFluxes == NULL)
    this->PrepareBoundaryFluxes();
 
  /* fill out Flux indicies */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < GridRank; i++) {
      ProjectedFluxes.LeftFluxStartGlobalIndex[dim][i] =
	BoundaryFluxes->LeftFluxStartGlobalIndex[dim][i]/ RefinementFactors[i];
      ProjectedFluxes.LeftFluxEndGlobalIndex[dim][i] =
	BoundaryFluxes->LeftFluxEndGlobalIndex[dim][i] / RefinementFactors[i];
      ProjectedFluxes.RightFluxStartGlobalIndex[dim][i] =
	BoundaryFluxes->RightFluxStartGlobalIndex[dim][i]/RefinementFactors[i];
      ProjectedFluxes.RightFluxEndGlobalIndex[dim][i] =
	BoundaryFluxes->RightFluxEndGlobalIndex[dim][i]/RefinementFactors[i];
    }

  if (CommunicationDirection != COMMUNICATION_POST_RECEIVE) {
 
    /* loop over all dimensions */
 
    for (dim = 0; dim < GridRank; dim++) {
 
      /* compute size of flux region */
 
      Dims[0] = Dims[1] = Dims[2] = 1;
 
      for (i = 0; i < GridRank; i++)
	Dims[i] = BoundaryFluxes->LeftFluxEndGlobalIndex[dim][i] -
	  BoundaryFluxes->LeftFluxStartGlobalIndex[dim][i] + 1;
 
      /* compute ProjectedFlux dimensions (and TotalRefinement). */
 
      int TotalRefinement = 1;
      for (i = 0; i < 3; i++) {
	TotalRefinement *= RefinementFactors[i];
	if (i != dim)
	  ProjectedDims[i] = Dims[i]/RefinementFactors[i];
	else
	  ProjectedDims[i] = 1;
      }
      size = ProjectedDims[0]*ProjectedDims[1]*ProjectedDims[2];
 
      /* compute the fraction of each (current) grid flux cell
         that each subgrid's flux cell occupies.  Er, whatever. */
 
      float dArea = 1.0/float(TotalRefinement);
 
      /* loop over fields */
 
      for (field = 0; field < NumberOfBaryonFields; field++) {
 
	/* Allocate and clear Fluxes */
 
	ProjectedFluxes.LeftFluxes[field][dim] = new float[size];
	ProjectedFluxes.RightFluxes[field][dim] = new float[size];
	for (i = 0; i < size; i++) {
	  ProjectedFluxes.LeftFluxes[field][dim][i] = 0.0;
	  ProjectedFluxes.RightFluxes[field][dim][i] = 0.0;
	}
 
	/* if this dim is of length 0, then there is no Flux. */
	
	if (GridDimension[dim] > 1 && MyProcessorNumber == ProcessorNumber) {
 
	  /* project (downsample by RefinementFactors[i] Fluxes */
 
	  for (k = 0; k < Dims[2]; k++) {
	    k1 = k/RefinementFactors[2];
	    for (j = 0; j < Dims[1]; j++) {
	      j1 = j/RefinementFactors[1];
	      for (i = 0; i < Dims[0]; i++) {
		i1 = i/RefinementFactors[0];
		*(ProjectedFluxes.LeftFluxes[field][dim] +
		  i1 + j1*ProjectedDims[0] +
		  k1*ProjectedDims[0]*ProjectedDims[1]) +=
		  (*(BoundaryFluxes->LeftFluxes[field][dim] +
		     i + j*Dims[0] + k*Dims[0]*Dims[1])) * dArea;
		*(ProjectedFluxes.RightFluxes[field][dim] +
		  i1 + j1*ProjectedDims[0] +
		  k1*ProjectedDims[0]*ProjectedDims[1]) +=
		  (*(BoundaryFluxes->RightFluxes[field][dim] +
		     i + j*Dims[0] + k*Dims[0]*Dims[1])) * dArea;
	      }
	    }
	  }
 
	}  // end: if Dims[dim] > 1
 
      }  // next field
 
      /* set unused flux pointers to null (makes cleanup easier) */
 
      for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	   field++) {
	ProjectedFluxes.LeftFluxes[field][dim]  = NULL;
	ProjectedFluxes.RightFluxes[field][dim] = NULL;
      }
 
    }  // next dimension
 
    /* set unused pointers to NULL */
 
    for (dim = GridRank; dim < MAX_DIMENSION; dim++)
      for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
	ProjectedFluxes.LeftFluxes[field][dim]  = NULL;
	ProjectedFluxes.RightFluxes[field][dim] = NULL;
      }

  } // ENDIF !COMMUNICATION_POST_RECEIVE
 
  /* If posting a receive, then record details of call. */

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
    CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
    CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = ParentGrid;
    CommunicationReceiveCallType[CommunicationReceiveIndex] = 11;
  }
#endif /* USE_MPI */

  /* If appropriate, receive data and exit. */

  if (ProcessorNumber != MyProcessorNumber) {
    if (CommunicationReceiveFluxes(&ProjectedFluxes, ProcessorNumber,
				   NumberOfBaryonFields, GridRank) == FAIL) {
      ENZO_FAIL("Error in CommunicationReceiveFluxes.\n");
    }
    return SUCCESS;
  }

  /* Send fluxes and delete. */
  
  if (ParentGrid->ProcessorNumber != ProcessorNumber) {
    if (CommunicationSendFluxes(&ProjectedFluxes, 
				ParentGrid->ProcessorNumber,
				NumberOfBaryonFields, GridRank) == FAIL) {
      ENZO_FAIL("Error in CommunicationSendFluxes.\n");

    }
    DeleteFluxes(&ProjectedFluxes);
  }

  return SUCCESS;
 
}
