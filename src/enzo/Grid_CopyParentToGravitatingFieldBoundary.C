/***********************************************************************
/
/  GRID CLASS (COPY PARENT DENSITY TO THE EXTRA BOUNDARY REGION OF
/              GRAVITATING MASS FIELD)
/
/  written by: Greg Bryan
/  date:       July, 1995
/  modified1:
/
/  PURPOSE:
/    The GravitatingMassField can have boundary points that will not be
/      set by the grid points on this level.  This routine sets those
/      points using a first-order interpolation from the parent grid.
/
/    Note: this must be the first routine called after clearing field.
/
/  RETURNS: FAIL or SUCCESS
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
#include "communication.h"
 
/* function prototypes */
 
extern "C" void FORTRAN_NAME(prolong)(float *source, float *dest, int *ndim,
				      int *sdim1, int *sdim2, int *sdim3,
				      int *ddim1, int *ddim2, int *ddim3,
				      int *start1, int *start2, int *start3,
				      int *refine1, int *refine2,int *refine3);
 
/* InterpolateBoundaryFromParent function */
 
int grid::CopyParentToGravitatingFieldBoundary(grid *ParentGrid)
{
  //  return SUCCESS;
  /* If this doesn't concern us, return. */
 
  if (this->CommunicationMethodShouldExit(ParentGrid))
    return SUCCESS;
 
  /* This routine is only required for subgrids. */
 
  if (GravityBoundaryType != SubGridIsolated)
    return SUCCESS;
 
  /* Declarations. */
 
  int ParentOffset[MAX_DIMENSION], ParentStartIndex[MAX_DIMENSION],
      ParentTempDim[MAX_DIMENSION], Refinement[MAX_DIMENSION],
      ParentDim[MAX_DIMENSION], SubGridExtra[MAX_DIMENSION],
      dim, i, j, k, gravityindex, size = 1;
 
  /* Compute refinement factors. */
 
  ParentGrid->ComputeRefinementFactors(this, Refinement);
 
  /* Set unused dims to zero. */
 
  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    ParentOffset[dim] = ParentStartIndex[dim] = SubGridExtra[dim] = 0;
    ParentTempDim[dim] = ParentDim[dim] = 1;
  }
 
  /* Compute the ParentOffset (in grid units) and ParentStartIndex and
     the region dim (in parent units). */
 
  for (dim = 0; dim < GridRank; dim++) {
    SubGridExtra[dim] = nint((GridLeftEdge[dim] -
			      GravitatingMassFieldLeftEdge[dim])
    			     /GravitatingMassFieldCellSize);
    // SubGridExtra[dim] = nint((CellLeftEdge[dim][0] -
    //  			      GravitatingMassFieldLeftEdge[dim])
    //  			     /GravitatingMassFieldCellSize);
    ParentOffset[dim] = nint((GravitatingMassFieldLeftEdge[dim] -
		  ParentGrid->GravitatingMassFieldLeftEdge[dim])/
			     GravitatingMassFieldCellSize);
    ParentStartIndex[dim] = ParentOffset[dim]/Refinement[dim] - 1;
    ParentTempDim[dim] = (ParentOffset[dim] +
		       GravitatingMassFieldDimension[dim]-1)/Refinement[dim] -
			 ParentStartIndex[dim] + 3;
    ParentDim[dim] = ParentGrid->GravitatingMassFieldDimension[dim];
    size *= GravitatingMassFieldDimension[dim];
    if (ParentStartIndex[dim] < 0 ||
	ParentStartIndex[dim]+ParentTempDim[dim] > ParentDim[dim]) {
      ENZO_VFAIL("ParentStartIndex[%"ISYM"] = %"ISYM" ParentTempDim = %"ISYM"(%"ISYM").\n",
	      dim, ParentStartIndex[dim], ParentTempDim[dim], ParentDim[dim])
    }
  }
 
  /* If posting a receive, then record details of call. */

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
    CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
    CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = ParentGrid;
    CommunicationReceiveCallType[CommunicationReceiveIndex] = 4;
  }
#endif /* USE_MPI */

  /* Copy data from other processor if needed (modify ParentDim and
     ParentStartIndex to reflect the fact that we are only coping part of
     the grid. */
 
  if (ProcessorNumber != ParentGrid->ProcessorNumber) {
    ParentGrid->CommunicationSendRegion(ParentGrid, ProcessorNumber,
	    GRAVITATING_MASS_FIELD, NEW_ONLY, ParentStartIndex, ParentTempDim);
    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	CommunicationDirection == COMMUNICATION_SEND)
      return SUCCESS;
    for (dim = 0; dim < GridRank; dim++) {
      ParentOffset[dim] -= Refinement[dim]*ParentStartIndex[dim];
      ParentDim[dim] = ParentTempDim[dim];
    }
  }
 
  /* Return if this is not our concern. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Interpolate (linear interpolation) */
 
#define NO_INTERPOLATE_LINEAR
 
#ifdef INTERPOLATE_LINEAR
 
  FORTRAN_NAME(prolong)(ParentGrid->GravitatingMassField,
			GravitatingMassField, &GridRank,
			ParentDim, ParentDim+1, ParentDim+2,
			GravitatingMassFieldDimension,
			GravitatingMassFieldDimension+1,
			GravitatingMassFieldDimension+2,
			ParentOffset, ParentOffset+1, ParentOffset+2,
			Refinement, Refinement+1, Refinement+2);
 
#else /* INTERPOLATE_LINEAR */
 
  /* Interpolate (nearest neighbour) */
 
  if(ParentGrid->GravitatingMassField == NULL) ENZO_FAIL("NO GMF in PARENT");
  int iparent, jparent, kparent, parentindex;
  for (k = 0; k < GravitatingMassFieldDimension[2]; k++) {
    kparent = nint((k + ParentOffset[2])/Refinement[2]);
    for (j = 0; j < GravitatingMassFieldDimension[1]; j++) {
      jparent = nint((j + ParentOffset[1])/Refinement[1]);
      parentindex = (kparent*ParentDim[1] + jparent)*ParentDim[0];
      gravityindex = (k*GravitatingMassFieldDimension[1] + j)*
	GravitatingMassFieldDimension[0];
      for (i = 0; i < GravitatingMassFieldDimension[0]; i++, gravityindex++) {
	iparent = nint((i+ParentOffset[0])/Refinement[0]);
	
	GravitatingMassField[gravityindex] =
	  ParentGrid->GravitatingMassField[parentindex+iparent];
	
      }
    }
  }
 
#endif /* INTERPOLATE_LINEAR */
 
  /* Clean up parent. */
 
  if (MyProcessorNumber != ParentGrid->ProcessorNumber) {
    delete [] ParentGrid->GravitatingMassField;
    ParentGrid->GravitatingMassField = NULL;
  }
 
  /* Add one to field to account for one subtracted in ComovingSourceTerm. */
 
  if (ComovingCoordinates)
    for (i = 0; i < size; i++)
      GravitatingMassField[i] += 1.0;
 
  /* Clear the region of GMF that will overlap with real grid points
     (i.e. clear the region that we shouldn't have set in the above loop). */
 
  for (k = SubGridExtra[2];
       k < GravitatingMassFieldDimension[2]-SubGridExtra[2]; k++)
    for (j = SubGridExtra[1];
	 j < GravitatingMassFieldDimension[1]-SubGridExtra[1]; j++) {
      gravityindex = (k*GravitatingMassFieldDimension[1] + j)*
	                GravitatingMassFieldDimension[0]
	           + SubGridExtra[0];
      //      if (j == GravitatingMassFieldDimension[1]/2 &&
      //	  k == GravitatingMassFieldDimension[2]/2)
      //	for (i = 0; i < GravitatingMassFieldDimension[0]; i++)
      //	  printf("%"ISYM" %"GSYM"\n", i,
      //		 GravitatingMassField[gravityindex-SubGridExtra[0]+i]);
      for (i = SubGridExtra[0];
	   i < GravitatingMassFieldDimension[0]-SubGridExtra[0];
	   i++, gravityindex++)
	GravitatingMassField[gravityindex] = 0;
    }
 
  return SUCCESS;
 
}
