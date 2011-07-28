/***********************************************************************
/
/  GRID CLASS (ADD OVERLAPPING PARTICLEMASSFIELD TO TARGET MASSFIELD)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
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
 
int grid::AddOverlappingParticleMassField(grid *OtherGrid,
					  FLOAT EdgeOffset[MAX_DIMENSION])
{
 
  if (MyProcessorNumber != ProcessorNumber &&
      MyProcessorNumber != OtherGrid->ProcessorNumber)
    return SUCCESS;

  //  fprintf(stderr, "AddOverlappiongMassField: %i %i", this, OtherGrid);

  /* declarations */

  int i, j, k, dim, thisindex, otherindex;
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
 
  /* Do a quick check to see if there is any overlap. */
 
  for (dim = 0; dim < GridRank; dim++) {
 
    Left[dim] = max(
          GravitatingMassFieldLeftEdge[dim] + EdgeOffset[dim],
	  OtherGrid->GravitatingMassFieldParticlesLeftEdge[dim]);
 
    Right[dim] = min(
	  //	  GridRightEdge[dim] +
	  //	  GridLeftEdge[dim] - GravitatingMassFieldLeftEdge[dim] +
	  GravitatingMassFieldLeftEdge[dim] +
	  GravitatingMassFieldCellSize * GravitatingMassFieldDimension[dim] +
	  EdgeOffset[dim],
	  OtherGrid->GravitatingMassFieldParticlesLeftEdge[dim] +
	  OtherGrid->GravitatingMassFieldParticlesCellSize *
	  OtherGrid->GravitatingMassFieldParticlesDimension[dim]);
 
    if (Left[dim]+0.5*GravitatingMassFieldCellSize >= Right[dim])
      return SUCCESS;
  }
 
  /* There is some overlap, so copy overlapping region */
 
  int Start[MAX_DIMENSION], End[MAX_DIMENSION], StartOther[MAX_DIMENSION];
  int ThisDim[MAX_DIMENSION], OtherDim[MAX_DIMENSION], Dim[MAX_DIMENSION];
 
  /* compute start and stop indicies of overlapping region for both this
     grid and the Other grid. */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Start[dim] = End[dim] = StartOther[dim] = 0;
    ThisDim[dim] = OtherDim[dim] = 1;
  }
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* Convert Left/Right to index positions in this grid. */
 
    Start[dim] = nint((Left[dim] -
			(GravitatingMassFieldLeftEdge[dim] + EdgeOffset[dim]))/
			        GravitatingMassFieldCellSize);
    End[dim] = nint((Right[dim]  -
		       (GravitatingMassFieldLeftEdge[dim] + EdgeOffset[dim]))/
		                GravitatingMassFieldCellSize) - 1;
 
    /* Compute index positions in the other grid. */
 
    StartOther[dim] = nint((Left[dim] -
			OtherGrid->GravitatingMassFieldParticlesLeftEdge[dim])/
		        OtherGrid->GravitatingMassFieldParticlesCellSize);
 
    /* copy dimensions into temporary space */
 
    ThisDim[dim] = GravitatingMassFieldDimension[dim];
    OtherDim[dim] = OtherGrid->GravitatingMassFieldParticlesDimension[dim];
  }
 
  /* Calculate dimensions */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Dim[dim] = End[dim] - Start[dim] + 1;

  /* If posting a receive, then record details of call. */

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_POST_RECEIVE &&
      MyProcessorNumber == ProcessorNumber) {
    CommunicationReceiveGridOne[CommunicationReceiveIndex]  = this;
    CommunicationReceiveGridTwo[CommunicationReceiveIndex]  = OtherGrid;
    CommunicationReceiveCallType[CommunicationReceiveIndex] = 6;
    for (dim = 0; dim < GridRank; dim++)
      CommunicationReceiveArgument[dim][CommunicationReceiveIndex] = 
	EdgeOffset[dim];
  }
#endif /* USE_MPI */
 
  /* Copy data from other processor if needed (modify ParentDim and
     ParentStartIndex to reflect the fact that we are only coping part of
     the grid. */
 
  if (ProcessorNumber != OtherGrid->ProcessorNumber) {
    OtherGrid->CommunicationSendRegion(OtherGrid, ProcessorNumber, 
	                GRAVITATING_MASS_FIELD_PARTICLES, NEW_ONLY, 
		        StartOther, Dim);
    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	CommunicationDirection == COMMUNICATION_SEND)
      return SUCCESS;    
    for (dim = 0; dim < GridRank; dim++) {
      OtherDim[dim] = Dim[dim];
      StartOther[dim] = 0;
    }
  }
   
  /* Return if this is not our concern. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Add zones */
 
  for (k = 0; k < Dim[2]; k++)
    for (j = 0; j < Dim[1]; j++) {
      thisindex = ((k + Start[2])*ThisDim[1] +
                   (j + Start[1])              )*ThisDim[0] +
	           (0 + Start[0]);
      otherindex = ((k + StartOther[2])*OtherDim[1] +
		    (j + StartOther[1])        )*OtherDim[0] +
		    (0 + StartOther[0]);
      for (i = 0; i < Dim[0]; i++, thisindex++, otherindex++)
	GravitatingMassField[thisindex] +=
	  OtherGrid->GravitatingMassFieldParticles[otherindex];
    }
 
  /* Clean up if we have transfered data. */
 
  if (MyProcessorNumber != OtherGrid->ProcessorNumber)
    OtherGrid->DeleteAllFields();
 
  return SUCCESS;
}
