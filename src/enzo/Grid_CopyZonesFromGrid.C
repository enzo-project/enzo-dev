/***********************************************************************
/
/  GRID CLASS (COPY OVERLAPPING ZONES FROM GRID IN ARGUMENT TO THIS GRID)
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
 
// This routine copies zones which overlap from the grid in the argument
//   to the current grid.  We use only the active region of the OtherGrid,
//   but copy into the entire region (including boundaries) of this grid.
//
// The input argument EdgeOffset is the amount the corner of this grid is
//   considered to have moved for grid comparison and copying purposes.
//   See Grid_CheckForOverlappingZones for more details.

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "communication.h"
 
int grid::CopyZonesFromGrid(grid *OtherGrid, FLOAT EdgeOffset[MAX_DIMENSION])
{
 
  /* Return if this doesn't involve us. */
 
  if (ProcessorNumber != MyProcessorNumber &&
      OtherGrid->ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
//  printf("CopyZonesFromGrid: %"ISYM"\n", NumberOfBaryonFields);
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  this->DebugCheck("CopyZonesFromGrid (before)");
 
  /* declarations */
 
  int dim;
 
  /* Compute the left and right edges of this grid (including ghost zones). */
 
  FLOAT GridLeft[MAX_DIMENSION], GridRight[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    GridLeft[dim]  = CellLeftEdge[dim][0] + EdgeOffset[dim];
    GridRight[dim] = CellLeftEdge[dim][GridDimension[dim]-1] +
                     CellWidth[dim][GridDimension[dim]-1]    +
                     EdgeOffset[dim];
  }
 
  /* Do a quick check to see if there is any overlap. */
 
  for (dim = 0; dim < GridRank; dim++)
    if (GridLeft[dim]  >= OtherGrid->GridRightEdge[dim] ||
        GridRight[dim] <= OtherGrid->GridLeftEdge[dim]   )
      return SUCCESS;
 
  /* There is some overlap, so copy overlapping region */
 
  FLOAT Left, Right;
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
  int StartOther[MAX_DIMENSION], Dim[MAX_DIMENSION];
  int OtherDim[MAX_DIMENSION];
 
  /* compute start and stop indicies of overlapping region for both this
     grid and the Other grid. */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Start[dim]      = 0;
    End[dim]        = 0;
    StartOther[dim] = 0;
    OtherDim[dim]   = 1;
  }
 
  for (dim = 0; dim < GridRank; dim++)
    if (GridDimension[dim] > 1) {
 
      /* Compute left and right positions in problem space.
	 note: include buffer zones of this grid but not the other grid. */
 
      Left  = max(GridLeft[dim], OtherGrid->GridLeftEdge[dim]);
      Right = min(GridRight[dim], OtherGrid->GridRightEdge[dim]);
 
      /* Convert this to index positions in this grid */
 
      Start[dim] = nint((Left  - GridLeft[dim]) / CellWidth[dim][0]);
      End[dim]   = nint((Right - GridLeft[dim]) / CellWidth[dim][0]) - 1;
 
      if (End[dim] - Start[dim] < 0)
	return SUCCESS;
 
      /* Compute index positions in the other grid */
 
      StartOther[dim] = nint((Left - OtherGrid->CellLeftEdge[dim][0])/
			     CellWidth[dim][0]);
 
      /* Copy dimensions into temporary space */
 
      OtherDim[dim] = OtherGrid->GridDimension[dim];
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
    CommunicationReceiveCallType[CommunicationReceiveIndex] = 2;
    for (dim = 0; dim < GridRank; dim++)
      CommunicationReceiveArgument[dim][CommunicationReceiveIndex] = 
	EdgeOffset[dim];
  }
#endif /* USE_MPI */

  /* Copy data from other processor if needed (modify OtherDim and
     StartOther to reflect the fact that we are only coping part of
     the grid. */
 
  if (traceMPI) 
    fprintf(tracePtr, "CopyZones SendRegion from %"ISYM" to %"ISYM"\n", 
	    ProcessorNumber, OtherGrid->ProcessorNumber);
 
  if (ProcessorNumber != OtherGrid->ProcessorNumber) {
    OtherGrid->CommunicationSendRegion(OtherGrid, ProcessorNumber,
				       ALL_FIELDS, NEW_ONLY, StartOther, Dim);
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
 
  /* Copy zones */
 
  int thisindex, otherindex;
  for (int field = 0; field < NumberOfBaryonFields; field++)
    for (int k = 0; k < Dim[2]; k++)
      for (int j = 0; j < Dim[1]; j++) {
	thisindex = (0 + Start[0]) + (j + Start[1])*GridDimension[0] +
                    (k + Start[2])*GridDimension[0]*GridDimension[1];
	otherindex = (0 + StartOther[0]) + (j + StartOther[1])*OtherDim[0] +
                     (k + StartOther[2])*OtherDim[0]*OtherDim[1];
	for (int i = 0; i < Dim[0]; i++, thisindex++, otherindex++)
	  BaryonField[field][thisindex] =
	    OtherGrid->BaryonField[field][otherindex];
      }
 
  /* Clean up if we have transfered data. */
 
  if (MyProcessorNumber != OtherGrid->ProcessorNumber)
    OtherGrid->DeleteAllFields();
 
  this->DebugCheck("CopyZonesFromGrid (after)");
 
  return SUCCESS;
 
}
