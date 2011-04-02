/***********************************************************************
/
/  GRID CLASS (COPY OVERLAPPING OVERLAPPING POTENTIAL FIELD)
/
/  written by: Greg Bryan
/  date:       January, 1998
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
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "communication.h"
 
int CopyPotentialFieldAverage = 0;
 
int grid::CopyPotentialField(grid *OtherGrid, FLOAT EdgeOffset[MAX_DIMENSION])
{
  /* Return if this doesn't involve us. */
 
  if (ProcessorNumber != MyProcessorNumber &&
      OtherGrid->ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  int dim, OnlyBoundary = FALSE;
 
  /* Compute the left and right edges of this grid (including ghost zones). */
#ifdef UNUSED
  FLOAT AnyOffset = EdgeOffset[0] + EdgeOffset[1] + EdgeOffset[2];
#endif /* UNUSED */
  FLOAT GridLeft[MAX_DIMENSION], GridRight[MAX_DIMENSION];
  FLOAT OtherGridLeft[MAX_DIMENSION], OtherGridRight[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    GridLeft[dim]  = GravitatingMassFieldLeftEdge[dim] + EdgeOffset[dim];
    GridRight[dim] = GridLeft[dim] + GravitatingMassFieldDimension[dim]*
                                     GravitatingMassFieldCellSize;
    if (CopyPotentialFieldAverage > 0) {
      OnlyBoundary = TRUE;
      OtherGridLeft[dim] = OtherGrid->GravitatingMassFieldLeftEdge[dim]
	 + OtherGrid->GravitatingMassFieldCellSize;
      OtherGridRight[dim] = OtherGrid->GravitatingMassFieldLeftEdge[dim] +
		(OtherGrid->GravitatingMassFieldDimension[dim] - 1)*
	//	(OtherGrid->GravitatingMassFieldDimension[dim])*
	OtherGrid->GravitatingMassFieldCellSize;
      //      OtherGridLeft[dim] = OtherGrid->CellLeftEdge[dim][0];
      //      OtherGridRight[dim] = OtherGrid->CellLeftEdge[dim][0] +
      //	OtherGrid->GridDimension[dim]*OtherGrid->CellWidth[dim][0];
    } else {
      OtherGridLeft[dim] = OtherGrid->GridLeftEdge[dim];
      OtherGridRight[dim] = OtherGrid->GridRightEdge[dim];
    }
 
    /* This little bit of logic makes sure that the outside edges of
       the cube get copied by modifying the OtherGrid edges.  Only do
       it if this is a periodic boundary copy in another dimension. */
 
#ifdef UNUSED
    if (fabs(AnyOffset) > 0 && EdgeOffset[dim] == 0) {
      if (nint((OtherGrid->GridLeftEdge[dim]-DomainLeftEdge[dim])/
	       CellWidth[dim][0]) == 0)
	OtherGridLeft[dim] = OtherGrid->CellLeftEdge[dim][0];
      FLOAT OtherGridRE = OtherGrid->CellLeftEdge[dim][0] +
	OtherGrid->GridDimension[dim]*OtherGrid->CellWidth[dim][0];
      if (nint((OtherGridRE - DomainRightEdge[dim])/CellWidth[dim][0]) == 0)
	OtherGridRight[dim] = OtherGridRE;
    }
#endif /* UNUSED */
  }
 
  /* Do a quick check to see if there is any overlap. */
 
  for (dim = 0; dim < GridRank; dim++)
    if (GridLeft[dim]  >= OtherGridRight[dim] ||
        GridRight[dim] <= OtherGridLeft[dim]   )
      return SUCCESS;
 
  /* There is some overlap, so copy overlapping region */
 
  FLOAT Left, Right;
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
  int StartOther[MAX_DIMENSION], Dim[MAX_DIMENSION];
  int OtherDim[MAX_DIMENSION];
 
  /* compute start and stop indicies of overlapping region for both this
     grid and the Other grid. */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Start[dim] = End[dim] = StartOther[dim] = 0;
    OtherDim[dim]   = 1;
  }
 
  for (dim = 0; dim < GridRank; dim++)
    if (GridDimension[dim] > 1) {
 
      /* Compute left and right positions in problem space.
	 note: include buffer zones of this grid but not the other grid. */
 
      Left  = max(GridLeft[dim], OtherGridLeft[dim]);
      Right = min(GridRight[dim], OtherGridRight[dim]);
 
      /* Convert this to index positions in this grid */
 
      Start[dim] = nint((Left  - GridLeft[dim]) / CellWidth[dim][0]);
      End[dim]   = nint((Right - GridLeft[dim]) / CellWidth[dim][0]) - 1;
 
      if (End[dim] - Start[dim] < 0)
	return SUCCESS;
 
      /* Compute index positions in the other grid */
 
      StartOther[dim] = nint((Left -
			      OtherGrid->GravitatingMassFieldLeftEdge[dim])/
			     CellWidth[dim][0]);
 
      /* Copy dimensions into temporary space */
 
      OtherDim[dim] = OtherGrid->GravitatingMassFieldDimension[dim];
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
    CommunicationReceiveCallType[CommunicationReceiveIndex] = 9;
    for (dim = 0; dim < GridRank; dim++)
      CommunicationReceiveArgument[dim][CommunicationReceiveIndex] = 
	EdgeOffset[dim];
  }
#endif /* USE_MPI */

  /* Copy data from other processor if needed (modify OtherDim and
     StartOther to reflect the fact that we are only coping part of
     the grid. */
 
  if (ProcessorNumber != OtherGrid->ProcessorNumber) {
    OtherGrid->CommunicationSendRegion(OtherGrid, ProcessorNumber,
				 POTENTIAL_FIELD, NEW_ONLY, StartOther, Dim);
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
 
  int i, j, k, thisindex, otherindex, OnBoundary;
  for (k = 0; k < Dim[2]; k++)
    for (j = 0; j < Dim[1]; j++) {
      thisindex = ((k + Start[2])*GravitatingMassFieldDimension[1] +
		   (j + Start[1]))*GravitatingMassFieldDimension[0] +
	                Start[0];
      otherindex = (0 + StartOther[0]) + (j + StartOther[1])*OtherDim[0] +
	           (k + StartOther[2])*OtherDim[0]*OtherDim[1];
 
      /* Determine if we are on the j or k boundary. */
 
      OnBoundary = FALSE;
      if (GridRank > 1 && (j+Start[1] == 0 ||
			   j+Start[1] == GravitatingMassFieldDimension[1]-1))
	OnBoundary = TRUE;
      if (GridRank > 2 && (k+Start[2] == 0 ||
			   k+Start[2] == GravitatingMassFieldDimension[2]-1))
	OnBoundary = TRUE;
 
      if (OnlyBoundary == TRUE && OnBoundary == FALSE) {
 
	/* Only copy the endpoints. */
 
	if (CopyPotentialFieldAverage == 2) {
	  if (Start[0] == 0)
	    PotentialField[thisindex] = 0.5*(PotentialField[thisindex] +
		  		   OtherGrid->PotentialField[otherindex]);
	  if (Start[0]+Dim[0] == GravitatingMassFieldDimension[0])
	    PotentialField[thisindex+Dim[0]-1] =
	      0.5*(PotentialField[thisindex+Dim[0]-1] +
		   OtherGrid->PotentialField[otherindex+Dim[0]-1]);
	} else {
	  if (Start[0] == 0)
	    PotentialField[thisindex] = OtherGrid->PotentialField[otherindex];
	  if (Start[0]+Dim[0] == GravitatingMassFieldDimension[0])
	    PotentialField[thisindex+Dim[0]-1] =
	      OtherGrid->PotentialField[otherindex+Dim[0]-1];
	}
 
      } else {
	
	/* Copy the whole line. */
 
	if (CopyPotentialFieldAverage == 2)
	  for (i = 0; i < Dim[0]; i++, thisindex++, otherindex++)
	    PotentialField[thisindex] = 0.5*(PotentialField[thisindex] +
		  		   OtherGrid->PotentialField[otherindex]);
	else
	  for (i = 0; i < Dim[0]; i++, thisindex++, otherindex++)
	    PotentialField[thisindex] = OtherGrid->PotentialField[otherindex];
 
      }
    }
 
  /* Clean up if we have transfered data. */
 
  if (MyProcessorNumber != OtherGrid->ProcessorNumber)
    OtherGrid->DeleteAllFields();
 
  return SUCCESS;
}
