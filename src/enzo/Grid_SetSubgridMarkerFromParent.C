/***********************************************************************
/
/  GRID CLASS (MARK SUBGRID)
/
/  written by: Tom Abel
/  date:       August 2004
/  modified1:  John Wise, May, 2010 -- modified the subgrid version to
/              accept parent grids on different processors.  Must use
/              grid IDs instead of grid pointers.
/
/  PURPOSE:
/
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
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

#define LEVEL_BIT_OFFSET 24

/* function prototypes */

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
			      int Target, int Tag, MPI_Comm CommWorld, 
			      int BufferSize);
#endif /* USE_MPI */


int grid::SetSubgridMarkerFromParent(grid *Parent, int level)
{

  /* Return if this grid is not on this processor. */

  if (this->CommunicationMethodShouldExit(Parent))
    return SUCCESS;

  /* declarations */
    
  int i, j, k, dim, field, index, pindex, size, activesize, gzsize;
  int pi, pj, pk, buffer_index, marker_level;
  int GZStart[MAX_DIMENSION], GZEnd[MAX_DIMENSION];
  FLOAT GridLeft[MAX_DIMENSION], GridRight[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  int Start[MAX_DIMENSION], End[MAX_DIMENSION];
  int ParentStart[MAX_DIMENSION], CellOffset[MAX_DIMENSION];
  int Refinement[MAX_DIMENSION], MarkerRefinement;
  Eint32 *buffer = NULL;
  grid *marker;

  /* check if the field has been allocated.  It should have been
     earlier! */ 

  if (MyProcessorNumber == ProcessorNumber)
    if (SubgridMarker == NULL)
      ENZO_VFAIL("SubgridMarker not allocated! level=%d", level)
  if (MyProcessorNumber == Parent->ProcessorNumber)
    if (Parent->SubgridMarker == NULL)
      ENZO_VFAIL("Parent SubgridMarker not allocated! level=%d", level)

  /* Calculate where the grid (including ghost zones) lies in the
     parent grid */

  Parent->ComputeRefinementFactors(this, Refinement);
  for (dim = 0; dim < GridRank; dim++) {
    ParentStart[dim] = 
      int((CellLeftEdge[dim][0] - Parent->CellLeftEdge[dim][0]) /
	  Parent->CellWidth[dim][0]);
    CellOffset[dim] = nint((CellLeftEdge[dim][0] -
			    Parent->CellLeftEdge[dim][ParentStart[dim]]) /
			   CellWidth[dim][0]);
  }
  for (dim = GridRank; dim < MAX_DIMENSION; dim++) {
    ParentStart[dim] = 1;
    CellOffset[dim] = 1;
  }

  /* If the parent grid is on the same processor, we can copy the
     SubgridMarker pointer directly.  Make sure not to overwrite any
     pre-existing pointers that were set from the siblings. */

  if (ProcessorNumber == Parent->ProcessorNumber) {
    
    for (dim = 0; dim < 3; dim++) {

      // Left face
      for (i = 0; i < 3; i++) {
	GZStart[i] = 0;
	GZEnd[i] = (i == dim) ? GridStartIndex[i]-1 : GridDimension[i]-1;
      }

      for (k = GZStart[2]; k <= GZEnd[2]; k++) {
	pk = ParentStart[2] + (k + CellOffset[2]) / Refinement[2];
	for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	  pj = ParentStart[1] + (j + CellOffset[1]) / Refinement[1];
	  index = (k*GridDimension[1] + j)*GridDimension[0] + GZStart[0];
	  for (i = GZStart[0]; i <= GZEnd[0]; i++, index++) {
	    pi = ParentStart[0] + (i + CellOffset[0]) / Refinement[0];
	    pindex = pi + 
	      Parent->GridDimension[0] * (pj + pk*Parent->GridDimension[1]);
	    SubgridMarker[index] = Parent->SubgridMarker[pindex];
	  }
	}
      }

      // Right face
      for (i = 0; i < 3; i++) {
	GZStart[i] = (i == dim) ? GridEndIndex[i]+1 : 0;
	GZEnd[i] = GridDimension[i]-1;
      }

      for (k = GZStart[2]; k <= GZEnd[2]; k++) {
	pk = ParentStart[2] + (k + CellOffset[2]) / Refinement[2];
	for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	  pj = ParentStart[1] + (j + CellOffset[1]) / Refinement[1];
	  index = (k*GridDimension[1] + j)*GridDimension[0] + GZStart[0];
	  for (i = GZStart[0]; i <= GZEnd[0]; i++, index++) {
	    pi = ParentStart[0] + (i + CellOffset[0]) / Refinement[0];
	    pindex = pi + 
	      Parent->GridDimension[0] * (pj + pk*Parent->GridDimension[1]);
	    SubgridMarker[index] = Parent->SubgridMarker[pindex];
	  }
	}
      }

    } // ENDFOR dim
    
  } // ENDIF same processor

#ifdef USE_MPI
  // Different processors
  else {

    // Calculate sizes of the ghost zones including overlap next to
    // the edges
    int dim1, dim2;
    for (dim = 0, gzsize = 0; dim < GridRank; dim++) {
      dim1 = (dim+1) % 3;
      dim2 = (dim+2) % 3;
      gzsize += 2 * GridStartIndex[dim] * GridDimension[dim1] * GridDimension[dim2];
    }

    if (CommunicationDirection == COMMUNICATION_RECEIVE)
      buffer = (Eint32*) CommunicationReceiveBuffer[CommunicationReceiveIndex];
    else
      buffer = new Eint32[gzsize];

    /* Copy grid IDs and levels to buffer on Parent processor.  IDs
       and levels are packed together in one integer with the level in
       the last 7 bits of the integer.  In future for simulations with
       more than 16M grids on a level, we should go to a 64-bit
       integer. */
    
    if (MyProcessorNumber == Parent->ProcessorNumber) {

      for (i = 0; i < gzsize; i++)
	buffer[i] = 0;

      buffer_index = 0;
      for (dim = 0; dim < 3; dim++) {

	// Left face
	for (i = 0; i < 3; i++) {
	  GZStart[i] = 0;
	  GZEnd[i] = (i == dim) ? GridStartIndex[i]-1 : GridDimension[i]-1;
	}

	for (k = GZStart[2]; k <= GZEnd[2]; k++) {
	  pk = ParentStart[2] + (k + CellOffset[2]) / Refinement[2];
	  for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	    pj = ParentStart[1] + (j + CellOffset[1]) / Refinement[1];
	    for (i = GZStart[0]; i <= GZEnd[0]; i++, buffer_index++) {
	      pi = ParentStart[0] + (i + CellOffset[0]) / Refinement[0];
	      pindex = pi + 
		Parent->GridDimension[0] * (pj + pk*Parent->GridDimension[1]);
	      marker = Parent->SubgridMarker[pindex];
	      if (marker != NULL) {
		buffer[buffer_index] = marker->GetGridID();

		// Determine the level of the grid pointed to by the
		// subgrid marker
		MarkerRefinement = 
		  nint( Parent->SubgridMarker[pindex]->CellWidth[0][0] /
			CellWidth[0][0] );
		marker_level = level;
		while (MarkerRefinement > 1) {
		  MarkerRefinement /= RefineBy;
		  marker_level--;
		}
		buffer[buffer_index] |= (marker_level << LEVEL_BIT_OFFSET);
	      } else {
		buffer[buffer_index] = INT_UNDEFINED;
	      }
	    }
	  }
	}

	// Right face
	for (i = 0; i < 3; i++) {
	  GZStart[i] = (i == dim) ? GridEndIndex[i]+1 : 0;
	  GZEnd[i] = GridDimension[i]-1;
	}

	for (k = GZStart[2]; k <= GZEnd[2]; k++) {
	  pk = ParentStart[2] + (k + CellOffset[2]) / Refinement[2];
	  for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	    pj = ParentStart[1] + (j + CellOffset[1]) / Refinement[1];
	    for (i = GZStart[0]; i <= GZEnd[0]; i++, buffer_index++) {
	      pi = ParentStart[0] + (i + CellOffset[0]) / Refinement[0];
	      pindex = pi + 
		Parent->GridDimension[0] * (pj + pk*Parent->GridDimension[1]);
	      marker = Parent->SubgridMarker[pindex];
	      if (marker != NULL) {
		buffer[buffer_index] = marker->GetGridID();

		// Determine the level of the grid pointed to by the
		// subgrid marker
		MarkerRefinement = 
		  nint( marker->CellWidth[0][0] / CellWidth[0][0] );
		marker_level = level;
		while (MarkerRefinement > 1) {
		  MarkerRefinement /= RefineBy;
		  marker_level--;
		}
		buffer[buffer_index] |= (marker_level << LEVEL_BIT_OFFSET);
	      } else {
		buffer[buffer_index] = INT_UNDEFINED;
	      }
	    }
	  }
	}

      } // ENDFOR dim

      /* Send the buffer */

//      printf("P%d: Send %d Subgrid markers from grid %d:%d to %d:%d\n",
//	     MyProcessorNumber, gzsize, level-1, Parent->ID, level, ID);
      CommunicationBufferedSend(buffer, gzsize, MPI_INT, ProcessorNumber,
				MPI_SENDMARKER_TAG, MPI_COMM_WORLD,
				BUFFER_IN_PLACE);

    } // ENDIF Parent processor

    else if (MyProcessorNumber == ProcessorNumber) {

      /* Post the receive call */

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
//	printf("P%d: Posting receive for %d Subgrid markers from "
//	       "grid %d:%d to %d:%d\n",
//	       MyProcessorNumber, gzsize, level-1, Parent->ID, level, ID);
	MPI_Irecv(buffer, gzsize, MPI_INT, Parent->ProcessorNumber,
		  MPI_SENDMARKER_TAG, MPI_COMM_WORLD,
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
	CommunicationReceiveBuffer[CommunicationReceiveIndex] = (float*) buffer;
	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = Parent;
	CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = level;
	CommunicationReceiveCallType[CommunicationReceiveIndex] = 19;
	CommunicationReceiveIndex++;
      } // ENDIF post receive

      /* Process the received data */

      if (CommunicationDirection == COMMUNICATION_RECEIVE) {

	/* Store the buffer in SubgridMarker, as we can't convert them
	   to a grid pointer here because we don't have LevelArray[].
	   Can't pass it without modifying
	   CommunicationReceiveHandler() everywhere */

//	printf("P%d: Processing receive for %d Subgrid markers from "
//	       "grid %d:%d to %d:%d\n",
//	       MyProcessorNumber, gzsize, level-1, Parent->ID, level, ID);

	buffer_index = 0;
	for (dim = 0; dim < 3; dim++) {

	  // Left face
	  for (i = 0; i < 3; i++) {
	    GZStart[i] = 0;
	    GZEnd[i] = (i == dim) ? GridStartIndex[i]-1 : GridDimension[i]-1;
	  }

	  for (k = GZStart[2]; k <= GZEnd[2]; k++)
	    for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	      index = GRIDINDEX_NOGHOST(GZStart[0], j, k);
	      for (i = GZStart[0]; i <= GZEnd[0]; i++, index++, buffer_index++)
		if (buffer[buffer_index] == INT_UNDEFINED)
		  SubgridMarker[index] = NULL;
		else
		  SubgridMarker[index] = (grid*) (long) buffer[buffer_index];
	    } // ENDFOR j

	  // Right face
	  for (i = 0; i < 3; i++) {
	    GZStart[i] = (i == dim) ? GridEndIndex[i]+1 : 0;
	    GZEnd[i] = GridDimension[i]-1;
	  }

	  for (k = GZStart[2]; k <= GZEnd[2]; k++)
	    for (j = GZStart[1]; j <= GZEnd[1]; j++) {
	      index = GRIDINDEX_NOGHOST(GZStart[0], j, k);
	      for (i = GZStart[0]; i <= GZEnd[0]; i++, index++, buffer_index++)
		if (buffer[buffer_index] == INT_UNDEFINED)
		  SubgridMarker[index] = NULL;
		else
		  SubgridMarker[index] = (grid*) (long) buffer[buffer_index];
	    } // ENDFOR j

	} // ENDFOR dim

	delete [] buffer;

      } // ENDIF receive

    } // ENDIF ProcessorNumber

  } // ENDELSE different processors
#endif /* USE_MPI */

  return SUCCESS;
  
}
