/***********************************************************************
/
/  GRID CLASS (SEND SUBGRID MARKER TO ANOTHER PROCESSOR)
/
/  written by: John Wise
/  date:       June, 2013
/  modified1:
/
/  NOTES: See SetSubgridMarker and the routines called from it for more
/         details on how the subgrid marker pointers are encoded and 
/         packed into an array of longs.
/
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <stdlib.h>
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
#include "CommunicationUtilities.h"

#define LEVEL_BIT_OFFSET 24
 
/* function prototypes */

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */

int grid::CommunicationSendSubgridMarker(grid *ToGrid, int ToProcessor)
{

#ifdef USE_MPI

  /* Note: No need to treat the case where ToProcessor ==
     ProcessorNumber.  This should never happen. */

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

  int i, dim, size, level, marker_level, MarkerRefinement;
  int Zero[] = {0, 0, 0};
  FLOAT FZero[] = {0.0, 0.0, 0.0};
  Eint32 *buffer;
  
  /* Calculate the grid size and level */

  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  level = nint(log(TopGridDx[0] / CellWidth[0][0]) / log(RefineBy));

  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = (Eint32*) CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
    buffer = new Eint32[size];

  /* If this is the from processor, pack the SubgridMarker */

  if (MyProcessorNumber == ProcessorNumber) {
    for (i = 0; i < size; i++) {
      if (SubgridMarker[i] != NULL) {

	buffer[i] = SubgridMarker[i]->GetGridID();
	if ( CellWidth[0][0] <= SubgridMarker[i]->CellWidth[0][0] ) {
	  // Coarser & sibling grids
	  MarkerRefinement = nint( SubgridMarker[i]->CellWidth[0][0] / 
				   CellWidth[0][0] );
	  marker_level = level;
	  while (MarkerRefinement > 1) {
	    MarkerRefinement /= RefineBy;
	    marker_level--;
	  }
	} else {
	  // Fine grids
	  MarkerRefinement = nint( CellWidth[0][0] / SubgridMarker[i]->CellWidth[0][0] );
	  marker_level = level;
	  while (MarkerRefinement > 1) {
	    MarkerRefinement /= RefineBy;
	    marker_level++;
	  }
	}

	buffer[i] |= (marker_level << LEVEL_BIT_OFFSET);

      } // ENDIF SubgridMarker != NULL
      else {
	buffer[i] = INT_UNDEFINED;
      }
    } // ENDFOR cells

    /* Send the buffer */

    CommunicationBufferedSend(buffer, size, MPI_INT, ToProcessor,
			      MPI_SENDMARKER_TAG, MPI_COMM_WORLD,
			      BUFFER_IN_PLACE);

  } // ENDIF from processor

  else if (MyProcessorNumber == ToProcessor) {

    /* Post the receive call */

    if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      MPI_Irecv(buffer, size, MPI_INT, ProcessorNumber,
		MPI_SENDMARKER_TAG, MPI_COMM_WORLD,
		CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
      CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
      CommunicationReceiveGridTwo[CommunicationReceiveIndex] = ToGrid;
      CommunicationReceiveCallType[CommunicationReceiveIndex] = 20;
      CommunicationReceiveBuffer[CommunicationReceiveIndex] = (float *) buffer;
      CommunicationReceiveDependsOn[CommunicationReceiveIndex] = 
	CommunicationReceiveCurrentDependsOn;
      CommunicationReceiveIndex++;
    }

    /* Process the received data.  This buffer is an encoded grid
       level/ID.  It will be decoded into a grid pointer after all of
       the communication is finished. */

    if (CommunicationDirection == COMMUNICATION_RECEIVE) {

      SubgridMarker = new grid*[size];
      for (i = 0; i < size; i++)
	SubgridMarker[i] = (grid*) buffer[i];
      
      delete [] buffer;

    }

  }

#endif /* USE_MPI */  

  return SUCCESS;
}
