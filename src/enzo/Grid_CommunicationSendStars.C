/***********************************************************************
/
/  GRID CLASS (SEND STARS FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  NOTES:  Adapted from grid::CommunicationSendParticles().
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
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

/* function prototypes */

#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_STAR;
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
                              int Target, int Tag, MPI_Comm CommWorld, 
			      int BufferSize);
#endif /* USE_MPI */
Star* StarBufferToList(StarBuffer *buffer, int n);
void InsertStarAfter(Star * &Node, Star * &NewNode);
void DeleteStarList(Star * &Node);

/* Send particle from this grid to ToGrid on processor ToProcessor, using
   FromNumber particles counting from FromStart.  Place into ToGrid at
   particle number ToStart. If ToStart = -1, then add to end. */

int grid::CommunicationSendStars(grid *ToGrid, int ToProcessor)
{

  int i, j, dim, index, TransferSize;
  StarBuffer *buffer = NULL;
  Star *RecvStars = NULL;
  Star *cstar;

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

  if (NumberOfStars == 0)
    return SUCCESS;

  /* Allocate buffer in ToProcessor.  This is automatically done in
     StarListToBuffer in the local processor. */

  TransferSize = this->NumberOfStars;
#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = (StarBuffer*) CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
#endif
    buffer = new StarBuffer[TransferSize];

  /* If this is the from processor, pack fields and delete stars. */

  if (MyProcessorNumber == ProcessorNumber) {
    buffer = this->Stars->StarListToBuffer(this->NumberOfStars);
    DeleteStarList(this->Stars);
  }
    
  /* Send buffer. */

#ifdef USE_MPI

  /* only send if processor numbers are not identical */

  if (ProcessorNumber != ToProcessor) {

    if (FirstTimeCalled) {
      MPI_Type_contiguous(sizeof(StarBuffer), MPI_BYTE, &MPI_STAR);
      MPI_Type_commit(&MPI_STAR);
      FirstTimeCalled = FALSE;
    }

    MPI_Status status;
    MPI_Arg PCount, Count = TransferSize;
    MPI_Arg Source = ProcessorNumber;
    MPI_Arg Dest = ToProcessor;
    MPI_Arg stat;

#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif
    if (MyProcessorNumber == ProcessorNumber)
      CommunicationBufferedSend(buffer, Count, MPI_STAR, 
				Dest, MPI_SENDSTAR_TAG, MPI_COMM_WORLD, 
				BUFFER_IN_PLACE);

    if (MyProcessorNumber == ToProcessor) {

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	MPI_Irecv(buffer, Count, MPI_STAR, Source,
		  MPI_SENDSTAR_TAG, MPI_COMM_WORLD,
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);

	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = ToGrid;
	CommunicationReceiveCallType[CommunicationReceiveIndex] = 18;
	CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = NumberOfStars;

	CommunicationReceiveBuffer[CommunicationReceiveIndex] = (float *) buffer;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] = 
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;
      }

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
	MPI_Recv(buffer, Count, MPI_STAR, Source,
		 MPI_SENDSTAR_TAG, MPI_COMM_WORLD, &status);

    } // ENDIF MyProcessorNumber == ToProcessor

#ifdef MPI_INSTRUMENTATION
    /* Zhiling Lan's instrumented part */
    endtime = MPI_Wtime();
    timer[7] += endtime-starttime;
    counter[7] ++;
    timer[8] += double(TransferSize);
    timer[28] += double(TransferSize*TransferSize);
    timer[27] += (endtime-starttime)*(endtime-starttime);
#endif /* MPI_INSTRUMENTATION */
  
  } // end: if (ProcessorNumber != ToProcessor)

#endif /* USE_MPI */

  /* If this is the to processor, unpack fields. */

  if (MyProcessorNumber == ToProcessor &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {

    RecvStars = StarBufferToList(buffer, TransferSize);
    InsertStarAfter(ToGrid->Stars, RecvStars);
    for (cstar = ToGrid->Stars; cstar; cstar = cstar->NextStar)
      cstar->CurrentGrid = ToGrid;
    delete [] buffer;
			  
  } // end: if (MyProcessorNumber...)

  return SUCCESS;
}

