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

  if (NumberOfStars == 0)
    return SUCCESS;

  /* Allocate buffer in ToProcessor.  This is automatically done in
     StarListToBuffer in the local processor. */

  TransferSize = this->NumberOfStars;
  if (MyProcessorNumber == ToProcessor)
    buffer = new StarBuffer[TransferSize];

  /* If this is the from processor, pack fields and delete stars. */

  if (MyProcessorNumber == ProcessorNumber) {
    buffer = this->Stars->StarListToBuffer(this->NumberOfStars);
    DeleteStarList(this->Stars);
  }
    
  /* Send buffer. */

#ifdef USE_MPI

  if (FirstTimeCalled) {
    MPI_Type_contiguous(sizeof(StarBuffer), MPI_BYTE, &MPI_STAR);
    MPI_Type_commit(&MPI_STAR);
    FirstTimeCalled = FALSE;
  }

  /* only send if processor numbers are not identical */

  if (ProcessorNumber != ToProcessor) {

    MPI_Status status;

#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif
    if (MyProcessorNumber == ProcessorNumber) {
      CommunicationBufferedSend(buffer, TransferSize, MPI_STAR, 
				ToProcessor, MPI_SENDSTAR_TAG, MPI_COMM_WORLD, 
				BUFFER_IN_PLACE);
    }

    if (MyProcessorNumber == ToProcessor) {
      if (MPI_Recv(buffer, TransferSize, MPI_STAR, ProcessorNumber,
		   MPI_SENDSTAR_TAG, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
	fprintf(stderr, "Proc %"ISYM" MPI_Recv error %"ISYM"\n", MyProcessorNumber,
		status.MPI_ERROR);
	fprintf(stderr, "P(%"ISYM"): TransferSize = %"ISYM" ProcessorNumber = %"ISYM"\n", 
		MyProcessorNumber, TransferSize, ProcessorNumber);
	char errstr[MPI_MAX_ERROR_STRING];
	Eint32 errlen;
	MPI_Error_string(status.MPI_ERROR, errstr, &errlen);
	fprintf(stderr, "MPI Error: %s\n", errstr);
	ENZO_FAIL("");
      }
    }

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

  if (MyProcessorNumber == ToProcessor) {

    RecvStars = StarBufferToList(buffer, TransferSize);
    InsertStarAfter(ToGrid->Stars, RecvStars);
			  
  } // end: if (MyProcessorNumber...)

  if (MyProcessorNumber == ToProcessor)
    delete [] buffer;

  return SUCCESS;
}

