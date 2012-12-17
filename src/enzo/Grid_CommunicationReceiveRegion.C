/***********************************************************************
/
/  GRID CLASS (RECEIVES FROM 'FAKE' GRID TO REAL GRID)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  Robert Harkness
/  date:       January, 2004
/
/  PURPOSE:
/
/  INPUTS:
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

// function prototypes
 
extern "C" void FORTRAN_NAME(copy3d)(float *source, float *dest,
                                   int *sdim1, int *sdim2, int *sdim3,
                                   int *ddim1, int *ddim2, int *ddim3,
                                   int *sstart1, int *sstart2, int *sstart3,
                                   int *dstart1, int *dstart2, int *dststart3);
 
#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */
 
 
int grid::CommunicationReceiveRegion(grid *FromGrid, int FromProcessor,
				     int SendField, int NewOrOld,
				     int RegionStart[], int RegionDim[],
				     int IncludeBoundary)
{
#ifdef USE_MPI 
  MPI_Request  RequestHandle;
  MPI_Status Status;
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
  MPI_Arg Source;
  MPI_Arg ZeroTag;

  int i, index, field, dim, Zero[] = {0, 0, 0};

  if (CommunicationShouldExit(FromProcessor, ProcessorNumber))
    return SUCCESS;

//  if (MyProcessorNumber != ProcessorNumber &&
//      MyProcessorNumber != FromProcessor)
//    return SUCCESS;

  // Compute size of region to transfer
  int NumberOfFields, RegionSize, TransferSize;
 
  NumberOfFields = ((SendField == ALL_FIELDS)? NumberOfBaryonFields : 1) *
    ((NewOrOld == NEW_AND_OLD)? 2 : 1);

  if (SendField == INTERPOLATED_FIELDS) {
    switch (OutputSmoothedDarkMatter) {
    case 1: NumberOfFields = 1; break;  // density
    case 2: NumberOfFields = 5; break;  // + rms velocity + 3-velocity
    }
  }

  RegionSize = RegionDim[0]*RegionDim[1]*RegionDim[2];
  TransferSize = RegionSize*NumberOfFields;
 
  // Allocate buffer
 
  float *buffer = NULL;
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
    buffer = new float[TransferSize];
 
  /* If this is the from processor, pack fields */
 
  int FromDim[MAX_DIMENSION], FromOffset[MAX_DIMENSION];
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    FromOffset[dim] = (dim < GridRank && IncludeBoundary == FALSE &&
		       SendField != INTERPOLATED_FIELDS) ?
      NumberOfGhostZones : 0;
    FromDim[dim] = RegionDim[dim] + 2*FromOffset[dim];
  }
 
  if (MyProcessorNumber == FromProcessor) {
 
    index = 0;
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < FromGrid->NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(FromGrid->BaryonField[field], &buffer[index],
			       FromDim, FromDim+1, FromDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       FromOffset, FromOffset+1, FromOffset+2);
	  index += RegionSize;
	}
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < FromGrid->NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  FORTRAN_NAME(copy3d)(FromGrid->OldBaryonField[field], &buffer[index],
			       FromDim, FromDim+1, FromDim+2,
			       RegionDim, RegionDim+1, RegionDim+2,
			       Zero, Zero+1, Zero+2,
			       FromOffset, FromOffset+1, FromOffset+2);
	  index += RegionSize;
      }

    if (SendField == INTERPOLATED_FIELDS) {
      for (field = 0; field < NumberOfFields; field++) {
	FORTRAN_NAME(copy3d)(FromGrid->InterpolatedField[field], &buffer[index],
			     FromDim, FromDim+1, FromDim+2,
			     RegionDim, RegionDim+1, RegionDim+2,
			     Zero, Zero+1, Zero+2,
			     FromOffset, FromOffset+1, FromOffset+2);
	index += RegionSize;
      }
    }
 
  } // ENDIF FromProcessor
 
  /* Send buffer */
 
  // Only send if processor numbers are not identical
 
  if (ProcessorNumber != FromProcessor) {
 
#ifdef MPI_INSTRUMENTATION
    starttime=MPI_Wtime();
#endif /* MPI_INSTRUMENTATION */
 
    if (MyProcessorNumber == FromProcessor) {
#ifdef MPI_INSTRUMENTATION
      if (traceMPI) 
	fprintf(tracePtr, "CRR RF: Sending %"ISYM" floats from %"ISYM" to %"ISYM"\n", 
		TransferSize, FromProcessor, ProcessorNumber);
#endif
      CommunicationBufferedSend(buffer, TransferSize, DataType, ProcessorNumber, 
				0, MPI_COMM_WORLD, BUFFER_IN_PLACE);
    } // ENDIF from processor
    
    if (MyProcessorNumber == ProcessorNumber) {

      /* Post the receive message without waiting for the message to
	 be received.  When the data arrives, this will be called again
	 in (the real) receive mode. */

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	MPI_Irecv(buffer, TransferSize, DataType, FromProcessor, 0, 
		  MPI_COMM_WORLD, 
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);
	CommunicationReceiveBuffer[CommunicationReceiveIndex] = buffer;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] =
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;      }

      /* If in send-receive mode, then wait for the message now. */

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE) {
	MPI_Recv(buffer, TransferSize, DataType, FromProcessor, 0, 
		 MPI_COMM_WORLD, &Status);
      }

    } // ENDIF grid processor
 
#ifdef MPI_INSTRUMENTATION
    endtime=MPI_Wtime();
    timer[5] += endtime-starttime;
    counter[5] ++;
    timer[6] += double(TransferSize);
    RecvComm += endtime-starttime;
    CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
  } // ENDIF different processors
 
  /* If this is the to processor, unpack fields */
 
  int GridSize = GridDimension[0]*GridDimension[1]*GridDimension[2];
  int ActiveDims[MAX_DIMENSION], ActiveSize = 1;
  for (dim = 0; dim < GridRank; dim++) {
    ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    ActiveSize *= ActiveDims[dim];
  }
 
  if (MyProcessorNumber == ProcessorNumber &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {
 
    index = 0;
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == NEW_ONLY)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  if (BaryonField[field] == NULL) {
	    BaryonField[field] = new float[GridSize];
	    for (i = 0; i < GridSize; i++)
	      BaryonField[field][i] = 0;
          }
	  FORTRAN_NAME(copy3d)(&buffer[index], BaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionStart, RegionStart+1, RegionStart+2,
			       Zero, Zero+1, Zero+2);
 
	  index += RegionSize;
	}
 
    if (NewOrOld == NEW_AND_OLD || NewOrOld == OLD_ONLY)
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (field == SendField || SendField == ALL_FIELDS) {
	  if (OldBaryonField[field] == NULL) {
	    OldBaryonField[field] = new float[GridSize];
	    for (i = 0; i < GridSize; i++)
	      BaryonField[field][i] = 0;
          }
	  FORTRAN_NAME(copy3d)(&buffer[index], OldBaryonField[field],
			       RegionDim, RegionDim+1, RegionDim+2,
			       GridDimension, GridDimension+1, GridDimension+2,
			       RegionStart, RegionStart+1, RegionStart+2,
			       Zero, Zero+1, Zero+2);
 
	  index += RegionSize;
	}

    if (SendField == INTERPOLATED_FIELDS)
      for (field = 0; field < NumberOfFields; field++) {
	if (InterpolatedField[field] == NULL) {
	  InterpolatedField[field] = new float[ActiveSize];
	  for (i = 0; i < ActiveSize; i++)
	    InterpolatedField[field][i] = 0;
	}
	FORTRAN_NAME(copy3d)(&buffer[index], InterpolatedField[field],
			     RegionDim, RegionDim+1, RegionDim+2,
			     ActiveDims, ActiveDims+1, ActiveDims+2,
			     RegionStart, RegionStart+1, RegionStart+2,
			     Zero, Zero+1, Zero+2);
 
	index += RegionSize;
      }

 
    /* Clean up */
 
    delete [] buffer;

  } // ENDIF unpack
 
#endif /* USE_MPI */ 

  return SUCCESS;
}
