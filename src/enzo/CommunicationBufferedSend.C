/***********************************************************************
/
/  COMMUNICATION ROUTINE: BUFFERED SEND WITH SELF-BUFFERING
/
/  written by: Greg Bryan
/  date:       January, 2001
/  modified1:
/
/  PURPOSE:
/    A replacement for MPI_Bsend, this routine allocates a buffer if
/     the BufferSize is non-negative and copies the data into that
/     buffer (otherwise the inbuffer is used directly).
/
************************************************************************/
 
#ifdef USE_MPI
 
#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h" 
void my_exit(int status);
/* Records the number of times we've been called. */
 
static int CallCount = 0;
 
/* Defines the number of calls to wait before scanning. */
 
#define NUMBER_OF_CALLS_BETWEEN_SCANS 30
 
/* The MPI Handle and buffer storage area. */
 
#define MAX_NUMBER_OF_MPI_BUFFERS 40000
 
static MPI_Request  RequestHandle[MAX_NUMBER_OF_MPI_BUFFERS];
static char        *RequestBuffer[MAX_NUMBER_OF_MPI_BUFFERS];
static int          LastActiveIndex = -1;
 
 
/* function prototypes */

int CommunicationBufferPurge(void) { 

  //fprintf(stderr,"CCO p%"ISYM" LastActive %"ISYM"!\n", MyProcessorNumber, LastActiveIndex);

  int i;
  MPI_Arg RequestDone;
  MPI_Arg stat;
  MPI_Status Status;
  
  int NewLastActiveIndex = -1;

  int BuffersPurged = 0;
  int BuffersActive = 0;

  for (i = 0; i < LastActiveIndex+1; i++) {
    if (RequestBuffer[i] != NULL) {
      stat = MPI_Test(RequestHandle+i, &RequestDone, &Status);
      if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
      if (RequestDone) {
	
	/* If the request is done, deallocate associated buffer. */

	//fprintf(stderr,"CCO p%"ISYM": mem- thread %"ISYM" finished\n",MyProcessorNumber, i);
	
	delete [] RequestBuffer[i];
	RequestBuffer[i] = NULL;
        BuffersPurged++;
        //fprintf(stderr, "CBP buffer %"ISYM" released\n", i);
	
      } else{

	//fprintf(stderr,"CCO p%"ISYM": mem- thread %"ISYM" active\n",MyProcessorNumber, i);

	NewLastActiveIndex = max(i, NewLastActiveIndex);
        BuffersActive++;
        //fprintf(stderr, "CBP buffer %"ISYM" remains active\n", i);

      }
    }
  } // end: loop over request handles

  LastActiveIndex = NewLastActiveIndex;

  // if (BuffersPurged != 0)
  // fprintf(stderr, "CBP %"ISYM": %"ISYM" buffers purged, %"ISYM" buffers remain active\n",
  //                  MyProcessorNumber, BuffersPurged, BuffersActive);

  return SUCCESS;
}

int CommunicationBufferedSendCancel(int Tag)
{
  
  /* Cancels all buffered sends with Tag */

  int i;
  MPI_Arg RequestDone, stat;
  MPI_Status Status;
  int NewLastActiveIndex = -1;
  int BuffersCancelled = 0;
  int BuffersActive = 0;

  for (i = 0; i < LastActiveIndex+1; i++) {
    if (RequestBuffer[i] != NULL) {
      stat = MPI_Test(RequestHandle+i, &RequestDone, &Status);
      if (stat != MPI_SUCCESS)
	ENZO_FAIL("Error in MPI_Test");
      if (Status.MPI_TAG == Tag && !RequestDone) {
	MPI_Cancel(RequestHandle+i);
	MPI_Wait(RequestHandle+i, MPI_STATUS_IGNORE);
	delete [] RequestBuffer[i];
	RequestBuffer[i] = NULL;
	BuffersCancelled++;
      } // ENDIF matching tag
      else {
	NewLastActiveIndex = max(i, NewLastActiveIndex);
	BuffersActive++;
      }
    } // ENDIF RequestBuffer[i] != NULL
  } // ENDFOR requests

  LastActiveIndex = NewLastActiveIndex;

  return SUCCESS;

}


int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize)
{
 
  int i;
  MPI_Arg RequestDone;
  MPI_Arg stat;
  MPI_Status Status;
  void *buffer_send;
 
  /* First, check to see if we should do a scan. */
 
  if (++CallCount % NUMBER_OF_CALLS_BETWEEN_SCANS == 0) {
 
    int NewLastActiveIndex = -1;
    for (i = 0; i < LastActiveIndex+1; i++) {
      if (RequestBuffer[i] != NULL) {
	stat = MPI_Test(RequestHandle+i, &RequestDone, &Status);
          if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
	if (RequestDone) {
 
	  /* If the request is done, deallocate associated buffer. */
 
	  delete [] RequestBuffer[i];
	  RequestBuffer[i] = NULL;
 
	} else
	  NewLastActiveIndex = max(i, NewLastActiveIndex);
      }
    } // end: loop over request handles
 
  }
 
  /* If necessary, allocate buffer. */
 
  if (BufferSize != BUFFER_IN_PLACE) {
    buffer_send = new char[BufferSize];
    memcpy(buffer_send, buffer, BufferSize);
  }
  else
    buffer_send = (void *) buffer;
 
  /* Find open spot. */
 
  int index = LastActiveIndex+1;
  for (i = 0; i < LastActiveIndex+1; i++)
    if (RequestBuffer[i] == NULL)
      index = i;
 
  /* Error check. */
 
  if (index >= MAX_NUMBER_OF_MPI_BUFFERS-1) {
    fprintf(stderr, "CommunicationBufferedSend: increase MAX_NUMBER_OF_MPI_BUFFERs\n");
    exit(EXIT_FAILURE);
  }
 
  /* call MPI send and store handle. */

  MPI_Arg Count = size;
  MPI_Arg Dest = Target;
  MPI_Arg Mtag = Tag;
 
  stat = MPI_Isend(buffer_send, Count, Type, Dest, Mtag, CommWorld, RequestHandle+index);
  if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
  // Uncommenting the next line can improve performance in some cases.
  // MPI_Wait(RequestHandle+index, &Status);
 
  /* Store buffer info. */
 
  RequestBuffer[index] = (char *) buffer_send;
  LastActiveIndex = max(LastActiveIndex, index);
 
  return SUCCESS;
}
 
#endif /* USE_MPI */
