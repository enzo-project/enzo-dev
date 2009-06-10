#define DEBUG 0
/***********************************************************************
/
/  NON-BLOCKING COMMUNICATION ROUTINES FOR PHOTONS
/
/  written by: John H. Wise
/  date:       December, 2007
/  modified1:  
/
/  PURPOSE:
/
/  TODO: Stop using MPI_BYTE for communication, and create an MPI
/        structure
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "GroupPhotonList.h"
#include "PhotonCommunication.h"

#define TOOLARGE 10000000

#ifdef USE_MPI
static MyRequest *PhotonNumberReceiveHandles = NULL;
static MyRequest *KeepTransportingHandles = NULL;
#endif /* USE_MPI */

int InitiateKeepTransportingCheck(int keep_transporting)
{

#ifdef USE_MPI
  int proc, index, i;
  MPI_Request dummy_req;

  /* If this is the first call, initialize the keep_transporting
     request array and buffer. */

  if (KeepTransportingHandles == NULL) {
    KeepTransportingHandles = new MyRequest[NumberOfProcessors];
    for (proc = 0; proc < NumberOfProcessors; index++) {
      KeepTransportingHandles[proc].nbuffer = 0;
      for (i = 0; i < MAX_PH_REQUESTS; i++) {
	KeepTransportingHandles[proc].requests[i] = NULL;
	KeepTransportingHandles[proc].buffer[i] = EMPTY_BUFFER;
      } // ENDFOR buffers
    } // ENDFOR processors
  } // ENDIF initialize KeepTransporting

  /* Now initiate the communication */

  index = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++) {

    if (proc != MyProcessorNumber) {

      // Find open spot in the receive queue
      for (i = 0; i < MAX_PH_REQUESTS; i++)
	if (KeepTransportingHandles[proc].buffer[i] == EMPTY_BUFFER) {
	  index = i;
	  break;
      }

      if (DEBUG)
	printf("EP[%"ISYM"]: Sending keep_transporting = %"ISYM" to P%"ISYM" (index %"ISYM")\n",
	       MyProcessorNumber, keep_transporting, proc, index);
      MPI_Irecv(&KeepTransportingHandles[proc].buffer[index], 1, MPI_INT, proc, 
		MPI_KEEPTRANSPORTING_TAG, MPI_COMM_WORLD,
		&KeepTransportingHandles[proc].requests[index]);
      KeepTransportingHandles[proc].nbuffer++;

      MPI_Isend(&keep_transporting, 1, MPI_INT, proc, MPI_KEEPTRANSPORTING_TAG,
		MPI_COMM_WORLD, &dummy_req);
      if (dummy_req != MPI_REQUEST_NULL)
	MPI_Request_free(&dummy_req);
    } // ENDIF other processor
    
  } // ENDFOR processors

#endif /* USE_MPI */
  
  return SUCCESS;

}

/**********************************************************************/
int StopKeepTransportingCheck()
{

#ifdef USE_MPI
  int proc, index, i;
  MPI_Request dummy_req;
  int HaltMessage = HALT_TRANSPORT;

  index = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++) {

    if (proc != MyProcessorNumber) {

      // Find open spot in the receive queue
      for (i = 0; i < MAX_PH_REQUESTS; i++)
	if (KeepTransportingHandles[proc].buffer[i] == EMPTY_BUFFER) {
	  index = i;
	  break;
      }
      
      if (DEBUG)
	printf("EP[%"ISYM"]: Sending keep_transporting = %"ISYM" to P%"ISYM" (index %"ISYM")\n",
	       MyProcessorNumber, HaltMessage, proc, index);

      MPI_Isend(&HaltMessage, 1, MPI_INT, proc, MPI_KEEPTRANSPORTING_TAG,
		MPI_COMM_WORLD, &dummy_req);
      if (dummy_req != MPI_REQUEST_NULL)
	MPI_Request_free(&dummy_req);
    } // ENDIF other processor
    
  } // ENDFOR processors

#endif /* USE_MPI */
  
  return SUCCESS;

}

/**********************************************************************/

int KeepTransportingCheck(int &keep_transporting)
{

#ifdef USE_MPI
  int i, index, proc;
  int *ReceivedMessage = new int[NumberOfProcessors];
  Eint32 ListOfIndices[MAX_PH_REQUESTS];
  Eint32 RequestDone;
  int ReceivedFromNProcs = 1;  // start with one ... the local processor
  int NumberOfReceives, NumberOfRequests, NumberOfBuffers;

  for (index = 0; index < NumberOfProcessors; index++)
    ReceivedMessage[index] = FALSE;
  ReceivedMessage[MyProcessorNumber] = TRUE;
  
  while (keep_transporting == 0 && ReceivedFromNProcs < NumberOfProcessors) {

    for (proc = 0; proc < NumberOfProcessors; proc++) {
    
      if (proc != MyProcessorNumber) {

	/* Check how many requests have been made in this processor
	   and adjust nbuffer accordingly */

	NumberOfBuffers = KeepTransportingHandles[proc].nbuffer;
	NumberOfRequests = 0;
	NumberOfReceives = 0;

	for (i = 0; i < NumberOfBuffers; i++)
	  if (KeepTransportingHandles[proc].requests[i] != NULL)
	    NumberOfRequests++;

	for (i = 0; i < NumberOfBuffers; i++) {
	  if (KeepTransportingHandles[proc].requests[i] != NULL) {

	    MPI_Test(&KeepTransportingHandles[proc].requests[i], 
		     &RequestDone, MPI_STATUSES_IGNORE);

	    if (RequestDone) {
	      keep_transporting = max(keep_transporting, 
				      KeepTransportingHandles[proc].buffer[i]);
	      if (DEBUG)
		printf("EP[%"ISYM"]: Received keep_transporting = %"ISYM" from P%"ISYM"-i%"ISYM"; "
		       "k_t = %"ISYM"\n", 
		       MyProcessorNumber, KeepTransportingHandles[proc].buffer[i], 
		       proc, index, keep_transporting);
	      KeepTransportingHandles[proc].nbuffer--;
	      KeepTransportingHandles[proc].buffer[i] = EMPTY_BUFFER;
	      if (KeepTransportingHandles[proc].requests[i] != MPI_REQUEST_NULL)
		MPI_Request_free(&KeepTransportingHandles[proc].requests[i]);
	      KeepTransportingHandles[proc].requests[i] = NULL;
	      NumberOfReceives++;

	    }  // ENDIF RequestDone
	    
	  } // ENDIF request != NULL

	} // ENDFOR buffers

	//	if (NumberOfReceives > 0 && ReceivedMessage[proc] == FALSE &&
	if ((NumberOfReceives == NumberOfRequests && NumberOfReceives > 0 &&
	    ReceivedMessage[proc] == FALSE) ||
	    keep_transporting == HALT_TRANSPORT) {
	  ReceivedFromNProcs++;
	  ReceivedMessage[proc] = TRUE;
	}

      } // ENDIF other processors

    } // ENDFOR processors

  } // ENDWHILE not(keep_transporting) and haven't received all of the
    // messages

  /* Cancel and free any unreceived messages */

//  for (index = 0; index < NumberOfProcessors-1; index++)
//    if (ReceivedMessage[index] == FALSE) {
//      MPI_Cancel(&KeepTransportingRequest[index]);
//      if (KeepTransportingRequest[index] != MPI_REQUEST_NULL)
//        MPI_Request_free(&KeepTransportingRequest[index]);
//    } // ENDIF no received message
  
#endif /* USE_MPI */
  return SUCCESS;

}

/**********************************************************************/

int InitializePhotonCommunication()
{

#ifdef USE_MPI
  int i, index, proc;

  /* Initialize the keep_transporting communication buffer */

  if (KeepTransportingHandles != NULL)
    delete [] KeepTransportingHandles;

  KeepTransportingHandles = new MyRequest[NumberOfProcessors];
  for (proc = 0; proc < NumberOfProcessors; proc++) {
    KeepTransportingHandles[proc].nbuffer = 0;
    for (i = 0; i < MAX_PH_REQUESTS; i++) {
      KeepTransportingHandles[proc].requests[i] = NULL;
      KeepTransportingHandles[proc].buffer[i] = EMPTY_BUFFER;
    } // ENDFOR buffers
  } // ENDFOR processors

  /* If this is the first call, initialize the request array */

  if (PhotonNumberReceiveHandles != NULL)
    delete [] PhotonNumberReceiveHandles;

  PhotonNumberReceiveHandles = new MyRequest[NumberOfProcessors];
  for (proc = 0; proc < NumberOfProcessors; proc++) {
    PhotonNumberReceiveHandles[proc].nbuffer = 0;
    for (i = 0; i < MAX_PH_REQUESTS; i++) {
      PhotonNumberReceiveHandles[proc].requests[i] = NULL;
      PhotonNumberReceiveHandles[proc].buffer[i] = EMPTY_BUFFER;
    } // ENDFOR requests
  } // ENDFOR processors
#endif /* USE_MPI */

  return SUCCESS;

}

/**********************************************************************/

int InitiatePhotonNumberSend(int *nPhoton)
{

#ifdef USE_MPI
  int i, index, proc;
  MPI_Request dummy_req;

  /* If this is the first call, initialize the request array */

  if (PhotonNumberReceiveHandles == NULL) {
    PhotonNumberReceiveHandles = new MyRequest[NumberOfProcessors];
    for (proc = 0; proc < NumberOfProcessors; proc++) {
      PhotonNumberReceiveHandles[proc].nbuffer = 0;
      for (i = 0; i < MAX_PH_REQUESTS; i++) {
	PhotonNumberReceiveHandles[proc].requests[i] = NULL;
	PhotonNumberReceiveHandles[proc].buffer[i] = EMPTY_BUFFER;
      } // ENDFOR requests
    } // ENDFOR processors
  } // ENDIF ReceiveHandles == NULL

  /* Post nonblocking receive and send calls for all processors */

  for (proc = 0; proc < NumberOfProcessors; proc++)
    if (proc != MyProcessorNumber) {

      // Find open spot in the receive queue
      for (i = 0; i < MAX_PH_REQUESTS; i++)
	if (PhotonNumberReceiveHandles[proc].buffer[i] == EMPTY_BUFFER) {
	  index = i;
	  break;
	}

      // Post receive for nPhoton_RECV
      MPI_Irecv(&PhotonNumberReceiveHandles[proc].buffer[index], 1, 
		MPI_INT, proc, MPI_NPHOTON_TAG, MPI_COMM_WORLD,
		&PhotonNumberReceiveHandles[proc].requests[index]);
      PhotonNumberReceiveHandles[proc].nbuffer++;

      // Send nPhoton
      if (nPhoton[proc] > TOOLARGE || nPhoton[proc] < 0) {
	printf("CTPh[P%"ISYM"]: WARNING -- Number of sent photons = %"ISYM" bad?\n",
	       MyProcessorNumber, nPhoton[proc]);
	ENZO_FAIL("Error in: "__FILE__);
      }
      MPI_Isend(nPhoton+proc, 1, MPI_INT, proc, MPI_NPHOTON_TAG, MPI_COMM_WORLD,
		&dummy_req);
      if (dummy_req != MPI_REQUEST_NULL)
	MPI_Request_free(&dummy_req);

    } // ENDIF other processor
#endif /* USE_MPI */

  return SUCCESS;

}

/**********************************************************************/

int InitializePhotonReceive(int group_size)
{

#ifdef USE_MPI
  int i, index, nPhoton_RECV, proc;
  Eint32 NumberOfReceives, NumberOfRequests;
  int tag;
  int FinishedAny = FALSE;
  Eint32 ListOfIndices[MAX_PH_REQUESTS];
  GroupPhotonList *RecvList = NULL;

  //  while (FinishedAny == FALSE) {
    for (proc = 0; proc < NumberOfProcessors; proc++) {
      if (proc != MyProcessorNumber) {

	/* Check if any this processor has received any nPhoton
	   messages */

//	printf("PHR[%"ISYM"]: Expected %"ISYM" receives from P%"ISYM"\n", MyProcessorNumber,
//	       PhotonNumberReceiveHandles[proc], proc);

	NumberOfRequests = 0;
	for (i = 0; i < PhotonNumberReceiveHandles[proc].nbuffer; i++) {
	  if (PhotonNumberReceiveHandles[proc].requests[i] != NULL)
	    NumberOfRequests++;
	}

	PhotonNumberReceiveHandles[proc].nbuffer = NumberOfRequests;
	MPI_Waitsome(NumberOfRequests,
		     PhotonNumberReceiveHandles[proc].requests,
		     &NumberOfReceives, ListOfIndices, MPI_STATUSES_IGNORE);

	for (i = 0; i < NumberOfReceives; i++) {

	  /* Get number of photons that we're receiving and free the
	     MPI request */

	  index = ListOfIndices[i];

	  if (PhotonNumberReceiveHandles[proc].requests[index] == MPI_REQUEST_NULL) {

	    PhotonNumberReceiveHandles[proc].nbuffer--;
	    nPhoton_RECV = PhotonNumberReceiveHandles[proc].buffer[index];
	    PhotonNumberReceiveHandles[proc].buffer[index] = EMPTY_BUFFER;
	    if (PhotonNumberReceiveHandles[proc].requests[index] != MPI_REQUEST_NULL)
	      MPI_Request_free(&PhotonNumberReceiveHandles[proc].requests[index]);

	    /* Now we know how many photons to receive, post the receive
	       call for them */
	
	    if (nPhoton_RECV > 0) {
	      if (nPhoton_RECV > TOOLARGE) {
		printf("CTPh[P%"ISYM"]: WARNING -- Number of received photons"
		       " = %"ISYM" too large?\n", MyProcessorNumber, nPhoton_RECV);
		ENZO_FAIL("Error in: "__FILE__);
	      }
	      RecvList = new GroupPhotonList[nPhoton_RECV];
	      tag = MPI_PHOTONGROUP_TAG*10+nPhoton_RECV;
	      if (DEBUG)
		printf("CTPh[P%"ISYM"]: Receiving %"ISYM" photons from P%"ISYM" "
		       "(TAG=%"ISYM", index = %"ISYM", %"ISYM"/%"ISYM")\n", 
		       MyProcessorNumber, nPhoton_RECV, proc, tag, index, i+1, 
		       NumberOfReceives);
	      MPI_Irecv(RecvList, group_size*nPhoton_RECV, MPI_BYTE, proc,
			tag, MPI_COMM_WORLD,
			PH_CommunicationReceiveMPI_Request + 
			PH_CommunicationReceiveIndex);
	      PH_CommunicationReceiveBuffer[PH_CommunicationReceiveIndex] =
		(char *) RecvList;
	      PH_CommunicationReceiveIndex++;
	    } else
	      RecvList = NULL;

	    FinishedAny = TRUE;

	  } // ENDIF request is NULL

	} // ENDFOR receives

      } // ENDIF other processor
    } // ENDFOR processors
    //  } // ENDWHILE received no nPhoton_RECV messages

#endif /* USE_MPI */

  return SUCCESS;

}
