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
#include "CommunicationUtilities.h"
#include "PhotonCommunication.h"

#ifdef USE_MPI
int CommunicationBufferedSendCancel(int Tag);
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
int CommunicationFindOpenRequest(MPI_Request *requests, Eint32 last_free,
				 Eint32 nrequests, Eint32 index, Eint32 &max_index);
void CommunicationCheckForErrors(int NumberOfStatuses, MPI_Status *statuses,
				 char *msg=NULL);
int KeepTransportingSend(int keep_transporting);
int KeepTransportingCheck(char* &kt_global, int &keep_transporting);
int InitializePhotonMessages(void);
static Eint32 PH_ListOfIndices[MAX_PH_RECEIVE_BUFFERS];
static MPI_Status PH_ListOfStatuses[MAX_PH_RECEIVE_BUFFERS];
#endif /* USE_MPI */

/******************************************************************/

int InitializePhotonCommunication(void)
{
#ifdef USE_MPI
  int i, proc;

  /* Receive any orphaned messages from the last RT call, so they
     don't interfere with this cycle. */

#ifdef NONBLOCKING_RT
  char *dummy = new char[NumberOfProcessors];
  KeepTransportingCheck(dummy, i);
  delete [] dummy;
  CommunicationBarrier();
#endif /* NONBLOCKING_RT */

  /* Initialize */

  PhotonMessageIndex = 0;
  PhotonMessageMaxIndex = 0;
  KeepTransMessageIndex = 0;
  KeepTransMessageMaxIndex = 0;
  PH_CommunicationReceiveIndex = 0;  
  PH_CommunicationReceiveMaxIndex = 0;

  for (i = 0; i < MAX_PH_RECEIVE_BUFFERS; i++)
    PH_CommunicationReceiveMPI_Request[i] = NULL;
  
  /* Anticipate the first call for the number of messages to receive */

#ifdef NONBLOCKING_RT
  InitializePhotonMessages();
#endif

#endif /* USE_MPI */
    return SUCCESS;
}

/**********************************************************************/
int InitializePhotonMessages(void)
{
#ifdef USE_MPI
  int proc;
  PhotonMessageIndex = 0;
  PhotonMessageMaxIndex = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++)
    if (proc != MyProcessorNumber) {
      MPI_Irecv(PhotonMessageBuffer+PhotonMessageIndex,
		1, MPI_INT, proc, MPI_NPHOTON_TAG, MPI_COMM_WORLD,
		PhotonMessageRequest+PhotonMessageIndex);
      PhotonMessageIndex++;
      PhotonMessageMaxIndex++;
    } // ENDIF other processor
#endif /* USE_MPI */
  return SUCCESS;
}
/**********************************************************************/

int KeepTransportingInitialize(char* &kt_global, bool initial_call)
{
#ifdef USE_MPI
  int proc;
  kt_global = new char[NumberOfProcessors];
  for (proc = 0; proc < NumberOfProcessors; proc++)
    kt_global[proc] = 1;
  for (proc = 0; proc < NumberOfProcessors; proc++)
    if (proc != MyProcessorNumber) {
      MPI_Irecv(KeepTransMessageBuffer+KeepTransMessageIndex, 
		1, MPI_CHAR, proc, 
		MPI_KEEPTRANSPORTING_TAG, MPI_COMM_WORLD,
		KeepTransMessageRequest+KeepTransMessageIndex);
      if (initial_call) {
	KeepTransMessageIndex++;
	KeepTransMessageMaxIndex++;
      } else {
	KeepTransMessageIndex =
	  CommunicationFindOpenRequest(KeepTransMessageRequest, NO_HINT,
				       10*MAX_PH_RECEIVE_BUFFERS,
				       KeepTransMessageIndex,
				       KeepTransMessageMaxIndex);
      }
      if (DEBUG)
	printf("P%"ISYM": Sending KT=%"ISYM" to P%"ISYM"\n", MyProcessorNumber, 
	       kt_global[MyProcessorNumber], proc);
      CommunicationBufferedSend(kt_global+MyProcessorNumber, 1, MPI_CHAR, proc,
				MPI_KEEPTRANSPORTING_TAG, MPI_COMM_WORLD, 1);

      } // ENDIF other processor
#endif /* USE_MPI */
  return SUCCESS;
}

int KeepTransportingFinalize(char* &kt_global, int keep_transporting)
{
#ifdef USE_MPI
  /* Send a halt message to all other processes to ensure that they
     exit the keep_transporting loop, only if this process exited the
     loop not from a halt signal. */

  if (keep_transporting != 2)
    KeepTransportingSend(HALT_TRANSPORT);

  CommunicationBarrier();
  KeepTransportingCheck(kt_global, keep_transporting);

  delete [] kt_global;

#endif /* USE_MPI */
  return SUCCESS;
}

int FinalizePhotonCommunication(void)
{
#ifdef USE_MPI
  int i;

  /* If there are any leftover MPI_Irecv calls from the photon
     "nPhoton" and keep_transporting calls, cancel them. */

#ifdef NONBLOCKING_RT
  for (i = 0; i < PhotonMessageMaxIndex; i++)
    if (PhotonMessageRequest[i] != MPI_REQUEST_NULL)
      MPI_Cancel(PhotonMessageRequest+i);

  for (i = 0; i < KeepTransMessageMaxIndex; i++)
    if (KeepTransMessageRequest[i] != MPI_REQUEST_NULL)
      MPI_Cancel(KeepTransMessageRequest+i);
#endif /* NONBLOCKING_RT */

  for (i = 0; i < PH_CommunicationReceiveMaxIndex; i++)
    if (PH_CommunicationReceiveMPI_Request[i] != MPI_REQUEST_NULL)
      MPI_Cancel(PH_CommunicationReceiveMPI_Request+i);

  /* Wait until all of the requests are cancelled */

  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

#ifdef NONBLOCKING_RT
  MPI_Waitall(PhotonMessageMaxIndex, PhotonMessageRequest, 
	      PH_ListOfStatuses);
  CommunicationCheckForErrors(PhotonMessageMaxIndex, PH_ListOfStatuses,
			      "Waitall photon message cancels");
  MPI_Waitall(KeepTransMessageMaxIndex, KeepTransMessageRequest, 
	      PH_ListOfStatuses);
  CommunicationCheckForErrors(KeepTransMessageMaxIndex, PH_ListOfStatuses,
			      "Waitall KT message cancels");
#endif /* NONBLOCKING_RT */

  MPI_Waitall(PH_CommunicationReceiveMaxIndex, 
	      PH_CommunicationReceiveMPI_Request,
	      PH_ListOfStatuses);
  CommunicationCheckForErrors(PH_CommunicationReceiveMaxIndex, PH_ListOfStatuses,
			      "Waitall KT message cancels");

  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

  /* Cancel all buffered header sends */

//  CommunicationBufferedSendCancel(MPI_KEEPTRANSPORTING_TAG);
//  CommunicationBufferedSendCancel(MPI_NPHOTON_TAG);

#endif /* USE_MPI */
  return SUCCESS;
}

/**********************************************************************/

int CommunicationNumberOfPhotonSends(int *nPhoton, int size)
{
#ifdef USE_MPI

  /* Send the number of messages that we're going to send, so the
     client can post the receives. */

  Eint32 NumberOfMessages, proc;

  for (proc = 0; proc < NumberOfProcessors; proc++)
#ifdef NONBLOCKING_RT
    if (proc != MyProcessorNumber && nPhoton[proc] > 0) {
#else
    if (proc != MyProcessorNumber) {
#endif
      NumberOfMessages = nPhoton[proc] / size;
      if (nPhoton[proc] % size > 0) NumberOfMessages++;
      if (DEBUG)
	printf("CRP[%"ISYM"]: Sending %"ISYM" messages, %"ISYM" photons to P%"ISYM"\n", MyProcessorNumber,
	       NumberOfMessages, nPhoton[proc], proc);
      CommunicationBufferedSend(&NumberOfMessages, 1, MPI_INT, proc,
				MPI_NPHOTON_TAG, MPI_COMM_WORLD, sizeof(Eint32));
    }
#endif
  return SUCCESS;
}


/**********************************************************************/

#ifdef USE_MPI
int PostPhotonReceives(Eint32 index, Eint32 proc, int size, MPI_Datatype type)
{
  int i;
  Eint32 NumberOfMessages;
  GroupPhotonList *ReceiveBuffer;
  
  NumberOfMessages = PhotonMessageBuffer[index];
  for (i = 0; i < NumberOfMessages; i++) {
    if (DEBUG)
      printf("CTPh[P%"ISYM"]: Receiving photons from P%"ISYM" "
	     "(message %"ISYM", %"ISYM"/%"ISYM", index %"ISYM")\n", 
	     MyProcessorNumber, proc, i, i+1, NumberOfMessages,
	     PH_CommunicationReceiveIndex);
    ReceiveBuffer = new GroupPhotonList[size];

    MPI_Irecv(ReceiveBuffer, size, type, proc,
	      MPI_PHOTONGROUP_TAG, MPI_COMM_WORLD,
	      PH_CommunicationReceiveMPI_Request +
	      PH_CommunicationReceiveIndex);

    PH_CommunicationReceiveBuffer[PH_CommunicationReceiveIndex] =
      (char *) ReceiveBuffer;

#ifdef NONBLOCKING_RT
    PH_CommunicationReceiveIndex = 
      CommunicationFindOpenRequest(PH_CommunicationReceiveMPI_Request,
				   NO_HINT, MAX_PH_RECEIVE_BUFFERS,
				   PH_CommunicationReceiveIndex,
				   PH_CommunicationReceiveMaxIndex);
#else
    PH_CommunicationReceiveIndex++;
    PH_CommunicationReceiveMaxIndex++;
#endif

    if (DEBUG)
      printf("P%d: PHCRIndex/Max = %d/%d\n", MyProcessorNumber, 
	     PH_CommunicationReceiveIndex, PH_CommunicationReceiveMaxIndex);

  } // ENDFOR i (messages)

  return SUCCESS;
}
#endif /* USE_MPI */

/**********************************************************************/
#ifdef USE_MPI
int CommunicationFindOpenRequest(MPI_Request *requests, Eint32 last_free,
				 Eint32 nrequests, Eint32 index,
				 Eint32 &max_index)
{
  int i, count;
  bool FoundOpenRequest = false;
  MPI_Status status;

  max_index = max(max_index, index);
  i = (last_free != NO_HINT) ? last_free : 0;
  i = i % nrequests;
  count = 0;

  while (!FoundOpenRequest) {
    if (requests[i] == NULL || requests[i] == MPI_REQUEST_NULL)
      FoundOpenRequest = true;
    else {
      i++;
      i = i % nrequests;
      if (count++ > MAX_PH_RECEIVE_BUFFERS)
	ENZO_VFAIL("Exceeded number (%"ISYM") of comm. buffers MAX_PH_RECEIVE_BUFFERS :: %"ISYM" %"ISYM" %"ISYM" %x", count, last_free, i, nrequests, requests[i]);
    }
  } // ENDWHILE

  max_index = max(max_index, i);
  return i;
}
#endif /* USE_MPI */
/**********************************************************************/

#ifdef USE_MPI
int InitializePhotonReceive(int max_size, bool local_transport,
			    MPI_Datatype MPI_PhotonType)
{

  MPI_Status status;
  MPI_Arg proc, MessageReceived;
  Eint32 NumberOfMessages, NumberOfReceives;
  GroupPhotonList *ReceiveBuffer = NULL;
  int i, j, index, RecvProc;

  /* Receive MPI messages that contain how many messages with the
     actual photon data that we'll be receiving from each process. */

  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  if (local_transport)
    MPI_Testsome(PhotonMessageMaxIndex, PhotonMessageRequest, &NumberOfReceives,
		 PH_ListOfIndices, PH_ListOfStatuses);
  else {
    MPI_Waitall(PhotonMessageMaxIndex, PhotonMessageRequest,
		PH_ListOfStatuses);
    NumberOfReceives = PhotonMessageMaxIndex;
  }
  if (NumberOfReceives > 0)
    CommunicationCheckForErrors(PhotonMessageMaxIndex, PH_ListOfStatuses,
				"Testsome InitializePhotonReceive");
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

  if (DEBUG && NumberOfReceives > 0)
    printf("P%"ISYM": Received %"ISYM" header messages, Index/MaxIndex = %"ISYM"/%"ISYM".\n", 
	   MyProcessorNumber, NumberOfReceives, PhotonMessageIndex, 
	   PhotonMessageMaxIndex);

  for (i = 0; i < NumberOfReceives; i++) {

#ifdef NONBLOCKING_RT    
    index = PH_ListOfIndices[i];
#else
    index = i;
#endif
    RecvProc = PH_ListOfStatuses[i].MPI_SOURCE;

    if (DEBUG)
      printf("P%"ISYM": Processing header message %"ISYM" (null=%"ISYM") from P%"ISYM".\n", 
	     MyProcessorNumber, index, 
	     (PhotonMessageRequest[index] == MPI_REQUEST_NULL),
	     RecvProc);

    PostPhotonReceives(index, RecvProc, max_size, MPI_PhotonType);

    /* Post another receive for the next loop.  Immediately test if
       there's already a message waiting.  If so, post another
       receive.  Repeat until no more messages. */

#ifdef NONBLOCKING_RT
    do {
      PhotonMessageIndex =
	CommunicationFindOpenRequest(PhotonMessageRequest, index,
				     MAX_PH_RECEIVE_BUFFERS, PhotonMessageIndex,
				     PhotonMessageMaxIndex);
      MPI_Irecv(PhotonMessageBuffer + PhotonMessageIndex,
		1, MPI_INT, RecvProc, MPI_NPHOTON_TAG, MPI_COMM_WORLD,
		PhotonMessageRequest + PhotonMessageIndex);
      MPI_Test(PhotonMessageRequest + PhotonMessageIndex, &MessageReceived,
	       MPI_STATUS_IGNORE);

      if (MessageReceived) {
//	CommunicationCheckForErrors(1, &status, 
//				    "InitializePhotonReceive immediate");
	PostPhotonReceives(PhotonMessageIndex, RecvProc, max_size, MPI_PhotonType);
      }

    } while (MessageReceived);
#endif /* NONBLOCKING_RT */

  } // ENDFOR i (receives)

  return SUCCESS;
}
#endif

/************************************************************************/
int KeepTransportingSend(int keep_transporting)
{
#ifdef USE_MPI
  int proc;
  char value = keep_transporting;
  for (proc = 0; proc < NumberOfProcessors; proc++)
    if (proc != MyProcessorNumber) {
      if (DEBUG)
	printf("P%"ISYM": Sending KT=%"ISYM" to P%"ISYM"\n", MyProcessorNumber,
	       value, proc);
      CommunicationBufferedSend(&value, 1, MPI_CHAR, proc,
				MPI_KEEPTRANSPORTING_TAG,
				MPI_COMM_WORLD, 1);
    } // ENDIF other processor
#endif /* USE_MPI */
  return SUCCESS;
}

/************************************************************************/
int KeepTransportingCheck(char* &kt_global, int &keep_transporting)
{
#ifdef USE_MPI
  int i, index, RecvProc;
  char value = keep_transporting;
  char received = RECV_DATA;
  bool PingRequired, PingReceived, AcceptMessage;
  MPI_Arg proc, NumberOfReceives, MessageReceived;
  MPI_Status status;

//  if (DEBUG)
//    printf("P%"ISYM": keep_transporting(before) = %"ISYM", KTMaxIndex = %"ISYM"\n", 
//	   MyProcessorNumber, keep_transporting, KeepTransMessageMaxIndex);

  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
  MPI_Testsome(KeepTransMessageMaxIndex, KeepTransMessageRequest,
	       &NumberOfReceives, PH_ListOfIndices, PH_ListOfStatuses);
  if (NumberOfReceives > 0)
    CommunicationCheckForErrors(KeepTransMessageMaxIndex, PH_ListOfStatuses,
				"KTCheck Testsome");
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

  if (DEBUG && NumberOfReceives > 0)
    printf("P%"ISYM": Received %"ISYM" KT messages, Index/MaxIndex = %"ISYM"/%"ISYM".\n", 
	   MyProcessorNumber, NumberOfReceives, KeepTransMessageIndex, 
	   KeepTransMessageMaxIndex);

  int second_recv, next_kt;
  for (i = 0; i < NumberOfReceives; i++) {
    AcceptMessage = true;
    index = PH_ListOfIndices[i];
    RecvProc = PH_ListOfStatuses[i].MPI_SOURCE;
    //if (RecvProc < 0) continue;   // Undefined rank
    PingRequired = (kt_global[RecvProc] == SENT_DATA);
    if (PingRequired) {
      AcceptMessage = (KeepTransMessageBuffer[index] == RECV_DATA);
      if (AcceptMessage) PingRequired = false;
      next_kt = -1;
    } else {
      if (KeepTransMessageBuffer[index] != RECV_DATA)
	next_kt = KeepTransMessageBuffer[index];
      else
	next_kt = -1;
    } // ENDELSE

    // Halt message overrides everything
    if (KeepTransMessageBuffer[KeepTransMessageIndex] == HALT_TRANSPORT)
      next_kt = 2;

    if (DEBUG)
      printf("P%"ISYM": Primary KT receive, P%"ISYM", = %"ISYM"\n",
	     MyProcessorNumber, RecvProc,
	     KeepTransMessageBuffer[KeepTransMessageIndex]);

    /* Post another receive for the next loop.  Immediately test if
       there's already a message waiting.  If so, post another
       receive.  Repeat until no more messages. */

    second_recv = 0;
    do {
      KeepTransMessageIndex =
	CommunicationFindOpenRequest(KeepTransMessageRequest, NO_HINT,
				     10*MAX_PH_RECEIVE_BUFFERS,
				     KeepTransMessageIndex,
				     KeepTransMessageMaxIndex);
      MPI_Irecv(KeepTransMessageBuffer + KeepTransMessageIndex, 
		1, MPI_CHAR, RecvProc,
		MPI_KEEPTRANSPORTING_TAG, MPI_COMM_WORLD,
		KeepTransMessageRequest + KeepTransMessageIndex);
      MPI_Test(KeepTransMessageRequest + KeepTransMessageIndex, 
	       &MessageReceived, MPI_STATUS_IGNORE);
      if (MessageReceived) {
//	CommunicationCheckForErrors(1, &status, 
//				    "KT immediate");
	if (PingRequired) {
	  AcceptMessage = 
	    (KeepTransMessageBuffer[KeepTransMessageIndex] == RECV_DATA);
	  if (AcceptMessage) PingRequired = false;
	}
	else if (AcceptMessage) {
	  if (KeepTransMessageBuffer[KeepTransMessageIndex] != RECV_DATA)
	    next_kt = KeepTransMessageBuffer[KeepTransMessageIndex];
	    //next_kt = max(next_kt, KeepTransMessageBuffer[KeepTransMessageIndex]);
	  else
	    next_kt = -1;   // Forget any older messages
	}

	// Halt message overrides everything
	if (KeepTransMessageBuffer[KeepTransMessageIndex] == HALT_TRANSPORT)
	  next_kt = 2;

	if (DEBUG)
	  printf("P%"ISYM": Secondary KT receive, P%"ISYM", Recv %"ISYM", = %"ISYM"\n",
		 MyProcessorNumber, RecvProc, second_recv,
		 KeepTransMessageBuffer[KeepTransMessageIndex]);
	second_recv++;
      }
    } while (MessageReceived);

    // Ping back the processor, saying that we've received this flag
    if (kt_global[RecvProc] == SENT_DATA) {
      CommunicationBufferedSend(&received, 1, MPI_CHAR, RecvProc,
				MPI_KEEPTRANSPORTING_TAG, MPI_COMM_WORLD, 1);
      CommunicationBufferedSend(&value, 1, MPI_CHAR, RecvProc,
				MPI_KEEPTRANSPORTING_TAG, MPI_COMM_WORLD, 1);
    }

    if (next_kt >= 0)
      kt_global[RecvProc] = next_kt;
    else
      kt_global[RecvProc] = TRANSPORT;

    if (DEBUG)
      printf("P%"ISYM": Setting kt_global[%"ISYM"] = %"ISYM".  %"ISYM" secondary receives\n", 
	     MyProcessorNumber, RecvProc, kt_global[RecvProc], second_recv);

  } // ENDFOR i (receives)

  // Find max in the global keep_transporting array
  kt_global[MyProcessorNumber] = keep_transporting;
  keep_transporting = 0;

  for (proc = 0; proc < NumberOfProcessors; proc++) {
    if (kt_global[proc] == SENT_DATA) 
      keep_transporting = max(keep_transporting, 1);
    else
      keep_transporting = max(keep_transporting, kt_global[proc]);
  }

  // Only exit if we're finished with all of our local work.
  if (keep_transporting == HALT_TRANSPORT && kt_global[MyProcessorNumber] == 1)
    keep_transporting = 1;

  if (DEBUG && NumberOfReceives > 0)
    printf("P%"ISYM": keep_transporting = %"ISYM"/%"ISYM", kt_global = %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",
	   MyProcessorNumber, value, keep_transporting, kt_global[0],
	   kt_global[1], kt_global[2], kt_global[3], kt_global[4], 
	   kt_global[5], kt_global[6], kt_global[7]);
#endif /* USE_MPI */
  return SUCCESS;
}

/**********************************************************************/

#ifdef USE_MPI
void CommunicationCheckForErrors(int NumberOfStatuses, MPI_Status *statuses,
				 char *msg)
{
  int i;
  char error_string[1024];
  MPI_Arg length, error_class, datasize;
  for (i = 0; i < NumberOfStatuses; i++)
    if (statuses[i].MPI_ERROR != 0) {
      MPI_Get_count(statuses+i, MPI_BYTE, &datasize);
      fprintf(stderr, "MPI Error %"ISYM" processor %"ISYM"\n",
	      statuses[i].MPI_ERROR, MyProcessorNumber);
      fprintf(stderr, "\t MPI_TAG = %"ISYM", MPI_SOURCE = %"ISYM", datasize = %"ISYM" bytes\n",
	      statuses[i].MPI_TAG, statuses[i].MPI_SOURCE, datasize);
      if (msg != NULL)
	fprintf(stderr, "P%"ISYM": error occurred at %s\n", MyProcessorNumber, msg);
      MPI_Error_class(statuses[i].MPI_ERROR, &error_class);
      MPI_Error_string(error_class, error_string, &length);
      fprintf(stderr, "P%"ISYM": %s\n", MyProcessorNumber, error_string);
      MPI_Error_string(statuses[i].MPI_ERROR, error_string, &length);
      fprintf(stderr, "P%"ISYM": %s\n", MyProcessorNumber, error_string);
      ENZO_FAIL("MPI communication error!");

    } // ENDIF MPI_ERROR
  return;
}
#endif /* USE_MPI */
