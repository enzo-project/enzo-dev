#ifdef DEFINE_STORAGE
#define EXTERN
#else
#define EXTERN extern
#endif /* DEFINE_STORAGE */

#define MAX_PH_RECEIVE_BUFFERS 50000
#define MAX_PH_REQUESTS 50000

// Maximum number of photons to send in one MPI message.
#define PHOTON_BUFFER_SIZE 100000

#define NO_TRANSPORT 0
#define TRANSPORT 1
#define HALT_TRANSPORT 2
#define SENT_DATA 3
#define RECV_DATA 4
#define EMPTY_BUFFER -1
#define BUFFER_TRASH -1
#define BUFFER_END -99999
#define NO_HINT -1

#ifdef USE_MPI

EXTERN char *PH_CommunicationReceiveBuffer[MAX_PH_RECEIVE_BUFFERS];
EXTERN MPI_Request PH_CommunicationReceiveMPI_Request[MAX_PH_RECEIVE_BUFFERS];
EXTERN Eint32 PH_CommunicationReceiveIndex;
EXTERN Eint32 PH_CommunicationReceiveMaxIndex;

EXTERN Eint32 PhotonMessageBuffer[MAX_PH_RECEIVE_BUFFERS];
EXTERN MPI_Request PhotonMessageRequest[MAX_PH_RECEIVE_BUFFERS];
EXTERN Eint32 PhotonMessageIndex;
EXTERN Eint32 PhotonMessageMaxIndex;

EXTERN char KeepTransMessageBuffer[10*MAX_PH_RECEIVE_BUFFERS];
EXTERN MPI_Request KeepTransMessageRequest[10*MAX_PH_RECEIVE_BUFFERS];
EXTERN Eint32 KeepTransMessageIndex;
EXTERN Eint32 KeepTransMessageMaxIndex;

#ifdef UNUSED
struct MyRequest {
  int nbuffer;
  Eint32 buffer[MAX_PH_REQUESTS];
  MPI_Request requests[MAX_PH_REQUESTS];
};
#endif /* UNUSED */

#endif /* USE_MPI */
