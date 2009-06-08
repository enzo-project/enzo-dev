#ifdef DEFINE_STORAGE
#define EXTERN
#else
#define EXTERN extern
#endif /* DEFINE_STORAGE */

#ifdef USE_MPI

#define MAX_PH_RECEIVE_BUFFERS 50000
#define MAX_PH_REQUESTS 50000

#define HALT_TRANSPORT 2
#define EMPTY_BUFFER -1
#define BUFFER_TRASH -1
#define BUFFER_END -99999

EXTERN char *PH_CommunicationReceiveBuffer[MAX_PH_RECEIVE_BUFFERS];
EXTERN MPI_Request PH_CommunicationReceiveMPI_Request[MAX_PH_RECEIVE_BUFFERS];
EXTERN Eint32 PH_CommunicationReceiveIndex;

struct MyRequest {
  int nbuffer;
  Eint32 buffer[MAX_PH_REQUESTS];
  MPI_Request requests[MAX_PH_REQUESTS];
};

#endif /* USE_MPI */
