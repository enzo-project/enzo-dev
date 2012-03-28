/***********************************************************************
/
/  COMMUNICATION ROUTINE: DETERMINE NUMBER OF NODES
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

int DetermineNumberOfNodes(void)
{

  if (NumberOfProcessors == 1) {
    CoresPerNode = 1;
    return NumberOfProcessors;
  }

  int i, j, err, NumberOfNodes = 0;
  char *AllHosts, *UniqueHosts, hostname[MAX_LINE_LENGTH];
  char host1[MAX_LINE_LENGTH], host2[MAX_LINE_LENGTH];
  bool unique;

  hostname[MAX_LINE_LENGTH-1] = '\0';
  err = gethostname(hostname, MAX_LINE_LENGTH-1);

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    AllHosts = new char[NumberOfProcessors * MAX_LINE_LENGTH];
    UniqueHosts = new char[NumberOfProcessors * MAX_LINE_LENGTH];
  } // ENDIF ROOT PROCESSOR

#ifdef USE_MPI
  MPI_Gather(hostname, MAX_LINE_LENGTH, MPI_BYTE, 
	     AllHosts, MAX_LINE_LENGTH, MPI_BYTE,
	     ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    // Count unique hostnames
    for (i = 0; i < NumberOfProcessors; i++) {

      // Extract hostname from array
      memcpy(host1, AllHosts+i*MAX_LINE_LENGTH, MAX_LINE_LENGTH);

      //printf("%"ISYM" :: %s\n", i, host1);
      unique = true;
      for (j = 0; j < NumberOfNodes; j++) {
	memcpy(host2, UniqueHosts+j*MAX_LINE_LENGTH, MAX_LINE_LENGTH);
	if (strcmp(host1, host2) == 0) {
	  unique = false;
	  break;
	}
      } // ENDFOR j

      if (unique) {
	memcpy(UniqueHosts+NumberOfNodes*MAX_LINE_LENGTH, 
	       host1, MAX_LINE_LENGTH);
	NumberOfNodes++;
      }
    } // ENDFOR processors

//    for (i = 0; i < NumberOfNodes; i++) {
//      memcpy(host1, UniqueHosts+i*MAX_LINE_LENGTH, MAX_LINE_LENGTH);
//      printf("Node %"ISYM" :: %s\n", i, host1);
//    }
    
    delete [] UniqueHosts;
    delete [] AllHosts;

  } // ENDIF ROOT_PROCESSOR

  CommunicationBroadcastValue(&NumberOfNodes, ROOT_PROCESSOR);
  CoresPerNode = NumberOfProcessors / NumberOfNodes;

  if (debug)
    printf("DetermineNumberOfNodes: %"ISYM" nodes, %"ISYM" cores per node\n", 
	   NumberOfNodes, CoresPerNode);

  return NumberOfNodes;

}
