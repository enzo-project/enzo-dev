/***********************************************************************
/
/  COMMUNICATION ROUTINE: BROADCAST A VALUE FROM ROOT TO OTHERS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
void my_exit(int status);
 
 
int CommunicationBroadcastValue(Eint32 *Value, int BroadcastProcessor)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Datatype DataTypeInt = MPI_INT;
  MPI_Arg Count = 1;
  MPI_Arg Root = BroadcastProcessor;
  MPI_Arg stat;

  stat = MPI_Bcast((void*) Value, Count, DataTypeInt, Root, MPI_COMM_WORLD);
    if( stat != MPI_SUCCESS ){ENZO_FAIL("");}
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[15]+= endtime-starttime;
  counter[15] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationBroadcastValue(Eint64 *Value, int BroadcastProcessor)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Datatype DataTypeInt = MPI_LONG_LONG_INT;
  MPI_Arg Count = 1;
  MPI_Arg Root = BroadcastProcessor;
  MPI_Arg stat;

  stat = MPI_Bcast((void*) Value, Count, DataTypeInt, Root, MPI_COMM_WORLD);
    if( stat != MPI_SUCCESS ){ENZO_FAIL("");}
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[15]+= endtime-starttime;
  counter[15] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationBroadcastValues(Eint32 *Values, int Number, int BroadcastProcessor)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Datatype DataTypeInt = MPI_INT;
  MPI_Arg Count = Number;
  MPI_Arg Root = BroadcastProcessor;
  MPI_Arg stat;

  stat = MPI_Bcast((void*) Values, Count, DataTypeInt, Root, MPI_COMM_WORLD);
    if( stat != MPI_SUCCESS ){ENZO_FAIL("");}
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[15]+= endtime-starttime;
  counter[15] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationBroadcastValues(Eint64 *Values, int Number, int BroadcastProcessor)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Datatype DataTypeInt = MPI_LONG_LONG_INT;
  MPI_Arg Count = Number;
  MPI_Arg Root = BroadcastProcessor;
  MPI_Arg stat;

  stat = MPI_Bcast((void*) Values, Count, DataTypeInt, Root, MPI_COMM_WORLD);
    if( stat != MPI_SUCCESS ){ENZO_FAIL("");}
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[15]+= endtime-starttime;
  counter[15] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
#endif /* USE_MPI */
 
  return SUCCESS;
}
