/***********************************************************************
/
/  STOCHASTIC FORCING CLASS METHOD: CommunicationBroadcastFlags
/
/  written by: Wolfram Schmidt
/  date:       October, 2005
/  modified1:  January, 2009
/  modified2: Sep, 2014: modified to support Enzo 2.4 // P. Grete
/
/  PURPOSE: broadcasts flags identifying non-zero modes 
/           from root to others
/
************************************************************************/


#include "preincludes.h"
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "StochasticForcing.h"
#include "global_data.h"
#include "error.h"

/* function prototypes */

void my_exit(int status);

void StochasticForcing::CommunicationBroadcastFlags(void)
{
  if (NumberOfProcessors == 1) return;
    
#ifdef USE_MPI
 
  if (debug) printf("Broadcasting flags, proc #%"ISYM"\n",MyProcessorNumber);
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  MPI_Arg stat;

  stat = MPI_Bcast(mask, NumModes, DataTypeInt, ROOT_PROCESSOR, MPI_COMM_WORLD);
    if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[15]+= endtime-starttime;
  counter[15] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
#endif /* USE_MPI */
}
