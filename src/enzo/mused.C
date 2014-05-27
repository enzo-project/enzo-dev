#ifdef USE_MPI
#include "mpi.h"
#endif

#if defined(SP2) || defined(__APPLE__)

// Memory usage from getrusage
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CommunicationUtilities.h"


Eint64 mused(void)
{

#define MEM_TRACE

#ifdef MEM_TRACE
  struct rusage temp;
#if defined (__APPLE__)
  Eint64 bytes = 1;
#else
  Eint64 bytes = 1024;
#endif
  int result;

  result = getrusage(RUSAGE_SELF, &temp);
  if( result == 0 ) {
    bytes *= (Eint64) temp.ru_maxrss;
  } else {
    bytes = ((Eint64) (0));
  }
  return(bytes);
#else
  return((Eint64) 0);
#endif

}

#elif defined(LINUX) || defined(IA64)

// Memory usage from proc tables

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CommunicationUtilities.h"

void my_exit(int status);

Eint64 mused(void)
{

#ifdef MEM_TRACE
  pid_t this_pid = getpid();
  int ps = getpagesize();

  char procname[120];
  Eint64 kb, res;
  FILE* ptr;

  sprintf(procname, "/proc/%lld/statm", ((Eint64) this_pid));
  ptr = fopen(procname, "r");
  fscanf(ptr, "%lld %lld", &kb, &res);
  fclose(ptr);
  //fprintf(stderr, "Proc statm Bytes: %lld\n", ((Eint64) kb*ps));
  return ((Eint64) ps*res);
#else
  return ((Eint64) 0);
#endif

}

#else

#include <stdio.h>
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CommunicationUtilities.h"

// Default case return zero bytes

Eint64 mused(void)
{
  return((Eint64) 0);
}

#endif

void PrintMemoryUsage(char *str)
{
#ifdef MEM_TRACE
  Eint64 MemInUse, MinMem, MaxMem, MeanMem;
  MemInUse = mused();
  MinMem = CommunicationMinValue(MemInUse);
  MaxMem = CommunicationMaxValue(MemInUse);
  MeanMem = MemInUse;
  CommunicationSumValues(&MeanMem, 1);
  MeanMem /= NumberOfProcessors;
  fprintf(memtracePtr, "%s  %16lld\n", str, MemInUse);
  if (debug) 
    printf("MEM %s :: mean/min/max = %lld/%lld/%lld\n", 
	   str, MeanMem, MinMem, MaxMem);
#endif
  return;
}
