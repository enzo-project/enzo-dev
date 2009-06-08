#if defined(SP2)

// Memory usage from getrusage

#include<stdio.h>
#include<sys/time.h>
#include<sys/resource.h>
#include<unistd.h>


long long int mused(void)
{

#ifdef MEM_TRACE
  struct rusage temp;
  long long int bytes;
  int result;

  result = getrusage(RUSAGE_SELF, &temp);
  if( result == 0 ) {
    bytes = ((long long int) (1024)) * ((long long int) temp.ru_maxrss);
  } else {
    bytes = ((long long int) (0));
  }
  return(bytes);
#else
  return((long long int) 0);
#endif

}

#elif defined(LINUX) || defined(IA64)

// Memory usage from proc tables

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
void my_exit(int status);

long long int mused(void)
{

#ifdef MEM_TRACE
  pid_t this_pid = getpid();
  int ps = getpagesize();

  char procname[120];
  long long int kb;
  FILE* ptr;

  sprintf(procname, "/proc/%lld/statm", ((long long int) this_pid));
  ptr = fopen(procname, "r");
  fscanf(ptr, "%lld", &kb);
  fclose(ptr);
  fprintf(stderr, "Proc statm Bytes: %lld\n", ((long long int) kb*ps));
  return ((long long int) kb*ps);
#else
  return ((long long int) 0);
#endif

}

#else

// Default case return zero bytes

long long int mused(void)
{
  return((long long int) 0);
}

#endif
