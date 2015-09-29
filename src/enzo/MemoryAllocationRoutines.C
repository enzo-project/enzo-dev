/***********************************************************************
/
/  OVERLOAD THE NEW/DELETE OPERATORS FOR FLOAT AND INT
/
/  written by: Greg Bryan
/  date:       March, 1996
/  modified1:  Robert Harkness
/  date:       March, 2004
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <new>
#ifdef USE_JEMALLOC
#define JEMALLOC_MANGLE
#include <jemalloc/jemalloc.h>
#undef JEMALLOC_MANGLE
#endif /* USE_JEMALLOC */
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#define NO_MEMORY_TRACE
#define NO_MALLOC_REPORT
#define MALLOC_REPORT_FREQUENCY 100

/************************************************************************
 *  NEW/DELETE OVERLOAD, ADDING MEMORY LOGGING AND/OR LOG2ALLOC
 ************************************************************************/

#if defined(MEMORY_TRACE) || defined(USE_LOG2ALLOC)

#ifdef USE_LOG2ALLOC
#define LOG2ALLOC(size) ((size_t)pow(2, (int)log2(size-1)+1))
#else
#define LOG2ALLOC(size) (size)
#endif

#ifdef MEMORY_TRACE
int CurrentMemoryUsage = 0;   // in words
int MaximumMemoryUsage = 0;
int NumberOfCalls      = 0;
void *FirstAddress     = NULL;
void *LargestAddress   = NULL;
#endif /* MEMORY_TRACE */
 
void* operator new(size_t NumberOfBytes)
{
 
  NumberOfBytes = LOG2ALLOC(NumberOfBytes);

  if (NumberOfBytes == 0)
    return NULL;

#ifdef MEMORY_TRACE
  CurrentMemoryUsage += NumberOfBytes;
  MaximumMemoryUsage = max(MaximumMemoryUsage, CurrentMemoryUsage);
  NumberOfCalls++;
 
#ifdef MALLOC_REPORT
  if (NumberOfCalls % MALLOC_REPORT_FREQUENCY == 0)
    printf("new_malloc: Current = %"GSYM"   Max = %"GSYM"\n",
	   float(CurrentMemoryUsage),
	   float(MaximumMemoryUsage));
#endif /* MALLOC_REPORT */
#endif /* MEMORY_TRACE */
 
  void *pointer = malloc(NumberOfBytes+sizeof(float));
 
  if (pointer == NULL) {
    fprintf(stderr, "Error allocating %"ISYM" bytes.\n", NumberOfBytes);
    exit(EXIT_FAILURE);
  }
 
#ifdef MEMORY_TRACE
  if (FirstAddress == NULL)
    FirstAddress = pointer;
  LargestAddress = max(LargestAddress, pointer);
#endif /* MEMORY_TRACE */

  *((float *) pointer) = float(NumberOfBytes);
 
  return (void *) (((float *) pointer) + 1);
 
}
 

#ifdef MEMORY_TRACE
void operator delete(void *pointer)
{
  if (pointer == NULL){
    return;}
 
  CurrentMemoryUsage -= *(((float *) pointer) - 1);
 
  free(((float *) pointer) - 1);
 
  return;
}
#endif /* MEMORY_TRACE */

#endif /* defined(MEMORY_TRACE) || defined(USE_LOG2ALLOC) */

/************************************************************************
 *  NEW/DELETE OVERLOAD, USING AN EXTERNAL LIBRARY, JEMALLOC
 ************************************************************************/

#ifdef USE_JEMALLOC

void* operator new(size_t NumberOfBytes) throw (std::bad_alloc) {
  if (NumberOfBytes == 0) return NULL;
  //void *pointer = jemalloc(NumberOfBytes + sizeof(float));
  void *pointer = jemalloc(NumberOfBytes);
  if (pointer == NULL)
    ENZO_VFAIL("Error allocation %d bytes", NumberOfBytes);
  return pointer;
  //*((float*) pointer) = float(NumberOfBytes);
  //return (void*) (((float*) pointer) + 1);
}

void* operator new[](size_t NumberOfBytes) throw (std::bad_alloc) {
  if (NumberOfBytes == 0) return NULL;
  //void *pointer = jemalloc(NumberOfBytes + sizeof(float));
  void *pointer = jemalloc(NumberOfBytes);
  if (pointer == NULL)
    ENZO_VFAIL("Error allocation %d bytes", NumberOfBytes);
  return pointer;
  //*((float*) pointer) = float(NumberOfBytes);
  //return (void*) (((float*) pointer) + 1);
}

void operator delete(void *pointer) throw() {
  if (pointer == NULL) return;
  //jefree( ((float*) pointer) - 1);
  jefree(pointer);
  return;
}

void operator delete[](void *pointer) throw() {
  if (pointer == NULL) return;
  //jefree( ((float*) pointer) - 1);
  jefree(pointer);
  return;
}

#endif
