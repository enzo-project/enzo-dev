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

#define NO_OVERLOAD_NEW
#define NO_MALLOC_REPORT
#define MALLOC_REPORT_FREQUENCY 100

/************************************************************************
 *  NEW/DELETE OVERLOAD, ADDING MEMORY LOGGING
 ************************************************************************/
 
#ifdef OVERLOAD_NEW
 
int CurrentMemoryUsage = 0;   // in words
int MaximumMemoryUsage = 0;
int NumberOfCalls      = 0;
void *FirstAddress     = NULL;
void *LargestAddress   = NULL;
 
void* operator new(size_t NumberOfBytes)
{
 
  if (NumberOfBytes == 0)
    return NULL;
 
  CurrentMemoryUsage += NumberOfBytes;
  MaximumMemoryUsage = max(MaximumMemoryUsage, CurrentMemoryUsage);
  NumberOfCalls++;
 
#ifdef MALLOC_REPORT
  if (NumberOfCalls % MALLOC_REPORT_FREQUENCY == 0)
    printf("new_malloc: Current = %"GSYM"   Max = %"GSYM"\n",
	   float(CurrentMemoryUsage),
	   float(MaximumMemoryUsage));
#endif /* MALLOC_REPORT */
 
  void *pointer = malloc(NumberOfBytes+sizeof(float));
 
  if (pointer == NULL) {
    fprintf(stderr, "Error allocating %"ISYM" bytes.\n", NumberOfBytes);
    exit(EXIT_FAILURE);
  }
 
  if (FirstAddress == NULL)
    FirstAddress = pointer;
  LargestAddress = max(LargestAddress, pointer);
 
  *((float *) pointer) = float(NumberOfBytes);
 
  return (void *) (((float *) pointer) + 1);
 
}
 
 
void operator delete(void *pointer)
{
  if (pointer == NULL){
    return;}
 
  CurrentMemoryUsage -= *(((float *) pointer) - 1);
 
  free(((float *) pointer) - 1);
 
  return;
}
 
#endif /* OVERLOAD_NEW */

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
