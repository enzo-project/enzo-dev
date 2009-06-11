/***********************************************************************
/
/  REPORT MEMORY USAGE (MALLOC, ETC)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness
/  date:       March, 2004
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
 
#define NO_OVERLOAD_NEW
 
#ifdef OVERLOAD_NEW
extern int CurrentMemoryUsage;   // in words
extern int MaximumMemoryUsage;
extern int NumberOfCalls     ;
extern void *FirstAddress;
extern void *LargestAddress;
#endif /* OVERLOAD_NEW */
 
 
int ReportMemoryUsage(char *header = NULL)
{
 
#ifdef OVERLOAD_NEW
 
  float Current = float(CurrentMemoryUsage),
        Maximum = float(MaximumMemoryUsage),
        Arena   = float(((int *) LargestAddress -
			 (int *) FirstAddress))*sizeof(int);
  printf("%s: P(%"ISYM"):CurrentMemoryUsage = %"GSYM"  Max = %"GSYM"  Arena = %"GSYM" (%"GSYM"%%)\n",
	 header, MyProcessorNumber, Current, Maximum, Arena,
	 Maximum/Arena*100.0);
 
#endif /* OVERLOAD_NEW */
 
  return SUCCESS;
}
