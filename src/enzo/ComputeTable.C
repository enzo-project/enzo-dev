/***********************************************************************
/
/  EVALUATE A TABLE OF VALUES OF THE GIVEN FUNCTION
/
/  written by: Greg Bryan
/  date:       June, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
 
#define DONT_USE_LOCK
#define UNLOCKED 0
#define LOCKED 1
 
static int Lock = UNLOCKED;
 
int ComputeTable(float Min, float Max, float Step, float (*Function)(float),
		 float **Table, float *TableStep, float *TableMin,
		 float *TableMax)
{
 
#ifdef USE_LOCK
  float a;
  while (Lock == LOCKED) {
    for (int l = 0; l < 10000; l++)
      a += POW(-1, l)*sin(l);
  }
 
  Lock = LOCKED;
#endif /* USE_LOCK */
 
  /* Check if the table that has previously been constructed (if there is one)
     is applicable for the requested table.  If it is, we're done. */
 
  if (*Table != NULL)
    if (*TableMin <= Min && *TableMax >= Max && Step <= *TableStep) {
      Lock = UNLOCKED;
      return SUCCESS;
    }
 
  /* If there was a previous table, delete it. */
 
  if (*Table != NULL)
    delete *Table;
  else {
    *TableMin  = Min;
    *TableMax  = Max;
    *TableStep = Step;
  }
 
  /* Set Table values to use the most demanding of the parameters. */
 
  *TableMin  = min(*TableMin, Min);
  *TableMax  = max(*TableMax, Max);
  *TableStep = min(*TableStep, Step);
 
  /* Allocate room for table. */
 
  int Space = int((*TableMax - *TableMin)/(*TableStep)) + 1;
  if (debug)
    printf("ComputeTable: building table with max = %"FSYM", n = %"ISYM".\n",
	   *TableMax, Space);
  *Table = new float[Space];
 
  /* Fill out table. */
 
  for (int i = 0; i < Space; i++)
    *(*Table + i) = (*Function)(Min + float(i)*Step);
 
  /* Done. */
 
  if (debug)
    printf("ComputeTable: done building table.\n");
 
#ifdef USE_LOCK
  Lock = UNLOCKED;
#endif /* USE_LOCK */
  return SUCCESS;
}
