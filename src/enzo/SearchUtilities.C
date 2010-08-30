/***********************************************************************
/
/  BASIC SEARCHING UTILITIES
/
/  written by: John Wise
/  date:       August, 2010
/  modified:   
/
/  PURPOSE:
/
************************************************************************/
// std::lower_bound has too much overhead

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"

int search_lower_bound(int *arr, int value, int low, int high, 
		       int total)
{
  int mid, width;
  if (high < low)
    ENZO_VFAIL("high (%d) < low (%d) when searching for lower bound.",
	       high, low);
  width = high-low;
  mid = low + width/2;
  // First catch here if it's the last recursive call
  //printf("low high width = %d %d %d\n", low, high, width);
  if (width <= 1) {
    if (mid == 0 || mid == total-1)
      return mid;
    else if (value >= arr[mid] && value < arr[mid+1])
      return mid;
    else // implied if (value >= arr[mid-1] && value < arr[mid])
      return mid-1;
  } // ENDIF width <= 1
  if (arr[mid] > value)
    return search_lower_bound(arr, value, low, mid-1, total);
  else if (arr[mid] < value)
    return search_lower_bound(arr, value, mid+1, high, total);
  else
    return mid;  // found it exactly.
}
