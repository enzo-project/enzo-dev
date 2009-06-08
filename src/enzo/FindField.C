/***********************************************************************
/
/  FIND FIELD FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: field index or -1 on failure
/
************************************************************************/
 
// Find field type field in array field_type, returning the index into the
//   field array or -1 if it is not there.
 
#include "macros_and_parameters.h"
#include "typedefs.h"
 
 
int FindField(int field, int farray[], int numfields)
{
  for (int i = 0; i < numfields; i++)
    if (field == farray[i])
      return i;
 
  /* not found */
 
  return -1;
}
