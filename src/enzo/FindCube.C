#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
 
int FindCube(char *cubename)
{
  int i;
 
//  if ( Unigrid ) {
    for( i = 0; i < MAX_CUBE_DUMPS; i++ )
    {
       if( CubeDumps[i] != NULL )
       {
         if( strcmp(cubename, CubeDumps[i]) == 0 )
         {
            //  fprintf(stderr, "Found %s, element %"ISYM"\n", CubeDumps[i], i);
            return i;
         }
       }
    }
//  }
 
//  fprintf(stderr, "Didn't find :%s:\n", cubename);
  return -1;
 
}
