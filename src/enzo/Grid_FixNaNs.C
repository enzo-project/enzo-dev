/***********************************************************************
/
/  GRID CLASS (Fix NaNs with crude interpolation)
/
/  written by:
/  date:
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
 
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
 
#define FIX_NANS
 
int grid::FixNaNs(const char *message)
{
  
#ifdef FIX_NANS
 
  // Return if this doesn't concern us
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  // Set this to zero so the compiler doesn't optimize everything away
 
  int ThisIsZero = GridStartIndex[0] - NumberOfGhostZones;
  int num, k, j, i, index;
  float x1, x2, y1, y2, z1, z2, newval;

  for (num = 0; num < NumberOfBaryonFields; num++) {
    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
        for (i = 0; i < GridDimension[0]; i++) {

          index = GRIDINDEX(i, j, k);
          if (BaryonField[num][index+ThisIsZero] != BaryonField[num][index]) {

            x1 = BaryonField[num][GRIDINDEX(i-1, j, k)];
            x2 = BaryonField[num][GRIDINDEX(i+1, j, k)];
            y1 = BaryonField[num][GRIDINDEX(i, j-1, k)];
            y2 = BaryonField[num][GRIDINDEX(i, j+1, k)];
            z1 = BaryonField[num][GRIDINDEX(i, j, k-1)];
            z2 = BaryonField[num][GRIDINDEX(i, j, k+1)];
            newval = ( x1 + x2 + y1 + y2 + z1 + z2 ) / 6.0;
            BaryonField[num][index] = newval;

            fprintf(stderr, "FixNaNs[%s](Proc %"ISYM"): gas %"ISYM" %"ISYM" %"GSYM"\n", message,
                    MyProcessorNumber, num, index, BaryonField[num][index]);
            fprintf(stderr, "    FixNaNs: x1, x2 - %"GSYM", %"GSYM"\n", x1, x2);
            fprintf(stderr, "    FixNaNs: y1, y2 - %"GSYM", %"GSYM"\n", y1, y2);
            fprintf(stderr, "    FixNaNs: z1, z2 - %"GSYM", %"GSYM"\n", z1, z2);

          }

        }
      }
    }
  }
  
#endif /* FIX_NANS */
		
  return SUCCESS;
}
