/***********************************************************************
/
/  POP III IMF LOOKUP TABLE INITIALIZATION
/
/  written by: John Wise
/  date:       April, 2010
/  modified1:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

void mt_init(unsigned_int seed);
unsigned_long_int mt_random(void);

int StarParticlePopIII_IMFInitialize(void)
{

  const float CutoffExponent = 1.6;

  if (IMFData != NULL)
    return SUCCESS;

  IMFData = new float[IMF_TABLE_ENTRIES];

  int i;
  float m, m0, dm, total_fn;

  dm = log10(PopIIIUpperMassCutoff / PopIIILowerMassCutoff) / 
    (float) (IMF_TABLE_ENTRIES-1);
  m0 = log10(PopIIILowerMassCutoff);
  total_fn = 0;
  for (i = 0; i < IMF_TABLE_ENTRIES; i++) {
    m = POW(10.0, m0 + i*dm);
    total_fn += POW(m, PopIIIInitialMassFunctionSlope) * 
      exp(-POW((PopIIIStarMass / m), CutoffExponent));
    IMFData[i] = total_fn;
  } // ENDFOR i

  // Normalize
  for (i = 0; i < IMF_TABLE_ENTRIES; i++)
    IMFData[i] /= IMFData[IMF_TABLE_ENTRIES-1];

  /* Initialize the random number generator.  Call it 1+NumberOfCalls
     to get the generator back to the same state as it was before the
     restart (if any). */

  if (PopIIIInitialMassFunctionSeed == INT_UNDEFINED)
    mt_init(time(NULL)); //+100*MyProcessorNumber);
  else
    mt_init(PopIIIInitialMassFunctionSeed); //+100*MyProcessorNumber);
  unsigned_long_int trash;
  for (i = 0; i < 1+PopIIIInitialMassFunctionCalls; i++)
    trash = mt_random();

  return SUCCESS;

}
