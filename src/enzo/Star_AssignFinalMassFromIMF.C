/***********************************************************************
/
/  ASSIGN FINAL MASS TO STAR OBJECT
/
/  written by: John Wise
/  date:       April, 2010
/  modified1: 
/
/  PURPOSE: When we're using an IMF, assign the this->FinalMass when we
/           first create the star object.
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

unsigned_long_int mt_random(void);

int Star::AssignFinalMassFromIMF(float TimeUnits)
{

  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  float x = (float) (random_int%max_random) / (float) (max_random);
  float dm = log10(PopIIIUpperMassCutoff / PopIIILowerMassCutoff) / 
    (float) (IMF_TABLE_ENTRIES-1);

  /* (binary) search for the mass bin corresponding to the random
     number */

  int width = IMF_TABLE_ENTRIES/2;
  int bin_number = IMF_TABLE_ENTRIES/2;
  
  while (width > 1) {
    width /= 2;
    if (x > IMFData[bin_number])
      bin_number += width;
    else if (x < IMFData[bin_number])
      bin_number -= width;
    else
      break;
  }
  
  this->FinalMass = PopIIILowerMassCutoff * POW(10.0, bin_number * dm);

  /* Adjust the lifetime (taken from the fit in Schaerer 2002) now we
     know the stellar mass.  It was set to the lifetime of a star with
     M=PopIIILowerMassCutoff in pop3_maker.src as a placeholder. */

  float logm = log10((float)this->FinalMass);

  // First in years, then convert to code units
  this->LifeTime = POW(10.0, (9.785 - 3.759*logm + 1.413*logm*logm - 
			      0.186*logm*logm*logm)) / (TimeUnits/3.1557e7);

//  printf("random_num = %f, mass = %f Msun, lifetime = %f Myr\n",
//	 x, this->FinalMass, this->LifeTime * TimeUnits / 3.1557e13);

  PopIIIInitialMassFunctionCalls++;

  return SUCCESS;

}
