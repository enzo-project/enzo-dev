/***********************************************************************
/
/  SEARCH FOR NEW PARTICLES THAT EXCEED THE MINIMUM STAR MASS AND 
/  ACTIVATE THEM
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "StarParticleData.h"

void Star::ActivateNewStar(FLOAT Time)
{
  int StarType;
  FILE *fptr;
  if (this->IsUnborn()) {  // unborn
    StarType = ABS(type);
    switch (StarType) {
    case PopII:
      if (Mass > StarClusterMinimumMass) {
	type = StarType;
	BirthTime = Time;
      }
      break;
    case PopIII:
      if (Mass >= PopIIIStarMass) {
	type = StarType;  // No minimum mass now.  User-specified mass.
	BirthTime = Time;
      }
      break;
    case PopIII_CF:
      type = StarType;
      BirthTime = Time;
      break;
    case BlackHole:
      // nothing to do
      break;
    } // ENDSWITCH type
  } // ENDIF FORMATION

  return;
}
