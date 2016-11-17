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

void Star::ActivateNewStar(FLOAT Time, float Timestep)
{
  int StarType;
  FILE *fptr;
  if (this->IsUnborn()) {  // unborn
    StarType = ABS(type);
    switch (StarType) {
    case SimpleSource:
      type = StarType;
      break;
    case PopII:
      if (Mass >= StarClusterMinimumMass) {
	type = StarType;
	if (StarClusterUnresolvedModel)
	  BirthTime = Time-Timestep;
	else
	  // slightly before to avoid round-off errors in comparisons
	  BirthTime = (1-1e-6)*Time;
      }
      break;
    case PopIII:
      if (Mass >= this->FinalMass) {
	type = StarType;
	BirthTime = (1-1e-6)*Time;
      }
      break;
    case PopIII_CF:
      type = StarType;
      break;
    case BlackHole:
      // nothing to do
      break;
    case MBH:  // probably you wouldn't need this activation routine anyway (see mbh_maker)
      type = StarType;  
      break;
    } // ENDSWITCH type
  } // ENDIF FORMATION

  return;
}
