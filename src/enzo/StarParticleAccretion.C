/***********************************************************************
/
/  ADDS ACCRETED MASS TO PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  NOTES:
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
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

int StarParticleAccretion(Star *&AllStars)
{

#define NOT_SEDOV_TEST
#define NOT_HII_REGION_TEST

#if defined(SEDOV_TEST) || defined(HII_REGION_TEST)
  return SUCCESS;
#endif

  Star *ThisStar;

  /* Add accreted mass to star particles */

  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar) {

    if (ThisStar->CalculateMassAccretion() == FAIL) {
      fprintf(stderr, "Error in star::CalculateMassAccretion.\n");
      ENZO_FAIL("");
    }

    if (ThisStar->Accrete() == FAIL) {
      fprintf(stderr, "Error in star::Accrete.\n");
      ENZO_FAIL("");
    }
    
  }

  return SUCCESS;

}
