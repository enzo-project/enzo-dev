/***********************************************************************
/
/  (WRAPPER) ADD X-RAY EMISSION FROM RADIATIVE PARTICLES
/
/  written by: John Wise
/  date:       March, 2012
/  modified1:
/
/  PURPOSE:
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
#include "CosmologyParameters.h"
#include "Star.h"

int grid::AddOpticallyThinXrays(Star *AllStars, int NumberOfSources)
{

  /* If we're merging rays, we already have a binary tree of the
     sources.  We can use that to speed up the calculation when we
     have more than 10 sources.  With smaller numbers, the overhead
     makes the direct calculation faster. */

#ifdef UNUSED // Difficult with multiple X-ray energy bins
  if (RadiativeTransferSourceClustering == TRUE && NumberOfSources >= 10)
    this->AddXraysFromTree();
  else
#endif
  this->AddXraysFromSources(AllStars);

  HasRadiation = TRUE;

  return SUCCESS;

}
