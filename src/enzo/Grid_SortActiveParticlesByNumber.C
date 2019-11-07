/***********************************************************************
/
/  GRID CLASS (SORT ACTIVE PARTICLES BY UNIQUE ID)
/
/  written by: Nathan Goldbaum
/  date:       Jul, 2012
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ActiveParticle.h"


void grid::SortActiveParticlesByNumber()
{
  ActiveParticles.sort_number(0, NumberOfActiveParticles);
}
