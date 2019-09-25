/***********************************************************************
/
/  GRID CLASS (DEBUG ACTIVE PARTICLE DATA)
/
/  written by: John Wise
/  date:       Februrary, 2012
/  modified1:  
/
/  PURPOSE:
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "fortran.def"
#include "CosmologyParameters.h"

#include "ActiveParticle.h"

int grid::DebugActiveParticles(int level)
{

  if (ProcessorNumber != MyProcessorNumber || NumberOfActiveParticles == 0)
    return SUCCESS;

  int i, inside;
  FLOAT *pos;

  /* Check if the active particles are within the grid */

  if (ActiveParticles.size() != NumberOfActiveParticles) {
    printf("Active particle count mismatch!");
    printf("NumberOfActiveParticles = %"GOUTSYM"\n", NumberOfActiveParticles);
    printf("ActiveParticles.size() = %"GOUTSYM"\n", ActiveParticles.size());
    ENZO_FAIL("")
  }

  for (i = 0; i < ActiveParticles.size(); i++) {
    pos = ActiveParticles[i]->ReturnPosition();
    inside = this->PointInGrid(pos);
    if (inside == FALSE) {
      printf("Active particle outside grid!  level %d, grid %d\n", level, this->ID);
      printf("\t pos       = %"GOUTSYM " %"GOUTSYM " %"GOUTSYM "\n", pos[0], pos[1], pos[2]);
      printf("\t left edge = %"GOUTSYM " %"GOUTSYM " %"GOUTSYM "\n",
	     GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
      printf("\t right edge = %"GOUTSYM " %"GOUTSYM " %"GOUTSYM "\n",
	     GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
      ENZO_FAIL("");
    }

  }

  return SUCCESS;

}
