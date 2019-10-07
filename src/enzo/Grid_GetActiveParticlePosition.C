/***********************************************************************
/
/  GRID CLASS (RETURN ACTIVE PARTICLE POSITIONS AS A 2D POINTER ARRAY)
/
/  written by: Nathan Goldbaum
/  date:       November, 2012
/  modified1:  
/
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

void grid::GetActiveParticlePosition(FLOAT *ActiveParticlePosition[]) 
{
  int i, dim;

  for (i = 0; i < NumberOfActiveParticles; i++) {
    FLOAT* pos = ActiveParticles[i]->ReturnPosition();
    for (dim = 0; dim < GridRank; dim++)
      ActiveParticlePosition[dim][i] = pos[dim];
  }

}
