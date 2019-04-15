/***********************************************************************
/
/  GRID CLASS (APPEND ACTIVE PARTICLE DATA TO AN ACTIVE PARTICLE ARRAY)
/
/  written by: Nathan Goldbaum
/  date:       March, 2012
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

int grid::AppendActiveParticlesToList(
    ActiveParticleList<ActiveParticleType>& APArray, int search_id) {
  
  // Return if this does not concern us
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int PNum;

  for (PNum = 0; PNum < NumberOfActiveParticles; PNum++) 
    if (search_id == ActiveParticles[PNum]->ReturnType()) 
      APArray.copy_and_insert(*ActiveParticles[PNum]);
      
  return SUCCESS;
} 
