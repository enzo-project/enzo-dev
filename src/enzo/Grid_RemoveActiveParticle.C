/***********************************************************************
/
/  REMOVE AN ACTIVE PARTICLE BY ITS ID
/
/  written by: Nathan Goldbaum
/  date:       December 2011
/  modified1:
/
/  RETURNS: 0 = not found; 1 = removed
/
************************************************************************/

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "ActiveParticle.h"

int grid::RemoveActiveParticle(PINT ID, int NewProcessorNumber)
{

  int i,j,found = FALSE;

  if (MyProcessorNumber != ProcessorNumber)
    return found;

  if (NumberOfActiveParticles == 0)
    return found;

  for (i=0; i < NumberOfActiveParticles; i++)
    if (this->ActiveParticles[i]->ReturnID() == ID) {
      found = TRUE;
      break;
    }
  
  if (found == FALSE)
    return found;

  this->NumberOfActiveParticles--;
  this->ActiveParticles.erase(i);

  return found;

}
