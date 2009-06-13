/***********************************************************************
/
/  CREATES STAR PARTICLES FROM EXISTING PARTICLES
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  NOTES:  negative types mark particles that have just been before 
/          and not been converted into a star particle.
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

int grid::UpdateStarParticles(int level)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0 || Stars == NULL)
    return SUCCESS;

  int i;
  Star *cstar;

  for (cstar = Stars; cstar; cstar = cstar->NextStar)
    if (cstar->type > 0)  // living stars only (<0 == waiting to be created)
      for (i = 0; i < NumberOfParticles; i++)
	if (cstar->Identifier == ParticleNumber[i]) {
	  cstar->CopyFromParticle(this, i, level);
	  break;
	} // ENDIF matched ID

  return SUCCESS;

}
