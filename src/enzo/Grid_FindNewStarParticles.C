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

void InsertStarAfter(Star * &Node, Star * &NewNode);

int grid::FindNewStarParticles(int level)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0)
    return SUCCESS;

  int i;
  Star *NewStar, *cstar;
  bool exists;

  for (i = 0; i < NumberOfParticles; i++)
    if (ParticleType[i] == -PARTICLE_TYPE_SINGLE_STAR ||
       (ParticleType[i] == PARTICLE_TYPE_STAR && UseSupernovaSeedFieldSourceTerms &&
	  (this->Time >= ParticleAttribute[0][i] &&
	   this->Time <= ParticleAttribute[0][i]+ParticleAttribute[1][i])) ||
	ParticleType[i] == -PARTICLE_TYPE_BLACK_HOLE ||
	ParticleType[i] == -PARTICLE_TYPE_CLUSTER ||
	ParticleType[i] == -PARTICLE_TYPE_COLOR_STAR ||
	ParticleType[i] == -PARTICLE_TYPE_SIMPLE_SOURCE ||
	ABS(ParticleType[i]) == PARTICLE_TYPE_MBH ||
	ParticleType[i] == -PARTICLE_TYPE_SUPERNOVA_SEEDFIELD ||
	(StarParticleRadiativeFeedback == TRUE &&
	 ParticleType[i] == PARTICLE_TYPE_STAR)) {

      // Check if it already exists (wasn't activated on the last
      // timestep, usually because of insufficient mass)
      exists = false;
      for (cstar = Stars; cstar; cstar = cstar->NextStar)
	if (cstar->Identifier == ParticleNumber[i]) {
	  cstar->SetLevel(level);
	  exists = true;
	  break;
	}

      if (!exists) {
	NewStar = new Star(this, i, level);

	/* If using an IMF for Pop III stars, assign the mass after
	   merging new (massless) star particles to avoid unnecessary
	   calls to the IMF. */

	if (ParticleType[i] == -PARTICLE_TYPE_SINGLE_STAR)
	  if (PopIIIInitialMassFunction == FALSE)
	    NewStar->AssignFinalMass(PopIIIStarMass);
	if (ParticleType[i] == -PARTICLE_TYPE_SIMPLE_SOURCE) 
	  NewStar->AssignFinalMass(PopIIIStarMass);
	InsertStarAfter(Stars, NewStar);
	NumberOfStars++;
      }

    }

  return SUCCESS;

}
