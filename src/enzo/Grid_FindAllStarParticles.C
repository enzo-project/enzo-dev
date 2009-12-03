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

#define RESET_BH_LIFETIMES

int grid::FindAllStarParticles(int level)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0)
    return SUCCESS;

  int i, StarType;
  Star *NewStar;
  
  /* Read only active star particles.  Unborn stars will be read later
     in grid::FindNewStarParticles. */
  
  for (i = 0; i < NumberOfParticles; i++) {
    //StarType = abs(ParticleType[i]);
    if (ParticleType[i] == PARTICLE_TYPE_SINGLE_STAR ||
	ParticleType[i] == PARTICLE_TYPE_BLACK_HOLE ||
	ParticleType[i] == PARTICLE_TYPE_CLUSTER ||
        ParticleType[i] == PARTICLE_TYPE_COLOR_STAR ||
	ParticleType[i] == PARTICLE_TYPE_MBH) {

#ifdef RESET_BH_LIFETIMES // Make BH lifetimes "infinite"
      if (ParticleType[i] == PARTICLE_TYPE_BLACK_HOLE &&
	  ParticleAttribute[1][i] < 1)
	ParticleAttribute[1][i] = huge_number;
#endif /* RESET_BH_LIFETIMES */

      NewStar = new Star(this, i, level);
      InsertStarAfter(Stars, NewStar);
      NumberOfStars++;

      /* For MBHFeedback=2 (FeedbackFlag=MBH_JETS), you need the angular 
	 momentum;  if no file to read in, assume zero angular momentum 
	 accreted so far.  -Ji-hoon Kim, Nov.2009 */

      if(MBHFeedback == 2 && ParticleType[i] == PARTICLE_TYPE_MBH) {
	NewStar->AssignAccretedAngularMomentum();
      	printf("MBH particle info (for MBHFedback=2, check angular momentum): \n"); 
	NewStar->PrintInfo();
      }
 
    } // ENDIF star
  } // ENDFOR particles

  return SUCCESS;

}
