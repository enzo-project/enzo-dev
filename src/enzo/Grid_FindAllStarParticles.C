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
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

#define RESET_BH_LIFETIMES
#define NO_RESET_MBH_MASS  //#####

int grid::FindAllStarParticles(int level)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0)
    return SUCCESS;

  int i, StarType;
  float LifetimeFactor;
  Star *NewStar;
  
  /* Read only active star particles.  Unborn stars will be read later
     in grid::FindNewStarParticles. */
  
  for (i = 0; i < NumberOfParticles; i++) {
    //StarType = abs(ParticleType[i]);
    if (ParticleType[i] == PARTICLE_TYPE_SINGLE_STAR ||
	ParticleType[i] == PARTICLE_TYPE_BLACK_HOLE ||
	(ParticleType[i] == PARTICLE_TYPE_STAR && UseSupernovaSeedFieldSourceTerms) ||
	ParticleType[i] == PARTICLE_TYPE_CLUSTER ||
        ParticleType[i] == PARTICLE_TYPE_COLOR_STAR ||
	ParticleType[i] == PARTICLE_TYPE_MBH ||
	ParticleType[i] == PARTICLE_TYPE_SIMPLE_SOURCE ||
        ParticleType[i] == PARTICLE_TYPE_SUPERNOVA_SEEDFIELD ||
	(StarParticleRadiativeFeedback == TRUE &&
	 ParticleType[i] == PARTICLE_TYPE_STAR)) {

      if (ParticleType[i] == PARTICLE_TYPE_STAR)
	LifetimeFactor = 12.0;
      else
	LifetimeFactor = 1.0;
      
      if (this->Time >= ParticleAttribute[0][i] &&
	  this->Time <= ParticleAttribute[0][i] + 
	  LifetimeFactor * ParticleAttribute[1][i]) {
	  
#ifdef RESET_BH_LIFETIMES // Make BH lifetimes "infinite"
	if (ParticleType[i] == PARTICLE_TYPE_BLACK_HOLE &&
	    ParticleAttribute[1][i] < 1)
	  ParticleAttribute[1][i] = huge_number;
#endif /* RESET_BH_LIFETIMES */

#ifdef RESET_MBH_MASS // Edit MBH Mass; only for test purpose
	const double Msun = 1.989e33;
	float MassConversion, DensityUnits, LengthUnits, TemperatureUnits, 
	  TimeUnits, VelocityUnits;
	GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, this->Time);

	double dx = LengthUnits * this->CellWidth[0][0];
	MassConversion = (float) (dx*dx*dx * double(DensityUnits) / Msun);

	if (ParticleType[i] == PARTICLE_TYPE_MBH)
	  ParticleMass[i] = 1.0e5 / MassConversion;
#endif /* RESET_MBH_MASS */  

	NewStar = new Star(this, i, level);
	InsertStarAfter(Stars, NewStar);
	NumberOfStars++;

	/* For MBHFeedback = 2 to 5 (FeedbackFlag=MBH_JETS), you need
	   the angular momentum; if no file to read in, assume zero
	   angular momentum accreted so far.  -Ji-hoon Kim, Nov.2009 */

	if((MBHFeedback >= 2 && MBHFeedback <= 5) && 
	   ParticleType[i] == PARTICLE_TYPE_MBH) {
	  NewStar->AssignAccretedAngularMomentum();
	  printf("MBH particle info (for MBHFeedback=2 to 5, check angular momentum): \n"); 
	  NewStar->PrintInfo();
	}
      } // ENDIF during lifetime
    } // ENDIF star
  } // ENDFOR particles
  
  return SUCCESS;

}
