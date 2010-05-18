/***********************************************************************
/
/  GRID CLASS (CLEAN UP THE PARTICLE LIST AFTER A MOVE)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
//
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::CleanUpMovedParticles()
{
 
  /* If there are no particles to clean up, we're done. */
 
  if (NumberOfParticles == 0)
    return SUCCESS;
 
  int i, j, n, dim, NumberOfParticlesRemaining = NumberOfParticles,
    *Type;
  PINT *Number;
  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION], *Mass,
        *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
 
  /* Count the number of unmoved particles. */
 
  for (i = 0; i < NumberOfParticles; i++)
    if (ParticleMass[i] == FLOAT_UNDEFINED)
      NumberOfParticlesRemaining--;
 
  if (NumberOfParticlesRemaining == NumberOfParticles)
    return SUCCESS;
 
  /* Get rid of moved particles in FromGrid. */
 
  if (NumberOfParticlesRemaining == 0)
    this->DeleteParticles();
 
  /* If there are unmoved particles left, then copy them to a new list and
     get rid of the old one. */
 
  else {
 
    /* Allocate space for the new set of particles. */
 
    Mass = new float[NumberOfParticlesRemaining];
    Number = new PINT[NumberOfParticlesRemaining];
    Type = new int[NumberOfParticlesRemaining];
    for (dim = 0; dim < GridRank; dim++) {
      Position[dim] = new FLOAT[NumberOfParticlesRemaining];
      Velocity[dim] = new float[NumberOfParticlesRemaining];
    }
    for (i = 0; i < NumberOfParticleAttributes; i++)
      Attribute[i] = new float[NumberOfParticlesRemaining];
 
    /* Copy unmoved particles to their new home. */
 
    j = 0;
    for (i = 0; i < NumberOfParticles; i++)
 
      if (ParticleMass[i] != FLOAT_UNDEFINED) {
 
	Mass  [j] = ParticleMass  [i];
	Number[j] = ParticleNumber[i];
	Type  [j] = ParticleType  [i];
	for (dim = 0; dim < GridRank; dim++) {
	  Position[dim][j] = ParticlePosition[dim][i];
	  Velocity[dim][j] = ParticleVelocity[dim][i];
	}
	for (n = 0; n < NumberOfParticleAttributes; n++)
	  Attribute[n][j] = ParticleAttribute[n][i];
 
	j++;   // increment moved particle counter
 
      }
 
    /* Delete FromGrid's particles (now copied). */
 
    this->DeleteParticles();
 
    /* Copy new pointers into their correct position. */
 
    this->SetParticlePointers(Mass, Number, Type, Position, Velocity,
			      Attribute);
 
  }
 
  /* Set the new number of particles. */
 
  NumberOfParticles = NumberOfParticlesRemaining;
 
  return SUCCESS;
}
