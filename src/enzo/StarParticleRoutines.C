/*-*-C++-*-*/
/***********************************************************************
/
/  ROUTINES FOR THE STAR PARTICLE STRUCTURE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE: Instead of restricting star particles to the typical 
/           particle attributes, this class gives more functionality 
/           to them.
/
************************************************************************/
#include <assert.h>
#include "typedefs.h"
#include "StarParticle.h"

StarParticle::StarParticle(void)
{
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    pos[dim] = vel[dim] = delta_vel[dim] = 0.0;
  accretion_rate = NULL;
  accretion_time = NULL;
  NextStar = NULL;
  CurrentGrid = NULL;
  Mass = FinalMass = DeltaMass = BirthTime = LifeTime = last_accretion_rate = 0.0;
  FeedbackFlag = Identifier = level = GridID = type = 0;
}

StarParticle::StarParticle(grid *_grid, int ParticleID)
{

  assert(ParticleID < _grid->NumberOfParticles);

  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    pos[dim] = _grid->ParticlePosition[dim][ParticleID];
    vel[dim] = _grid->ParticleVelocity[dim][ParticleID];
    delta_vel[dim] = 0.0;
  }
  accretion_rate = NULL;
  accretion_time = NULL;
  NextStar = NULL;
  CurrentGrid = _grid;
  DeltaMass = 0.0;
  last_accretion_rate = 0.0;
  level = 0;

  GridID = _grid->ID;
  type = _grid->ParticleType[ParticleID];
  Identifier = _grid->ParticleNumber[ParticleID];
  Mass = FinalMass = _grid->ParticleMass[ParticleID];
  BirthTime = _grid->ParticleAttribute[0][ParticleID];
  LifeTime = _grid->ParticleAttribute[1][ParticleID];
}

StarParticle::~StarParticle(void)
{
  if (accretion_rate != NULL)
    delete [] accretion_rate;
  if (accretion_time != NULL)
    delete [] accretion_time;
}

