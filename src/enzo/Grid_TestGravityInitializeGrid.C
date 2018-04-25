/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GRAVITY TEST)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

/* Random number generator */
 
void mt_init(unsigned_int seed);
unsigned_long_int mt_random();
 
int grid::TestGravityInitializeGrid(float CentralDensity,
				    int NumberOfNewParticles,
				    int UseBaryons)
{
  /* declarations */
 
  int dim, i, size, field, vel;
  float phi, r, theta;
  int mt_random_seed = 123456789;
  int max_random = (1<<16);
  mt_init((unsigned_int) mt_random_seed);
 
  if (UseBaryons) {
 
    /* create fields */
 
    NumberOfBaryonFields = 0;
    FieldType[NumberOfBaryonFields++] = Density;
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;
    vel = NumberOfBaryonFields;
    FieldType[NumberOfBaryonFields++] = Velocity1;
    if (GridRank > 1)
      FieldType[NumberOfBaryonFields++] = Velocity2;
    if (GridRank > 2)
      FieldType[NumberOfBaryonFields++] = Velocity3;
  }
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Set particles. */
 
  if (NumberOfNewParticles > 0) {
 
    /* Set number of particles for this grid. */
 
    NumberOfParticles = NumberOfNewParticles;
 
    /* Allocate space. */
 
    this->AllocateNewParticles(NumberOfParticles);
 
    /* Create random particle positions and set velocities to zero. */
 
    for (dim = 0; dim < GridRank; dim++)
      for (i = 1; i < NumberOfParticles; i++) {
 
#ifdef UNUSED
 
	ParticlePosition[dim][i] = FLOAT(mt_random() % max_random)/FLOAT(max_random);
 
#endif /* UNUSED */
	
	ParticleVelocity[dim][i] = 0.0;
 
	ParticleMass[i] = tiny_number;
	
 
      }
 
    /* Set positions randomly distributed in log r from center. */
 
    float r1 = log10(0.2*(*CellWidth[0]));
    float r2 = log10(0.45*(DomainRightEdge[0] - DomainLeftEdge[0]));
    for (i = 0; i < NumberOfParticles; i++) {
 
      /* Compute random r, phi and theta. */
 
      r     = POW(10, (r2 - r1)*(float(mt_random() % max_random) /  float(max_random)) + r1);
      phi   = 2.0*pi*float(mt_random() % max_random) /  float(max_random);
      theta =     pi*float(mt_random() % max_random) /  float(max_random);
 
      /* Turn these into x/y/z. */
 
      if (GridRank == 1)
	ParticlePosition[0][i] = r*sign(phi - pi);
 
      if (GridRank == 2) {
	ParticlePosition[0][i] = r*cos(phi);
	ParticlePosition[1][i] = r*sin(phi);
      }
 
      if (GridRank == 3) {
	ParticlePosition[0][i] = r*sin(theta)*cos(phi);
	ParticlePosition[1][i] = r*sin(theta)*sin(phi);
	ParticlePosition[2][i] = r*cos(theta);
      }
 
      /* Shift center from 0,0,0 to the middle of the volume. */
 
      for (dim = 0; dim < GridRank; dim++)
	ParticlePosition[dim][i] +=
	  0.5*(DomainLeftEdge[dim] + DomainRightEdge[dim]);
 
      /* Set particle identifier. */
 
      ParticleNumber[i] = i;
      ParticleType[i]   = PARTICLE_TYPE_DARK_MATTER;
 
    }
 
    /* Set central particle (0). */
 
    if (!UseBaryons) {
 
      ParticleMass[0] = CentralDensity;
 
      for (dim = 0; dim < GridRank; dim++) {
	ParticlePosition[dim][0] = 0.5*(DomainLeftEdge[dim] +
		  		        DomainRightEdge[dim]);
	ParticleVelocity[dim][0] = 0.0;
      }
    }
 
  }
 
  /* If requested, set up the baryon field. */
 
  if (UseBaryons) {
 
    /* compute size of fields */
 
    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    /* allocate fields */
 
    this->AllocateGrids();
 
    /* set density, total energy */
 
    for (i = 0; i < size; i++) {
      BaryonField[0][i] = tiny_number;
      BaryonField[1][i] = tiny_number;
    }
 
    /* Set central density. */
 
    int Middle[MAX_DIMENSION], Size[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      Middle[dim] = max(GridDimension[dim]/2 - 1, 0);
      Size[dim]   = min(GridDimension[dim], 2);
    }
    float SpikeDensity = CentralDensity/float(Size[0]*Size[1]*Size[2]);
 
    for (int k = 0; k < Size[2]; k++)
      for (int j = 0; j < Size[1]; j++)
	for (int i = 0; i < Size[0]; i++)
	  BaryonField[0][
	    (k + Middle[2])*GridDimension[0]*GridDimension[1] +
	    (j + Middle[1])*GridDimension[0]                  +
	    (i + Middle[0])                                   ]
	  += SpikeDensity;
 
    /* set velocities */
 
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < size; i++)
	BaryonField[vel+dim][i] = 0.0;
 
    /* Set internal energy if necessary. */
 
    if (DualEnergyFormalism)
      for (i = 0; i < size; i++)
	BaryonField[2][i] = tiny_number;
 
  }
 
  return SUCCESS;
}
