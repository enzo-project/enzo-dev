/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A SPHERE GRAVITY TEST)
/
/  written by: Greg Bryan
/  date:       September, 1995
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
 
 
int grid::TestGravitySphereInitializeGrid(float SphereInteriorDensity,
					  float SphereExteriorDensity,
					  float SphereRadius,
					  int   SphereType,
					  int   UseBaryons,
					  FLOAT SphereCenter[])
{
  /* declarations */
 
  int dim, i, j, k, size, vel, field;
  float phi, r, theta, pi = 3.14159;
 
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

    if (WritePotential)
      FieldType[NumberOfBaryonFields++] = GravPotential;

  }
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Set particles. */
 
  if (NumberOfParticles > 0) {
 
    /* Allocate space. */
 
    this->AllocateNewParticles(NumberOfParticles);
 
    /* Create random particle positions and set velocities to zero. */
 
    for (dim = 0; dim < GridRank; dim++)
      for (i = 1; i < NumberOfParticles; i++) {
 
#ifdef UNUSED
 
	ParticlePosition[dim][i] = FLOAT(rand())/FLOAT(RAND_MAX);
 
#endif /* UNUSED */
	
	ParticleVelocity[dim][i] = 0.0;
 
	ParticleMass[i] = tiny_number;
 
      }
 
    /* Set positions randomly distributed in log r from center. */
 
    float r1 = log10(0.2*CellWidth[0][0]);
    float r2 = log10(0.45*(DomainRightEdge[0] - DomainLeftEdge[0]));
    for (i = 0; i < NumberOfParticles; i++) {
 
      /* Compute random r, phi and theta. */
 
      r     = POW(10, (r2 - r1)*float(rand())/float(RAND_MAX) + r1);
      phi   = 2.0*pi*float(rand())/float(RAND_MAX);
      theta =     pi*float(rand())/float(RAND_MAX);
 
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
	ParticlePosition[dim][i] += SphereCenter[dim];
 
      /* Set particle identifier. */
 
      ParticleNumber[i] = i;
      ParticleType[i]   = PARTICLE_TYPE_DARK_MATTER;
 
    }
 
  }
 
  /* If requested, set up the baryon field. */
 
  if (UseBaryons) {
 
    /* compute size of fields */
 
    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    /* allocate fields */
 
    for (field = 0; field < NumberOfBaryonFields; field++)
      if (BaryonField[field] == NULL)
	BaryonField[field] = new float[size];
 
    /* Set densities. */
 
    float density, r, x, y = 0, z = 0;
 
    for (k = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++)
	for (i = 0; i < GridDimension[0]; i++) {
 
	  /* Compute position */
 
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  if (GridRank > 1)
	    y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  if (GridRank > 2)
	    z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
 
	  /* Find distance from center. */
 
	  r = sqrt((x - SphereCenter[0]) * (x - SphereCenter[0]) +
	           (y - SphereCenter[1]) * (y - SphereCenter[1]) +
	           (z - SphereCenter[2]) * (z - SphereCenter[2]) );
 
	  /* set density */
 
	  if (r < max(SphereRadius, CellWidth[0][0])) {
	    if (SphereType == 0)
	      density = SphereInteriorDensity;
	    if (SphereType == 1)
	      density = SphereInteriorDensity*POW(r/SphereRadius, -2);
	    if (SphereType == 2)
	      density = SphereInteriorDensity*POW(r/SphereRadius,
						  float(-2.25));
	  }
	  else
	    density = SphereExteriorDensity;
 
	  BaryonField[0][k*GridDimension[0]*GridDimension[1] +
	                 j*GridDimension[0]                  + i] =
	    density;
	}
 
    /* Set velocities */
 
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < size; i++)
	BaryonField[vel+dim][i] = 0.0;
 
    /* set total energy */
 
    for (i = 0; i < size; i++)
      BaryonField[1][i] = tiny_number;
 
    /* Set internal energy if necessary. */
 
    if (DualEnergyFormalism)
      for (i = 0; i < size; i++)
	BaryonField[2][i] = tiny_number;
 
  }
 
  return SUCCESS;
}
