/***********************************************************************
/
/  GRID CLASS (ADD A PARTICLE BELONGING TO THIS GRID FROM A LIST)
/
/  written by: Peng Wang
/  date:       January, 2009
/  modified1:  March, 2012 (JHW) -- removed boundary check and only 
/              moving one particle.
/
/  PURPOSE:
/
************************************************************************/

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

int grid::AddOneParticleFromList(ParticleEntry *List, const int place)
{

  /* Return if this doesn't involve us. */

  NumberOfParticles++;

  if (MyProcessorNumber != ProcessorNumber) return SUCCESS;

  /* Copy the old and new ones to a new list and
     get rid of the old one. */

  int i, n, dim;
  PINT *Number;
  int *Type;
  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION];
  float *Mass;
  float *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];

  Number = new PINT[NumberOfParticles];
  Type = new int[NumberOfParticles];
  Mass = new float[NumberOfParticles];
  for (dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[NumberOfParticles];
    Velocity[dim] = new float[NumberOfParticles];
  }
  for (i = 0; i < NumberOfParticleAttributes; i++)
    Attribute[i] = new float[NumberOfParticles];

  /* copy old particles to their new home. */

  for (i = 0; i < NumberOfParticles - 1; i++) {

    Mass[i]   = ParticleMass[i];
    Number[i] = ParticleNumber[i];
    Type[i]   = ParticleType[i];
    for (dim = 0; dim < GridRank; dim++) {
      Position[dim][i] = ParticlePosition[dim][i];
      Velocity[dim][i] = ParticleVelocity[dim][i];
    }
    for (n = 0; n < NumberOfParticleAttributes; n++)
      Attribute[n][i] = ParticleAttribute[n][i];
  
  }

  /* append new particles */

  i = NumberOfParticles - 1;
  Mass[i]   = List[place].Mass/pow(CellWidth[0][0],3);
  Number[i] = List[place].Number;
  Type[i]   = List[place].Type;
  for (dim = 0; dim < GridRank; dim++) {
    Position[dim][i] = List[place].Position[dim];
    Velocity[dim][i] = List[place].Velocity[dim];
  }
  for (n = 0; n < NumberOfParticleAttributes; n++)
    Attribute[n][i] = List[place].Attribute[n];

  /* Delete FromGrid's particles (now copied). */

  this->DeleteParticles();
  
  /* Copy new pointers into their correct position. */
  
  this->SetParticlePointers(Mass, Number, Type, Position, Velocity, 
			    Attribute);

  return SUCCESS;
}

