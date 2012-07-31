/***********************************************************************
/
/  GRID CLASS (ADD PARTICLES BELONGING TO THIS GRID FROM A LIST)
/
/  written by: Peng Wang
/  date:       January, 2009
/  modified1:
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

int grid::AddParticlesFromList(ParticleEntry *List, const int &Size,
			       int *AddedNewParticleNumber)
{

  /* Return if this doesn't involve us. */

  if (MyProcessorNumber != ProcessorNumber) return SUCCESS;

  //AddedNewParticleNumber = 0;
  if (Size < 1) return SUCCESS;

  /* Count the number of unmoved particles. */

  int *NewIndex = new int[Size];
  //printf("Grid: x: %g %g y: %g %g z: %g %g\n", GridLeftEdge[0], GridRightEdge[0],
  // GridLeftEdge[1],GridRightEdge[1],
  // GridLeftEdge[2],GridRightEdge[2]);
  int Count = 0;
  for (int i = 0; i < Size; i++) {
    //printf("%d: %g %g %g\n", List[i].Number, List[i].Position[0], List[i].Position[1], List[i].Position[2]);
    if (List[i].Position[0] > GridLeftEdge[0] &&
	List[i].Position[0] < GridRightEdge[0] &&
	List[i].Position[1] > GridLeftEdge[1] &&
	List[i].Position[1] < GridRightEdge[1] &&
	List[i].Position[2] > GridLeftEdge[2] &&
	List[i].Position[2] < GridRightEdge[2]) {
      NewIndex[Count++] = i;
      AddedNewParticleNumber[i] = 1;
    }
  }

  if (!Count) {
    delete [] NewIndex;
    return SUCCESS;
  }

  NumberOfParticles += Count;

  /* Copy the old and new ones to a new list and
     get rid of the old one. */

  PINT *Number;
  int *Type;
  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION];
  float *Mass;
  float *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];

  Number = new PINT[NumberOfParticles];
  Type = new int[NumberOfParticles];
  Mass = new float[NumberOfParticles];
  for (int dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[NumberOfParticles];
    Velocity[dim] = new float[NumberOfParticles];
  }
  for (int i = 0; i < NumberOfParticleAttributes; i++)
    Attribute[i] = new float[NumberOfParticles];

  /* copy old particles to their new home. */

  for (int i = 0; i < NumberOfParticles - Count; i++) {

    Mass[i]   = ParticleMass[i];
    Number[i] = ParticleNumber[i];
    Type[i]   = ParticleType[i];
    for (int dim = 0; dim < GridRank; dim++) {
      Position[dim][i] = ParticlePosition[dim][i];
      Velocity[dim][i] = ParticleVelocity[dim][i];
    }
    for (int n = 0; n < NumberOfParticleAttributes; n++)
      Attribute[n][i] = ParticleAttribute[n][i];
  
  }

  /* append new particles */
  
  for (int i = NumberOfParticles - Count, j = 0; i < NumberOfParticles; i++, j++) {
    Mass[i]   = List[NewIndex[j]].Mass/pow(CellWidth[0][0],3);
    Number[i] = List[NewIndex[j]].Number;
    Type[i]   = List[NewIndex[j]].Type;
    for (int dim = 0; dim < GridRank; dim++) {
      Position[dim][i] = List[NewIndex[j]].Position[dim];
      Velocity[dim][i] = List[NewIndex[j]].Velocity[dim];
    }
    for (int n = 0; n < NumberOfParticleAttributes; n++)
      Attribute[n][i] = List[NewIndex[j]].Attribute[n];
  }

  /* Delete FromGrid's particles (now copied). */

  this->DeleteParticles();
  
  /* Copy new pointers into their correct position. */
  
  this->SetParticlePointers(Mass, Number, Type, Position, Velocity, 
			    Attribute);

  delete [] NewIndex;

  return SUCCESS;
}

/************************************************************************/

int grid::CheckGridBoundaries(FLOAT *Position)
{
  //printf("x: %g %g y: %g %g z: %g %g\n", GridLeftEdge[0], GridRightEdge[0],
  // GridLeftEdge[1],GridRightEdge[1],
  // GridLeftEdge[2],GridRightEdge[2]);

  if (Position[0] >= GridLeftEdge[0] &&
      Position[0] <= GridRightEdge[0] &&
      Position[1] >= GridLeftEdge[1] &&
      Position[1] <= GridRightEdge[1] &&
      Position[2] >= GridLeftEdge[2] &&
      Position[2] <= GridRightEdge[2])
    return 1;

  return 0;

}
