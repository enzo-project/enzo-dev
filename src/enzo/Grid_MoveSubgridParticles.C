/***********************************************************************
/
/  GRID CLASS (MOVE APPROPRIATE PARTICLES FROM SPECIFIED GRID TO THIS GRID)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness, Jan 2002
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
 
int grid::MoveSubgridParticles(grid* FromGrid,
                               int *Counter,
                               PINT *X_Number,
                               int *X_Type,
                               float *X_Mass,
                               FLOAT *X_Position[],
                               float *X_Velocity[],
                               float *X_Attribute[])
{
 
int start;
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* Error check. */
 
  if (ProcessorNumber != FromGrid->ProcessorNumber) {
    ENZO_FAIL("This routine not parallelized.\n");
  }
 
  /* If there are no particles to move, we're done. */
 
  if (FromGrid->NumberOfParticles == 0)
    return SUCCESS;
 
  int i, j, k, dim;
 
  /* To begin, set all particles to move. */
 
  int *MoveParticle = new int[FromGrid->NumberOfParticles];
 
  for (i = 0; i < FromGrid->NumberOfParticles; i++)
    if (FromGrid->ParticleMass[i] == FLOAT_UNDEFINED)
      MoveParticle[i] = FALSE;
    else
      MoveParticle[i] = TRUE;
 
  /* If a particle is outside of this grid (subgrid of FromGrid), unmark it. */
 
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < FromGrid->NumberOfParticles; i++)
      if (FromGrid->ParticlePosition[dim][i] < GridLeftEdge[dim] ||
	  FromGrid->ParticlePosition[dim][i] > GridRightEdge[dim] )
	MoveParticle[i] = FALSE;
 
  /* Compute the number of particles left. */
 
  int TotalNumberOfParticles = 0;
  for (i = 0; i < FromGrid->NumberOfParticles; i++)
    TotalNumberOfParticles += MoveParticle[i];
 
  /* If there are no particles to move, clean up and exit. */
 
  if (TotalNumberOfParticles == 0) {
    delete MoveParticle;
    return SUCCESS;
  }
 
  /* Compute the number of particles left in the old grid (needed later). */
 
  int FromNumberOfParticles = FromGrid->NumberOfParticles -
                              TotalNumberOfParticles;
 
  /* Debugging info. */
 
  if (debug)
    printf("MoveSubgridParticles: %"ISYM" particles (after: Top = %"ISYM", Sub = %"ISYM").\n",
	   TotalNumberOfParticles, FromNumberOfParticles,
	   TotalNumberOfParticles + NumberOfParticles);
 
  /* Add in this grid's particles. */
 
  TotalNumberOfParticles += NumberOfParticles;
 
  /* Compute the increase in mass for particles moving to the subgrid. */
 
  float RefinementFactors[MAX_DIMENSION];
  FromGrid->ComputeRefinementFactorsFloat(this, RefinementFactors);
  float MassIncrease = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    MassIncrease *= RefinementFactors[dim];
 
  /* (1) Move Particles from FromGrid to this grid. */
 
  /* Allocate space for the particles. */
 
/* KILL KILL KILL
 
  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION], *Mass,
        *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  int   *Number, *Type;
 
  Mass = new float[TotalNumberOfParticles];
  Number = new int[TotalNumberOfParticles];
  Type = new int[TotalNumberOfParticles];
  for (dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[TotalNumberOfParticles];
    Velocity[dim] = new float[TotalNumberOfParticles];
  }
  for (i = 0; i < NumberOfParticleAttributes; i++)
    Attribute[i] = new float[TotalNumberOfParticles];
 
*/
 
  /* Copy this grid's particles to the new space. */
 
  /*
  for (i = 0; i < NumberOfParticles; i++) {
    Mass  [i] = ParticleMass  [i];
    Number[i] = ParticleNumber[i];
    Type  [i] = ParticleType  [i];
  }
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < NumberOfParticles; i++) {
      Position[dim][i] = ParticlePosition[dim][i];
      Velocity[dim][i] = ParticleVelocity[dim][i];
    }
  for (j = 0; j < NumberOfParticleAttributes; j++)
    for (i = 0; i < NumberOfParticles; i++)
      Attribute[j][i] = ParticleAttribute[j][i];
  */
 
  /* Copy FromGrid's particles to new space (starting at NumberOfParticles). */
 
  start = *Counter;
  printf("Counter %"ISYM"\n",start);
 
  j = 0;
  for (i = 0; i < FromGrid->NumberOfParticles; i++)
 
    if (MoveParticle[i] == TRUE) {

 
//      Mass[j + NumberOfParticles] = (FromGrid->ParticleMass[i]) * MassIncrease;
//      Number[j + NumberOfParticles] = FromGrid->ParticleNumber[i];
//      Type[j + NumberOfParticles] = FromGrid->ParticleType[i];
 
//      for (dim = 0; dim < GridRank; dim++) {
//	Position[dim][j + NumberOfParticles] =
//	  FromGrid->ParticlePosition[dim][i];
//	Velocity[dim][j + NumberOfParticles] =
//	  FromGrid->ParticleVelocity[dim][i];
//      }
//      for (k = 0; k < NumberOfParticleAttributes; k++)
//	Attribute[k][j + NumberOfParticles] =
//	  FromGrid->ParticleAttribute[k][i];
 
      X_Number[start] = FromGrid->ParticleNumber[i];
      X_Type[start] = FromGrid->ParticleType[i];
      X_Mass[start] = (FromGrid->ParticleMass[i]) * MassIncrease;
      for (dim = 0; dim < GridRank; dim++) {
        X_Position[dim][start] = FromGrid->ParticlePosition[dim][i];
        X_Velocity[dim][start] = FromGrid->ParticleVelocity[dim][i];
      }
 
      for (k = 0; k < NumberOfParticleAttributes; k++) {
        X_Attribute[k][start] = FromGrid->ParticleAttribute[k][i];
      }
 
      j++;   // increment moved particle counter
      start++;
 
      FromGrid->ParticleMass[i] = FLOAT_UNDEFINED; // erase old one
 
    }
 
  printf("Counter %"ISYM"\n",start);
  *Counter = start;
 
  /* Delete this grid's particles (now copied). */
 
  this->DeleteParticles();
 
  /* Copy new pointers into their correct position. */
 
//  this->SetParticlePointers(Mass, Number, Type, Position, Velocity, Attribute);
 
  /* Set this's grid's new NumberOfParticles. */
 
  NumberOfParticles = TotalNumberOfParticles;
 
  /* Clean up. */
 
  delete MoveParticle;
 
  return SUCCESS;
}
