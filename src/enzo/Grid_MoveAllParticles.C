/***********************************************************************
/
/  GRID CLASS (MOVE ALL PARTICLES FROM SPECIFIED GRID TO THIS GRID)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  May, 2009 by John Wise
/                Move particles to "empty" grid on the processor local
/                to the subgrid.  This distributes memory usage during
/                RebuildHierarchy in simulations with nested grids.
/
/  PURPOSE:
/
/    NOTE: We assume all the from grids are at the same level!
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

int grid::MoveAllParticles(int NumberOfGrids, grid* FromGrid[])
{

  if (NumberOfGrids < 1) {
    ENZO_VFAIL("NumberOfGrids(%"ISYM") must be > 0.\n", NumberOfGrids)
  }
 
  /* Determine total number of local particles. */

  int NumberOfSubgridParticles = 0;
  int TotalNumberOfParticles = NumberOfParticles;
  int i, j, grid, dim, *Type;
  PINT *Number;
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (MyProcessorNumber == FromGrid[grid]->ProcessorNumber)
      NumberOfSubgridParticles += FromGrid[grid]->NumberOfParticles;
  if (NumberOfSubgridParticles == 0) 
    return SUCCESS;
  
  TotalNumberOfParticles += NumberOfSubgridParticles;
 
  /* Debugging info. */

  if (debug1) printf("MoveAllParticles: %"ISYM" (before: ThisGrid = %"ISYM").\n",
		     TotalNumberOfParticles, NumberOfParticles);
 
  /* Allocate space for the particles. */
 
  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION], *Mass,
        *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
 
  Mass = new float[TotalNumberOfParticles];
  Number = new PINT[TotalNumberOfParticles];
  Type = new int[TotalNumberOfParticles];
  for (int dim = 0; dim < GridRank; dim++) {
    Position[dim] = new FLOAT[TotalNumberOfParticles];
    Velocity[dim] = new float[TotalNumberOfParticles];
  }
  for (int i = 0; i < NumberOfParticleAttributes; i++)
    Attribute[i] = new float[TotalNumberOfParticles];
  
  if (Velocity[GridRank-1] == NULL) {
    ENZO_FAIL("malloc error (out of memory?)\n");
  }
 
  /* Compute the decrease in mass for particles moving to this grid
     (We assume here all grids are from the same level). */
 
  float RefinementFactors[MAX_DIMENSION];
  this->ComputeRefinementFactorsFloat(FromGrid[0], RefinementFactors);
  float MassDecrease = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    MassDecrease *= RefinementFactors[dim];
  MassDecrease = 1.0/MassDecrease;
 
  /* Copy this grid's particles to the new space. */
 
  for (i = 0; i < NumberOfParticles; i++) {
    Mass[i]   = ParticleMass[i];
    Number[i] = ParticleNumber[i];
    Type[i]   = ParticleType[i];
  }
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < NumberOfParticles; i++) {
      Position[dim][i] = ParticlePosition[dim][i];
      Velocity[dim][i] = ParticleVelocity[dim][i];
    }
  for (j = 0; j < NumberOfParticleAttributes; j++)
    for (i = 0; i < NumberOfParticles; i++)
      Attribute[j][i] = ParticleAttribute[j][i];
 
  /* Delete this grid's particles (now copied). */
 
  this->DeleteParticles();
 
  /* Copy new pointers into their correct position. */
 
  this->SetParticlePointers(Mass, Number, Type, Position, Velocity,
			    Attribute);
 
  /* Copy FromGrids' particles to new space on local "fake" grid. */
 
  int Index = NumberOfParticles;
  for (grid = 0; grid < NumberOfGrids; grid++) {

    for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++) {
      Mass[Index+i] = FromGrid[grid]->ParticleMass[i] * MassDecrease;
      Number[Index+i] = FromGrid[grid]->ParticleNumber[i];
      Type[Index+i] = FromGrid[grid]->ParticleType[i];
    }
    
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++) {
	Position[dim][Index+i] = FromGrid[grid]->ParticlePosition[dim][i];
	Velocity[dim][Index+i] = FromGrid[grid]->ParticleVelocity[dim][i];
      }

    for (j = 0; j < NumberOfParticleAttributes; j++)
      for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++)
	Attribute[j][Index+i] = FromGrid[grid]->ParticleAttribute[j][i];
    
    Index += FromGrid[grid]->NumberOfParticles;

  } // ENDFOR grids 
  
  /* Set new number of particles in this grid. */
 
  NumberOfParticles = TotalNumberOfParticles;
 
  /* Delete FromGrid's particles (and set number of particles to zero). */
 
  for (grid = 0; grid < NumberOfGrids; grid++) {
    FromGrid[grid]->NumberOfParticles = 0;
    FromGrid[grid]->DeleteParticles();
  }
 
  return SUCCESS;
}











/* Old version of MoveAllParticles : 
   might just work for less computationally intensive runs */

int grid::MoveAllParticlesOld(int NumberOfGrids, grid* FromGrid[])  
{

  if (NumberOfGrids < 1) {
    fprintf(stderr, "NumberOfGrids(%d) must be > 0.\n", NumberOfGrids);
    return FAIL;
  }

  /* Determine total number of particles. */

  int TotalNumberOfParticles = NumberOfParticles;
  int i, j, grid, dim, *Type;
  PINT *Number;

  for (grid = 0; grid < NumberOfGrids; grid++)
    TotalNumberOfParticles += FromGrid[grid]->NumberOfParticles;
  if (TotalNumberOfParticles == 0)
    return SUCCESS;

  /* Debugging info. */

  if (debug) printf("MoveAllParticles: %d (before: ThisGrid = %d).\n",
		    TotalNumberOfParticles, NumberOfParticles);

  /* Allocate space for the particles. */

  FLOAT *Position[MAX_DIMENSION];
  float *Velocity[MAX_DIMENSION], *Mass,
        *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];

  if (MyProcessorNumber == ProcessorNumber) {
     Mass = new float[TotalNumberOfParticles];
     Number = new PINT[TotalNumberOfParticles]; 
     Type = new int[TotalNumberOfParticles];
     for (int dim = 0; dim < GridRank; dim++) {
       Position[dim] = new FLOAT[TotalNumberOfParticles];
       Velocity[dim] = new float[TotalNumberOfParticles];
     }
     for (int i = 0; i < NumberOfParticleAttributes; i++)
       Attribute[i] = new float[TotalNumberOfParticles];

     if (Velocity[GridRank-1] == NULL) {
       fprintf(stderr, "malloc error (out of memory?)\n");
       return FAIL;
     }
  }

  /* Compute the decrease in mass for particles moving to this grid
     (We assume here all grids are from the same level). */

  float RefinementFactors[MAX_DIMENSION];
  this->ComputeRefinementFactorsFloat(FromGrid[0], RefinementFactors);
  float MassDecrease = 1.0;
  for (dim = 0; dim < GridRank; dim++)
    MassDecrease *= RefinementFactors[dim];
  MassDecrease = 1.0/MassDecrease;

  /* Copy this grid's particles to the new space. */

   if (MyProcessorNumber == ProcessorNumber) {
     for (i = 0; i < NumberOfParticles; i++) {
       Mass[i]   = ParticleMass[i];
       Number[i] = ParticleNumber[i];
       Type[i]   = ParticleType[i];
     }
     for (dim = 0; dim < GridRank; dim++)
       for (i = 0; i < NumberOfParticles; i++) {
	 Position[dim][i] = ParticlePosition[dim][i];
	 Velocity[dim][i] = ParticleVelocity[dim][i];
       }
     for (j = 0; j < NumberOfParticleAttributes; j++)
       for (i = 0; i < NumberOfParticles; i++)
	 Attribute[j][i] = ParticleAttribute[j][i];
   }

  /* Delete this grid's particles (now copied). */  

  if (MyProcessorNumber == ProcessorNumber) {
    this->DeleteParticles();

    /* Copy new pointers into their correct position. */

    this->SetParticlePointers(Mass, Number, Type, Position, Velocity,  
			      Attribute);
  }

  /* Copy FromGrids' particles to new space (starting at NumberOfParticles). */

  int Index = NumberOfParticles;
  for (grid = 0; grid < NumberOfGrids; grid++) {

   /* If on the same processor, just copy. */

    if (MyProcessorNumber == ProcessorNumber &&
        MyProcessorNumber == FromGrid[grid]->ProcessorNumber) {

      //      fprintf(stderr, "P(%d) copying %d particles\n", MyProcessorNumber,
      //	     FromGrid[grid]->NumberOfParticles);

      for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++) {
	Mass[Index+i] = FromGrid[grid]->ParticleMass[i] * MassDecrease;
	Number[Index+i] = FromGrid[grid]->ParticleNumber[i];
	Type[Index+i] = FromGrid[grid]->ParticleType[i];
      }

      for (dim = 0; dim < GridRank; dim++)
	for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++) {
	  Position[dim][Index+i] = FromGrid[grid]->ParticlePosition[dim][i];
	  Velocity[dim][Index+i] = FromGrid[grid]->ParticleVelocity[dim][i];
	}
      for (j = 0; j < NumberOfParticleAttributes; j++)
	for (i = 0; i < FromGrid[grid]->NumberOfParticles; i++)
	  Attribute[j][Index+i] = FromGrid[grid]->ParticleAttribute[j][i];
    }

    /* Otherwise, communicate. */

    else {
      if (MyProcessorNumber == ProcessorNumber ||
          MyProcessorNumber == FromGrid[grid]->ProcessorNumber)
	if (FromGrid[grid]->CommunicationSendParticles(this, ProcessorNumber,
              0, FromGrid[grid]->NumberOfParticles, Index) == FAIL) {
	  fprintf(stderr, "Error in grid->CommunicationSendParticles.\n");
	  return FAIL;
        }

      /* Change mass, as required. */

      if (MyProcessorNumber == ProcessorNumber)
	for (i = Index; i < Index+FromGrid[grid]->NumberOfParticles; i++)
	  Mass[i] *= MassDecrease;

    }

    Index += FromGrid[grid]->NumberOfParticles;

  } // end: loop over grids.

  NumberOfParticles = TotalNumberOfParticles; 

  /* Delete FromGrid's particles (and set number of particles to zero). */

  for (grid = 0; grid < NumberOfGrids; grid++) {
    FromGrid[grid]->NumberOfParticles = 0;
    if (MyProcessorNumber == FromGrid[grid]->ProcessorNumber)

      FromGrid[grid]->DeleteParticles();
  }

  return SUCCESS;
}
