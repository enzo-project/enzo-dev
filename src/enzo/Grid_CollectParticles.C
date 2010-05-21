/***********************************************************************
/
/  GRID CLASS (COLLECT PARTICLES INTO ONE PROCESSOR)
/
/  written by: John Wise
/  date:       May, 2009
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
 
 
int grid::CollectParticles(int GridNum, int* &NumberToMove, 
			   int &StartIndex, int &EndIndex, 
			   particle_data *&List, int CopyDirection)
{
 
  /* Declarations. */

  int i, j, dim, n1, grid, proc;

  /* ----------------------------------------------------------------- */
  /* Copy particle out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If there are no particles to move, we're done. */

    if (NumberOfParticles == 0)
      return SUCCESS;

    /* If this is the correct processor, no copy-outs required. */

    if (MyProcessorNumber == ProcessorNumber)
      return SUCCESS;

    /* Add to the particle count to move */

    // NumberOfParticles is still the number of local particles, not
    // the actual total!
    NumberToMove[ProcessorNumber] += NumberOfParticles;
 
    /* Move and delete particles */

    n1 = StartIndex;

    for (i = 0; i < NumberOfParticles; i++) {
      for (dim = 0; dim < GridRank; dim++) {
	List[n1].pos[dim] = ParticlePosition[dim][i];
	List[n1].vel[dim] = ParticleVelocity[dim][i];
      }
      List[n1].mass = ParticleMass[i];
      List[n1].id = ParticleNumber[i];
      List[n1].type = ParticleType[i];
      for (j = 0; j < NumberOfParticleAttributes; j++)
	List[n1].attribute[j] = ParticleAttribute[j][i];
      List[n1].grid = GridNum;
      List[n1].proc = ProcessorNumber;
      //ParticleMass[i] = FLOAT_UNDEFINED;
      n1++;
    } // ENDFOR particles

    StartIndex = n1;
    this->DeleteParticles();

  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy particle back into grid. */
 
  else {

    if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* Count up total number. */
 
    int TotalNumberOfParticles;
    int NumberOfNewParticles = EndIndex - StartIndex;

    TotalNumberOfParticles = NumberOfParticles + NumberOfNewParticles;
//    for (i = 0; i < NumberOfParticles; i++)
//      if (ParticleMass[i] == FLOAT_UNDEFINED)
//	TotalNumberOfParticles--;
 
    /* Allocate space for the particles. */

    FLOAT *Position[MAX_DIMENSION];
    float *Velocity[MAX_DIMENSION], *Mass,
          *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    int *Type;
    PINT *Number;
 
    if (TotalNumberOfParticles > 0) {

      Mass = new float[TotalNumberOfParticles];
      Number = new PINT[TotalNumberOfParticles];
      Type = new int[TotalNumberOfParticles];
      for (dim = 0; dim < GridRank; dim++) {
	Position[dim] = new FLOAT[TotalNumberOfParticles];
	Velocity[dim] = new float[TotalNumberOfParticles];
      }
      for (i = 0; i < NumberOfParticleAttributes; i++)
	Attribute[i] = new float[TotalNumberOfParticles];

      if (Velocity[GridRank-1] == NULL && TotalNumberOfParticles != 0) {
	ENZO_FAIL("malloc error (out of memory?)\n");
      }
 
      /* Copy this grid's particles to the new space. */

      int n;

      for (i = 0; i < NumberOfParticles; i++) {
	Mass[i] = ParticleMass[i];
	Number[i] = ParticleNumber[i];
	Type[i] = ParticleType[i];
      }

      for (dim = 0; dim < GridRank; dim++)
	for (i = 0; i < NumberOfParticles; i++) {
	  Position[dim][i] = ParticlePosition[dim][i];
	  Velocity[dim][i] = ParticleVelocity[dim][i];
	}
	
      for (j = 0; j < NumberOfParticleAttributes; j++)
	for (i = 0; i < NumberOfParticles; i++)
	  Attribute[j][i] = ParticleAttribute[j][i];
 
      /* Copy new particles */

      n = NumberOfParticles;
      for (i = StartIndex; i < EndIndex; i++) {
	Mass[n] = List[i].mass;
	Number[n] = List[i].id;
	Type[n] = List[i].type;
	n++;
      }

      for (dim = 0; dim < GridRank; dim++) {
	n = NumberOfParticles;
	for (i = StartIndex; i < EndIndex; i++) {
	  Position[dim][n] = List[i].pos[dim];
	  Velocity[dim][n] = List[i].vel[dim];
	  n++;
	}
      }
      
      for (j = 0; j < NumberOfParticleAttributes; j++) {
	n = NumberOfParticles;
	for (i = StartIndex; i < EndIndex; i++) {
	  Attribute[j][n] = List[i].attribute[j];
	  n++;
	}
      }
      
    } // ENDIF TotalNumberOfParticles > 0
 
    /* Set new number of particles in this grid. */
 
    NumberOfParticles = TotalNumberOfParticles;
 
    /* Copy new ones into place. */
 
    this->DeleteParticles();
    if (TotalNumberOfParticles > 0)
      this->SetParticlePointers(Mass, Number, Type, Position, Velocity,
				Attribute);
 
  } // end: if (COPY_IN)


 
  return SUCCESS;
}
