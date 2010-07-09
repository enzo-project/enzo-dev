/***********************************************************************
/
/  GRID CLASS (COPY PARTICLES INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/
/  PURPOSE:
/
************************************************************************/
#ifndef OPTIMIZED_CTP 
//
 
#ifdef USE_MPI
#include <mpi.h>
#endif
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
 
 
int grid::CommunicationTransferParticles(grid* Grids[], int NumberOfGrids,
		 int ToGrid[6], int NumberToMove[6],
		 float_int *ParticleData[6], int CopyDirection)
{
 
  /* Declarations. */
 
  int i, j, k, dim, grid;
 
  /* ----------------------------------------------------------------- */
  /* Copy particle out of grid. */
 
  if (CopyDirection == COPY_OUT) {

    /* If there are no particles to move, we're done. */
 
    if (NumberOfParticles == 0)
      return SUCCESS;
 
    /* Count particles to move (mark already counted by setting mass < 0). */
 
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < NumberOfParticles; i++)
	if (ParticleMass[i] >= 0) {
	  if (ParticlePosition[dim][i] < GridLeftEdge[dim]) {
	    NumberToMove[dim*2+0]++;
	    ParticleMass[i] = -ParticleMass[i];
	  }
	  if (ParticlePosition[dim][i] > GridRightEdge[dim]) {
	    NumberToMove[dim*2+1]++;
	    ParticleMass[i] = -ParticleMass[i];
	  }
	}
 
    /* Allocate space. */
 
    int TotalToMove = 0, NumberOfParticleFields = 9+NumberOfParticleAttributes;
    for (i = 0; i < 6; i++) {
      TotalToMove += NumberToMove[i];
      if (NumberToMove[i] > 0) {
//	fprintf(stdout, "P%"ISYM", grid %x: side %"ISYM", NumberToMove = %"ISYM"\n",
//	       ProcessorNumber, this, i, NumberToMove[i]);
	fflush(stdout);
	ParticleData[i] = new float_int[NumberToMove[i]*
				        NumberOfParticleFields];
      } else
	ParticleData[i] = NULL;
    }
 
    /* Set ToGrid. */
 
    for (dim = 0; dim < GridRank; dim++) {
      int DimSize = nint((DomainRightEdge[dim] -
		      DomainLeftEdge[dim])/CellWidth[dim][0]);
 
      /* Find Left grid */
 
      for (grid = 0; grid < NumberOfGrids; grid++) {
	int FoundIt = TRUE;
	for (i = 0; i < GridRank; i++) {
	  if (dim != i && nint(GridLeftEdge[i]/CellWidth[i][0]) !=
	      nint(Grids[grid]->GridLeftEdge[i]/CellWidth[i][0]))
	    FoundIt = FALSE;
	  if (dim == i && (nint(
               Grids[grid]->GridRightEdge[i]/CellWidth[i][0]) % DimSize)
	      != nint(GridLeftEdge[i]/CellWidth[i][0]))
	    FoundIt = FALSE;
	}
	if (FoundIt) {
	  ToGrid[dim*2+0] = grid;
	  break;
	}
      }
 
      /* Find right grid */
 
      for (grid = 0; grid < NumberOfGrids; grid++) {
	int FoundIt = TRUE;
	for (i = 0; i < GridRank; i++) {
	  if (dim != i && nint(GridLeftEdge[i]/CellWidth[i][0]) !=
	      nint(Grids[grid]->GridLeftEdge[i]/CellWidth[i][0]))
	    FoundIt = FALSE;
	  if (dim == i && (nint(
               GridRightEdge[i]/CellWidth[i][0]) % DimSize)
	      != nint(Grids[grid]->GridLeftEdge[i]/CellWidth[i][0]))
	    FoundIt = FALSE;
	}
	if (FoundIt) {
	  ToGrid[dim*2+1] = grid;
	  break;
	}
      }
 
    } // end loop over dims
 
    if (TotalToMove == 0)
      return SUCCESS;
 
    /* Move particles (mark those moved by setting mass = FLOAT_UNDEFINED). */
 
    for (dim = 0; dim < GridRank; dim++) {
      int n1 = 0, n2 = 0;
      for (i = 0; i < NumberOfParticles; i++)
	if (ParticleMass[i] != FLOAT_UNDEFINED) {
 
	  /* shift left. */
 
	  if (ParticlePosition[dim][i] < GridLeftEdge[dim]) {
	    for (j = 0; j < GridRank; j++) {
	      ParticleData[dim*2][n1++].FVAL = ParticlePosition[j][i];
	      ParticleData[dim*2][n1++].fval = ParticleVelocity[j][i];
	    }
	    ParticleData[dim*2][n1++].fval = -ParticleMass[i];
	    ParticleData[dim*2][n1++].IVAL = ParticleNumber[i];
	    ParticleData[dim*2][n1++].ival = ParticleType[i];
	    for (j = 0; j < NumberOfParticleAttributes; j++)
	      ParticleData[dim*2][n1++].fval = ParticleAttribute[j][i];
	    ParticleMass[i] = FLOAT_UNDEFINED;
	  }
 
	  /* shift right. */
 
	  if (ParticlePosition[dim][i] > GridRightEdge[dim]) {
	    for (j = 0; j < MAX_DIMENSION; j++) {
	      ParticleData[dim*2+1][n2++].FVAL = ParticlePosition[j][i];
	      ParticleData[dim*2+1][n2++].fval = ParticleVelocity[j][i];
	    }
	    ParticleData[dim*2+1][n2++].fval = -ParticleMass[i];
	    ParticleData[dim*2+1][n2++].IVAL = ParticleNumber[i];
	    ParticleData[dim*2+1][n2++].ival = ParticleType[i];
	    for (j = 0; j < NumberOfParticleAttributes; j++)
	      ParticleData[dim*2+1][n2++].fval = ParticleAttribute[j][i];
	    ParticleMass[i] = FLOAT_UNDEFINED;
	  }
 
	} // end: if (ParticleMass[i] != FLOAT_UNDEFINED)
    } // end: loop over dims
 
  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy particle back into grid. */
 
  else {

    double t00, t01, ttt;
    double T1, T2, T3, T4;

    t00 = 0.0;
    t01 = 0.0;
    ttt = 0.0;

    T1 = 0.0;
    T2 = 0.0;
    T3 = 0.0;
    T4 = 0.0;

#ifdef USE_MPI
    t00 = MPI_Wtime();
#endif
 
    /* Count up total number. */
 
    int TotalNumberOfParticles = NumberOfParticles;
 
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleMass[i] == FLOAT_UNDEFINED)
	TotalNumberOfParticles--;
 
    for (i = 0; i < 6; i++)
      TotalNumberOfParticles += NumberToMove[i];
 
    if (TotalNumberOfParticles == 0 && NumberOfParticles == 0)
      return SUCCESS;
 
    /* Allocate space for the particles. */

#ifdef USE_MPI
    T1 = MPI_Wtime();
#endif
 
    FLOAT *Position[MAX_DIMENSION];
    float *Velocity[MAX_DIMENSION], *Mass,
          *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    PINT *Number;
    int *Type;
 
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

#ifdef USE_MPI
    T2 = MPI_Wtime();
#endif
 
    /* Copy this grid's particles to the new space (don't copy any erased
       particles, those with Mass == FLOAT_UNDEFINED). */

/* 
    int n = 0;
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleMass[i] != FLOAT_UNDEFINED) {
	Mass[n]   = ParticleMass[i];
	Number[n] = ParticleNumber[i];
	Type[n] = ParticleType[i];
 
	for (dim = 0; dim < GridRank; dim++) {
	  Position[dim][n] = ParticlePosition[dim][i];
	  Velocity[dim][n] = ParticleVelocity[dim][i];
	}
	for (j = 0; j < NumberOfParticleAttributes; j++)
	  Attribute[j][n] = ParticleAttribute[j][i];
	n++;
      }
*/

    int n;

    n = 0;
    for (i = 0; i < NumberOfParticles; i++) {
      if (ParticleMass[i] != FLOAT_UNDEFINED) {
        Mass[n]   = ParticleMass[i];
        Number[n] = ParticleNumber[i];
        Type[n] = ParticleType[i];
        n++;
      }
    }

    for (dim = 0; dim < GridRank; dim++) {
      n = 0;
      for (i = 0; i < NumberOfParticles; i++) {
        if (ParticleMass[i] != FLOAT_UNDEFINED) {
          Position[dim][n] = ParticlePosition[dim][i];
          Velocity[dim][n] = ParticleVelocity[dim][i];
          n++;
        }
      }
    }

    for (j = 0; j < NumberOfParticleAttributes; j++) {
      n = 0;
      for (i = 0; i < NumberOfParticles; i++) {
        if (ParticleMass[i] != FLOAT_UNDEFINED) {
          Attribute[j][n] = ParticleAttribute[j][i];
          n++;
        }
      }
    }

#ifdef USE_MPI
    T3 = MPI_Wtime();
#endif
 
    /* Copy new particles (and wrap around edges). */
 
    for (j = 0; j < 6; j++)
      if (NumberToMove[j] > 0) {
	int n2 = 0;
	for (i = 0; i < NumberToMove[j]; i++) {
 
	  for (dim = 0; dim < GridRank; dim++) {
	    Position[dim][n] = ParticleData[j][n2++].FVAL;
	    if (Position[dim][n] > DomainRightEdge[dim])
	      Position[dim][n] -= DomainRightEdge[dim] - DomainLeftEdge[dim];
	    if (Position[dim][n] < DomainLeftEdge[dim])
	      Position[dim][n] += DomainRightEdge[dim] - DomainLeftEdge[dim];
	    Velocity[dim][n] = ParticleData[j][n2++].fval;
	  }
 
	  Mass[n]   = ParticleData[j][n2++].fval;
	  Number[n] = ParticleData[j][n2++].IVAL;
	  Type[n] = ParticleData[j][n2++].ival;
	  for (k = 0; k < NumberOfParticleAttributes; k++)
	    Attribute[k][n] = ParticleData[j][n2++].fval;
 
	  n++;
	}
      }

#ifdef USE_MPI
    T4 = MPI_Wtime();
#endif
 
    /* Set new number of particles in this grid. */
 
    NumberOfParticles = TotalNumberOfParticles;
 
    /* Delete old pointers and copy new ones into place. */
 
    this->DeleteParticles();
    this->SetParticlePointers(Mass, Number, Type, Position, Velocity,
			      Attribute);
 
#ifdef USE_MPI
    t01 = MPI_Wtime();
#endif
    ttt = t01-t00;
    // fprintf(stderr, "COPY IN %"ISYM" : %16.6e : %16.6e %16.6e %16.6e %16.6e %16.6e\n", MyProcessorNumber, ttt,
    //         T1-t00, T2-T1, T3-T2, T4-T3, t01-T4);

  } // end: if (COPY_IN)


 
  return SUCCESS;
}
#endif /* OPTIMIZED_CTP */
