/***********************************************************************
/
/  GRID CLASS (COPY PARTICLES INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/  modified2:  May, 2009 by John Wise: optimized version to transfer
/                particles in one sweep with collective calls.
/
/  PURPOSE:
/
************************************************************************/
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <string.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"

int search_lower_bound(int *arr, int value, int low, int high, 
		       int total);
 
int grid::CommunicationTransferParticles(grid* Grids[], int NumberOfGrids,
					 int ThisGridNum, int TopGridDims[],
					 int *&NumberToMove, 
					 int StartIndex, int EndIndex, 
					 particle_data *&List,
					 int *Layout, int *GStartIndex[], 
					 int *GridMap, 
					 int CopyDirection)
{
 
  /* Declarations. */
 
  int i, j, k, dim, grid, proc, grid_num, bin, CenterIndex;
  int GridPosition[MAX_DIMENSION];
  int *ToGrid, *pbin;
  FLOAT r[MAX_DIMENSION];
 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    GridPosition[dim] = 0;

  /* ----------------------------------------------------------------- */
  /* Copy particle out of grid. */
 
  if (CopyDirection == COPY_OUT) {

    /* If there are no particles to move, we're done. */
 
    if (NumberOfParticles == 0)
      return SUCCESS;

//    if (MyProcessorNumber != ProcessorNumber)
//      return SUCCESS;
 
    /* Count the number of particles already moved */

    int PreviousTotalToMove = 0;
    for (i = 0; i < NumberOfProcessors; i++)
      PreviousTotalToMove += NumberToMove[i];

    /* Count particles to move.  Apply perioidic wrap to the
       particles. */
 
    ToGrid = new int[NumberOfParticles];

    float DomainWidth[MAX_DIMENSION], DomainWidthInv[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
      DomainWidthInv[dim] = 1.0/DomainWidth[dim];
    }

    // Periodic boundaries
    for (dim = 0; dim < GridRank; dim++) 
      for (i = 0; i < NumberOfParticles; i++) {
	if (ParticlePosition[dim][i] > DomainRightEdge[dim])
	  ParticlePosition[dim][i] -= DomainWidth[dim];
	else if (ParticlePosition[dim][i] < DomainLeftEdge[dim])
	  ParticlePosition[dim][i] += DomainWidth[dim];
      }

    for (i = 0; i < NumberOfParticles; i++) {

      for (dim = 0; dim < GridRank; dim++) {

	if (Layout[dim] == 1) {
	  GridPosition[dim] = 0;
	} else {

	  CenterIndex = 
	  (int) (TopGridDims[dim] * 
		 (ParticlePosition[dim][i] - DomainLeftEdge[dim]) *
		 DomainWidthInv[dim]);

	  GridPosition[dim] = 
	    search_lower_bound(GStartIndex[dim], CenterIndex, 0, Layout[dim],
			       Layout[dim]);
	  GridPosition[dim] = min(GridPosition[dim], Layout[dim]-1);
	  
	} // ENDELSE

      } // ENDFOR dim

      grid_num = GridPosition[0] + 
	Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);

      grid = GridMap[grid_num];
      if (grid != ThisGridNum) {
	proc = Grids[grid]->ReturnProcessorNumber();
	NumberToMove[proc]++;
//	printf("grid %d->%d: Pos = %f %f %f\n", ThisGridNum, grid, 
//	       ParticlePosition[0][i], ParticlePosition[1][i], 
//	       ParticlePosition[2][i]);
      }

      ToGrid[i] = grid;

    } // ENDFOR particles

    /* Allocate space. */
 
    int TotalToMove = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++)
      TotalToMove += NumberToMove[proc];

    if (TotalToMove > PreviousTotalToMove) {
 
      // Increase the size of the list to include the particles from
      // this grid

      particle_data *NewList = new particle_data[TotalToMove];
      memcpy(NewList, List, PreviousTotalToMove * sizeof(particle_data));
      delete [] List;
      List = NewList;
      //      particle_data *TempList = List;
      //      List = new particle_data[TotalToMove];
      //      for (i = 0; i < PreviousTotalToMove; i++)
      //	List[i] = TempList[i];
      //      delete [] TempList;
 
      /* Move particles (mark those moved by setting mass = FLOAT_UNDEFINED). */

      int n1 = PreviousTotalToMove;

      for (i = 0; i < NumberOfParticles; i++) {
	grid = ToGrid[i];
	if (grid != ThisGridNum) {
	  proc = Grids[grid]->ReturnProcessorNumber();
	  for (dim = 0; dim < GridRank; dim++) {
	    List[n1].pos[dim] = ParticlePosition[dim][i];
	    List[n1].vel[dim] = ParticleVelocity[dim][i];
	  }
	  List[n1].mass = ParticleMass[i];
	  List[n1].id = ParticleNumber[i];
	  List[n1].type = ParticleType[i];
	  for (j = 0; j < NumberOfParticleAttributes; j++)
	    List[n1].attribute[j] = ParticleAttribute[j][i];
	  List[n1].grid = grid;
	  List[n1].proc = MyProcessorNumber;
	  ParticleMass[i] = FLOAT_UNDEFINED;
	  n1++;
	} // ENDIF different grid
      } // ENDFOR particles

    } // ENDIF TotalToMove > PreviousTotalToMove

    delete [] ToGrid;
 
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
 
    int TotalNumberOfParticles;
    int NumberOfNewParticles = EndIndex - StartIndex;

    TotalNumberOfParticles = NumberOfParticles + NumberOfNewParticles;
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleMass[i] == FLOAT_UNDEFINED)
	TotalNumberOfParticles--;
 
    /* Allocate space for the particles. */

#ifdef USE_MPI
    T1 = MPI_Wtime();
#endif
 
    FLOAT *Position[MAX_DIMENSION];
    float *Velocity[MAX_DIMENSION], *Mass,
      *Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    PINT *Number;
    int *Type;

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

#ifdef USE_MPI
      T2 = MPI_Wtime();
#endif
 
      /* Copy this grid's particles to the new space. */

      int n = 0;

      for (i = 0; i < NumberOfParticles; i++)
	if (ParticleMass[i] >= 0) {
	  Mass[n] = ParticleMass[i];
	  Number[n] = ParticleNumber[i];
	  Type[n] = ParticleType[i];
	  n++;
	}

      for (dim = 0; dim < GridRank; dim++) {
	n = 0;
	for (i = 0; i < NumberOfParticles; i++)
	  if (ParticleMass[i] >= 0) {
	    Position[dim][n] = ParticlePosition[dim][i];
	    Velocity[dim][n] = ParticleVelocity[dim][i];
	    n++;
	  }
      }
	
      for (j = 0; j < NumberOfParticleAttributes; j++) {
	n = 0;
	for (i = 0; i < NumberOfParticles; i++)
	  if (ParticleMass[i] >= 0) {
	    Attribute[j][n] = ParticleAttribute[j][i];
	    n++;
	  }
      }

#ifdef USE_MPI
      T3 = MPI_Wtime();
#endif
 
      /* Copy new particles. Periodic wrap is now done in the COPY_OUT. */

      for (i = StartIndex; i < EndIndex; i++) {

	for (dim = 0; dim < GridRank; dim++) {
	  Position[dim][n] = List[i].pos[dim];
	  Velocity[dim][n] = List[i].vel[dim];
	} // ENDFOR dim
      
	Mass[n] = List[i].mass;
	Number[n] = List[i].id;
	Type[n] = List[i].type;
	for (j = 0; j < NumberOfParticleAttributes; j++)
	  Attribute[j][n] = List[i].attribute[j];
      
	n++;

      } // ENDFOR particles

#ifdef USE_MPI
      T4 = MPI_Wtime();
#endif

    } // ENDIF TotalNumberOfParticles > 0
 
    /* Set new number of particles in this grid. */
 
    NumberOfParticles = TotalNumberOfParticles;
 
    /* Copy new ones into place. */
 
    this->DeleteParticles();
    if (NumberOfParticles > 0)
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
