/***********************************************************************
/
/  COMMUNICATION ROUTINE: TRANSFER PARTICLES
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified:   May, 2009 by John Wise: optimized version to transfer
/                particles in one sweep with collective calls.
/
/  PURPOSE:
/
************************************************************************/
#ifdef OPTIMIZED_CTP
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
 
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
void my_exit(int status);
 
// function prototypes
 
#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_ParticleMoveList;
#endif
 
Eint32 compare_grid(const void *a, const void *b);
int Enzo_Dims_create(int nnodes, int ndims, int *dims);
int CommunicationAllSumIntegerValues(int *Values, int Number);
int CommunicationShareParticles(int *NumberToMove, particle_data* &SendList,
				int &NumberOfReceives,
				particle_data* &SharedList);

#define KEEP_PARTICLES_LOCAL
 
int CommunicationTransferParticles(grid *GridPointer[], int NumberOfGrids)
{

  if (NumberOfGrids == 1)
    return SUCCESS;
 
  int i, j, jstart, jend, dim, grid, proc;

  /* Assume that the grid was split by Enzo_Dims_create, and create a
     map from grid number to an index that is determined from (i,j,k)
     of the grid partitions. */

  int *GridMap = new int[NumberOfGrids];
  int Layout[MAX_DIMENSION], LayoutTemp[MAX_DIMENSION];
  int Rank, grid_num, GridPosition[MAX_DIMENSION], Dims[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];

  GridPointer[0]->ReturnGridInfo(&Rank, Dims, Left, Right); // need rank
  if (Enzo_Dims_create(NumberOfGrids, Rank, LayoutTemp) == FAIL) {
    fprintf(stderr, "Error in Enzo_Dims_create.\n");
    return FAIL;
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Layout[MAX_DIMENSION-1-dim] = LayoutTemp[dim];

  for (grid = 0; grid < NumberOfGrids; grid++) {
    GridPointer[grid]->ReturnGridInfo(&Rank, Dims, Left, Right);
    for (dim = 0; dim < Rank; dim++)
      GridPosition[dim] = 
	int(Layout[dim] * (Left[dim] - DomainLeftEdge[dim]) /
	    (DomainRightEdge[dim] - DomainLeftEdge[dim]));
    grid_num = GridPosition[0] + 
      Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
    GridMap[grid_num] = grid;
  } // ENDFOR grids
 
  int *NumberToMove = new int[NumberOfProcessors];
  particle_data *SendList = NULL;
  particle_data *SharedList = NULL;
 
  for (i = 0; i < NumberOfProcessors; i++)
    NumberToMove[i] = 0;

#ifdef DEBUG_CTP
  if (MyProcessorNumber == ROOT_PROCESSOR) 
  for (j = 0; j < NumberOfGrids; j++)
    fprintf(stderr, "Pa(%"ISYM") CTP grid[%"ISYM"] = %"ISYM"\n",
	    MyProcessorNumber, j, GridPointer[j]->ReturnNumberOfParticles());
#endif
 
  /* Generate the list of particle moves. */

  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("CTP: Before COPY_OUT for %"ISYM" grids.\n", NumberOfGrids);

  int Zero = 0;
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (GridPointer[grid]->
	CommunicationTransferParticles(GridPointer, NumberOfGrids, grid, 
				       NumberToMove, Zero, Zero, SendList, 
				       Layout, GridMap, COPY_OUT) == FAIL) {
      fprintf(stderr, "Error in grid->CommunicationTransferParticles(OUT).\n");
      return FAIL;
    }

  int TotalNumberToMove = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++)
    TotalNumberToMove += NumberToMove[proc];

  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("CTP: After COPY_OUT for %"ISYM" grids. Counted %"ISYM" particles.\n", 
	   NumberOfGrids, TotalNumberToMove);

  int NumberOfReceives = 0;
#ifdef KEEP_PARTICLES_LOCAL

  /* Sort by grid number.  We no longer share with other processors
     because the particles are stored locally until
     CommunicationCollectParticles(ALLGRIDS). */

  SharedList = SendList;
  NumberOfReceives = TotalNumberToMove;
  int particle_data_size = sizeof(particle_data);
  qsort(SharedList, TotalNumberToMove, particle_data_size, compare_grid);

#else

  /* Share the Particle moves. */

  if (CommunicationShareParticles(NumberToMove, SendList, NumberOfReceives,
				  SharedList) == FAIL) {
    fprintf(stderr, "Error in CommunicationShareParticles.\n");
    return FAIL;
  }

#endif

  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("CTP: sorted particle list.\n");

  /* Copy particles back to grids */

  jstart = 0;
  jend = 0;

  // Copy shared particles to grids, if any
  if (NumberOfReceives > 0) {
    for (j = 0; j < NumberOfGrids && jend < NumberOfReceives; j++) {
      while (SharedList[jend].grid <= j) {
	jend++;
	if (jend == NumberOfReceives) break;
      }
      if (GridPointer[j]->CommunicationTransferParticles
	  (GridPointer, NumberOfGrids, j, NumberToMove, jstart, jend, 
	   SharedList, Layout, GridMap, COPY_IN) == FAIL) {
	fprintf(stderr, "Error in grid->CommunicationTransferParticles(IN).\n");
	return FAIL;
      }
      jstart = jend;
    } // ENDFOR grids
  } // ENDIF NumberOfRecieves > 0

  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("CTP: after COPY_IN.\n");

  /* Even if we have no receives, we still have to remove the sent
     particles. */

#ifndef KEEP_PARTICLES_LOCAL
  else {

    for (j = 0; j < NumberOfGrids; j++)
      if (GridPointer[j]->ReturnProcessorNumber() == MyProcessorNumber)
	if (GridPointer[j]->CleanUpMovedParticles() == FAIL) {
	  fprintf(stderr, "Error in grid->CleanUpMovedParticles.\n");
	  return FAIL;
	}

  } // ENDELSE NumberOfReceives > 0
#else

  /* If the loop with grid::CTP stopped before we covered all of the
     grids, be sure to remove the particles on the remaining grids on
     all processors. If NumberOfReceives == 0, then no particles were
     moved. */

  if (NumberOfReceives > 0)
    for (j = SharedList[NumberOfReceives-1].grid; j < NumberOfGrids; j++)
      if (GridPointer[j]->CleanUpMovedParticles() == FAIL) {
	fprintf(stderr, "Error in grid->CleanUpMovedParticles.\n");
	return FAIL;
      }

  if (MyProcessorNumber == ROOT_PROCESSOR && NumberOfReceives > 0)
    printf("CTP: cleaned up particles in grids %"ISYM"->%"ISYM"\n",
	   SharedList[NumberOfReceives-1].grid, NumberOfGrids-1);
    
#endif /* KEEP_PARTICLES_LOCAL */
 
  /* Set number of particles so everybody agrees. */

  // (JHW, May 2009) No longer needed because we synchronize in
  // CommunicationCollectParticles now.
 
//  if (NumberOfProcessors > 1) {
//    int *AllNumberOfParticles = new int[NumberOfGrids];
//    for (j = 0; j < NumberOfGrids; j++)
//      if (GridPointer[j]->ReturnProcessorNumber() == MyProcessorNumber)
//	AllNumberOfParticles[j] = GridPointer[j]->ReturnNumberOfParticles();
//      else
//	AllNumberOfParticles[j] = 0;
//
//    CommunicationAllSumIntegerValues(AllNumberOfParticles, NumberOfGrids);
//    for (j = 0; j < NumberOfGrids; j++)
//      GridPointer[j]->SetNumberOfParticles(AllNumberOfParticles[j]);
//
//    delete [] AllNumberOfParticles;
//  }

  /* Cleanup. */

  if (SendList != SharedList)
    delete [] SendList;
  delete [] SharedList;
  delete [] NumberToMove;
  delete [] GridMap;

#ifdef USE_MPI
  int temp_int = TotalNumberToMove;
  MPI_Reduce(&temp_int, &TotalNumberToMove, 1, IntDataType, MPI_SUM, 
	     ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif
  if (debug)
    printf("CommunicationTransferParticles: moved = %"ISYM"\n",
  	   TotalNumberToMove);
 
  return SUCCESS;
}

#endif /* OPTIMIZED_CTP */
