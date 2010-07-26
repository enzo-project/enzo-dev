/***********************************************************************
/
/  GRID CLASS (MOVE APPROPRIATE PARTICLES FROM SPECIFIED GRID TO THIS GRID)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
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
#include "Hierarchy.h"
 
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
 
 
int grid::MoveSubgridParticlesFast(int NumberOfSubgrids, grid* ToGrids[],
				   int AllLocal)
{

  if (debug1) 
    printf("MoveSubgridParticlesFast: %"ISYM"\n", NumberOfParticles);
 
  /* If there are no particles to move, we're done. */
 
  if (NumberOfParticles == 0 || NumberOfSubgrids == 0)
    return SUCCESS;
 
  int i, j, dim, index, subgrid, n;
 
  /* Initialize. */
 
  int *ParticlesToMove = new int[NumberOfSubgrids];
  for (i = 0; i < NumberOfSubgrids; i++)
    ParticlesToMove[i] = 0;
 
  /* Error check. */
 
  if (BaryonField[NumberOfBaryonFields] == NULL &&
      MyProcessorNumber == ProcessorNumber) {
    ENZO_FAIL("Subgrid field not present.\n");
  }
 
  /* Loop over particles and count the number in each subgrid. */
 
  int i0 = 0, j0 = 0, k0 = 0;
  if (MyProcessorNumber == ProcessorNumber) {
    for (i = 0; i < NumberOfParticles; i++) {
 
      /* Compute index of particle position. */
 
      i0 = int((ParticlePosition[0][i] - CellLeftEdge[0][0])/CellWidth[0][0]);
      if (GridRank > 0)
       j0 = int((ParticlePosition[1][i] - CellLeftEdge[1][0])/CellWidth[1][0]);
      if (GridRank > 1)
       k0 = int((ParticlePosition[2][i] - CellLeftEdge[2][0])/CellWidth[2][0]);
 
      i0 = max(min(GridEndIndex[0], i0), GridStartIndex[0]);
      j0 = max(min(GridEndIndex[1], j0), GridStartIndex[1]);
      k0 = max(min(GridEndIndex[2], k0), GridStartIndex[2]);
 
      index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;
 
      /* Find subgrid number of this particle, and add to count. */
 
      subgrid = nint(BaryonField[NumberOfBaryonFields][index])-1;
      if (subgrid >= 0)
	ParticlesToMove[subgrid]++;
      if (subgrid < -1 || subgrid > NumberOfSubgrids-1) {
	ENZO_VFAIL("particle subgrid (%"ISYM"/%"ISYM") out of range\n", subgrid,
		NumberOfSubgrids)
      }
 
    }  // end: loop over particles
  } // end: if (MyProcessorNumber)
 
  /* Communicate number of send particles to subgrids */
 
  if (AllLocal == FALSE)
    for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
      if (CommunicationBroadcastValue(&ParticlesToMove[subgrid],
				      ProcessorNumber) == FAIL) {
	ENZO_FAIL("Error in CommunicationBroadcastValue.\n");
      }
/*
    if ((MyProcessorNumber == ProcessorNumber ||
	 MyProcessorNumber == ToGrids[subgrid]->ProcessorNumber) &&
	ProcessorNumber != ToGrids[subgrid]->ProcessorNumber) {
      ENZO_FAIL("this routine not parallelized.\n");
      if (CommunicationSendInt(MyProcessorNumber,
			       ToGrids[subgrid]->ProcessorNumber,
			       &ParticlesToMove[subgrid]) == FAIL) {
        ENZO_FAIL("Error in CommunicationSendInt.\n");
      }
    }
*/
 
  /* Allocate space on all the subgrids with particles. */
 
  if (MyProcessorNumber == ProcessorNumber)
    for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
 
      if (ParticlesToMove[subgrid] > 0) {
 
	if (ToGrids[subgrid]->ParticlePosition[0] != NULL ||
	    ToGrids[subgrid]->NumberOfParticles != 0) {
	  ENZO_VFAIL("Particles already in subgrid %"ISYM" (n=%"ISYM", nm=%"ISYM")\n",
		  subgrid, ToGrids[subgrid]->NumberOfParticles,
		  ParticlesToMove[subgrid])
	}
 
	ToGrids[subgrid]->AllocateNewParticles(ParticlesToMove[subgrid]);
 
	if (debug1) printf("MoveSubgridParticles: subgrid[%"ISYM"] = %"ISYM"\n",
			  subgrid, ParticlesToMove[subgrid]);
 
      } // end: if (ParticlesToMove > 0)
 
  /* Compute the increase in mass for particles moving to the subgrid. */
 
  float RefinementFactors[MAX_DIMENSION], MassIncrease = 1.0;
  this->ComputeRefinementFactorsFloat(ToGrids[0], RefinementFactors);
  for (dim = 0; dim < GridRank; dim++)
    MassIncrease *= RefinementFactors[dim];
 
  if (MyProcessorNumber == ProcessorNumber) {
 
    /* Loop over particles and move them to the appropriate ToGrid, depending
       on their position. */
 
    for (i = 0; i < NumberOfParticles; i++) {
 
      /* Compute index of particle position. */
 
      i0 = int((ParticlePosition[0][i] - CellLeftEdge[0][0])/CellWidth[0][0]);
      if (GridRank > 0)
       j0 = int((ParticlePosition[1][i] - CellLeftEdge[1][0])/CellWidth[1][0]);
      if (GridRank > 1)
       k0 = int((ParticlePosition[2][i] - CellLeftEdge[2][0])/CellWidth[2][0]);
 
      i0 = max(min(GridEndIndex[0], i0), GridStartIndex[0]);
      j0 = max(min(GridEndIndex[1], j0), GridStartIndex[1]);
      k0 = max(min(GridEndIndex[2], k0), GridStartIndex[2]);
 
      index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;
 
      /* Find subgrid number of this particle, and move it. */
 
      subgrid = nint(BaryonField[NumberOfBaryonFields][index])-1;
 
      if (subgrid >= 0) {
	n = ToGrids[subgrid]->NumberOfParticles;
	ToGrids[subgrid]->ParticleMass[n] = ParticleMass[i] * MassIncrease;
	ToGrids[subgrid]->ParticleNumber[n] = ParticleNumber[i];
	ToGrids[subgrid]->ParticleType[n] = ParticleType[i];
	for (dim = 0; dim < GridRank; dim++) {
	 ToGrids[subgrid]->ParticlePosition[dim][n] = ParticlePosition[dim][i];
	 ToGrids[subgrid]->ParticleVelocity[dim][n] = ParticleVelocity[dim][i];
        }
	for (j = 0; j < NumberOfParticleAttributes; j++)
	  ToGrids[subgrid]->ParticleAttribute[j][n] = ParticleAttribute[j][i];
 
	ToGrids[subgrid]->NumberOfParticles++;
 
	/* Mark this particle removed. */
 
	ParticleMass[i] = FLOAT_UNDEFINED;
 
      } // end: if (subgrid >= 0)
 
    } // end: loop over particles
 
    /* Clean up the moved particles. */
 
    this->CleanUpMovedParticles();
 
    delete [] BaryonField[NumberOfBaryonFields];
    BaryonField[NumberOfBaryonFields] = NULL;
 
  } // end: if (MyProcessorNumber)
 
  /* Transfer particles from fake to real grids (and clean up). */
 
  for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
    if ((MyProcessorNumber == ProcessorNumber ||
         MyProcessorNumber == ToGrids[subgrid]->ProcessorNumber) &&
	ProcessorNumber != ToGrids[subgrid]->ProcessorNumber)
      if (ParticlesToMove[subgrid] != 0) {
	if (this->CommunicationSendParticles(ToGrids[subgrid],
             ToGrids[subgrid]->ProcessorNumber, 0, ParticlesToMove[subgrid], 0)
	    == FAIL) {
	  ENZO_FAIL("Error in grid->CommunicationSendParticles.\n");
	}
	if (MyProcessorNumber == ProcessorNumber)

	  ToGrids[subgrid]->DeleteAllFields();
      }
 
  delete [] ParticlesToMove;
 
  return SUCCESS;
}
