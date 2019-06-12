/***********************************************************************
/
/  GRID CLASS (MOVE APPROPRIATE ACTIVE PARTICLES FROM SPECIFIED GRID 
/              TO THIS GRID)
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  John Wise -- re-purposed for active particles
/  date:       December, 2011
/
/  NOTES: Adapted from grid::MoveSubgridParticlesFast()
/
************************************************************************/

//
#ifdef USE_MPI
#endif

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CommunicationUtilities.h"
#include "ActiveParticle.h"

int grid::MoveSubgridActiveParticles(int NumberOfSubgrids, grid* ToGrids[],
				     int AllLocal)
{

  /* If there are no stars to move, we're done. */

  if (NumberOfActiveParticles == 0 || NumberOfSubgrids == 0)
    return SUCCESS;

  int i, j, dim, index, index1, index2, n;

  /* Initialize. */

  int *NumberToMove = new int[NumberOfSubgrids];
  for (i = 0; i < NumberOfSubgrids; i++)
    NumberToMove[i] = 0;

  /* Error check. */

  if (BaryonField[NumberOfBaryonFields] == NULL &&
      MyProcessorNumber == ProcessorNumber) {
    ENZO_FAIL("Subgrid field not present.\n");
  }

  int *subgrid = new int[NumberOfActiveParticles];

  /* Loop over particles and count the number in each subgrid. */

  int i0 = 0, j0 = 0, k0 = 0;
  int NumberToMoveLocal = 0;
  ActiveParticleType *np;
  ActiveParticleList<ActiveParticleType> *OldParticles = NULL;

  if (MyProcessorNumber == ProcessorNumber) {
    for (i = 0; i < NumberOfActiveParticles; i++) {

      /* Compute index of particle position. */

      i0 = int((this->ActiveParticles[i]->pos[0] - CellLeftEdge[0][0]) / 
          CellWidth[0][0]);
      if (GridRank > 1)
        j0 = int((this->ActiveParticles[i]->pos[1] - CellLeftEdge[1][0]) / 
            CellWidth[1][0]);
      if (GridRank > 2)
        k0 = int((this->ActiveParticles[i]->pos[2] - CellLeftEdge[2][0]) / 
            CellWidth[2][0]);

      i0 = max(min(GridEndIndex[0], i0), GridStartIndex[0]);
      j0 = max(min(GridEndIndex[1], j0), GridStartIndex[1]);
      k0 = max(min(GridEndIndex[2], k0), GridStartIndex[2]);

      index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;

      /* Find subgrid number of this particle, add to count, and move
	 star to "fake" grid. */

      subgrid[i] = nint(BaryonField[NumberOfBaryonFields][index])-1;
      if (subgrid[i] >= 0) NumberToMoveLocal++;

      if (subgrid[i] < -1 || subgrid[i] > NumberOfSubgrids-1) {
	ENZO_VFAIL("particle subgrid (%"ISYM"/%"ISYM") out of range\n", subgrid[i],
		NumberOfSubgrids)
      }

    } // ENDFOR particles

    float RefinementFactors[MAX_DIMENSION], MassIncrease = 1.0;
    this->ComputeRefinementFactorsFloat(ToGrids[0], RefinementFactors);
    for (dim = 0; dim < GridRank; dim++)
      MassIncrease *= RefinementFactors[dim];

    index = 0;
    for (i = 0; i < NumberOfActiveParticles; i++) {
      if (subgrid[i] >= 0) {
        np = this->ActiveParticles[i]->clone();
        this->ActiveParticles.erase(i);

        np->CurrentGrid = ToGrids[subgrid[i]];
        np->IncreaseLevel();
        np->AdjustMassByFactor(MassIncrease);
        np->GridID = ToGrids[subgrid[i]]->ID;
        
        ToGrids[subgrid[i]]->AddActiveParticle(np);
        delete np;

        this->NumberOfActiveParticles--;
        ToGrids[subgrid[i]]->NumberOfActiveParticles++;
        NumberToMove[subgrid[i]]++;
        
      } // end: if (subgrid[i] >= 0)

    }  // end: loop over stars
    
  } // end: if (MyProcessorNumber)

  /* Communicate number of send stars to subgrids */

#ifdef USE_MPI
  if (AllLocal == FALSE)
    CommunicationAllReduceValues(NumberToMove, NumberOfSubgrids, MPI_MAX);
#endif

  /* Transfer stars to other processors (and clean up). */

  int sg;
  for (sg = 0; sg < NumberOfSubgrids; sg++)
    if ((MyProcessorNumber == ProcessorNumber ||
            MyProcessorNumber == ToGrids[sg]->ProcessorNumber) &&
        ProcessorNumber != ToGrids[sg]->ProcessorNumber)
      if (NumberToMove[sg] > 0) {
        if (this->CommunicationSendActiveParticles(ToGrids[sg], 
                ToGrids[sg]->ProcessorNumber) == FAIL) {
          ENZO_FAIL("Error in grid->CommunicationSendStars.\n");
          
        }
      }

  delete[] NumberToMove;
  delete[] subgrid;

  return SUCCESS;
}
