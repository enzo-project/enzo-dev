/***********************************************************************
/
/  GRID CLASS (MOVE APPROPRIATE STARS FROM SPECIFIED GRID TO THIS GRID)
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  NOTES: Adapted from grid::MoveSubgridParticlesFast()
/
************************************************************************/

//
#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
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
#include "CommunicationUtilities.h"

Star *PopStar(Star * &Node);
void InsertStarAfter(Star * &Node, Star * &NewNode);


int grid::MoveSubgridStars(int NumberOfSubgrids, grid* ToGrids[],
			   int AllLocal)
{

  /* If there are no stars to move, we're done. */

  if (NumberOfStars == 0 || NumberOfSubgrids == 0)
    return SUCCESS;

  Star *cstar, *MoveStar;
  int i, j, dim, index, subgrid, n;

  /* Initialize. */

  int *StarsToMove = new int[NumberOfSubgrids];
  for (i = 0; i < NumberOfSubgrids; i++)
    StarsToMove[i] = 0;

  /* Error check. */

  if (BaryonField[NumberOfBaryonFields] == NULL &&
      MyProcessorNumber == ProcessorNumber) {
    ENZO_FAIL("Subgrid field not present.\n");
  }

  /* Loop over particles and count the number in each subgrid. */

  int i0 = 0, j0 = 0, k0 = 0;
  if (MyProcessorNumber == ProcessorNumber) {
    cstar = Stars;
    Stars = NULL;
    while (cstar != NULL) {

      /* Compute index of particle position. */

      i0 = int((cstar->pos[0] - CellLeftEdge[0][0])/CellWidth[0][0]);
      if (GridRank > 1)
	j0 = int((cstar->pos[1] - CellLeftEdge[1][0])/CellWidth[1][0]);
      if (GridRank > 2)
	k0 = int((cstar->pos[2] - CellLeftEdge[2][0])/CellWidth[2][0]);

      i0 = max(min(GridEndIndex[0], i0), GridStartIndex[0]);
      j0 = max(min(GridEndIndex[1], j0), GridStartIndex[1]);
      k0 = max(min(GridEndIndex[2], k0), GridStartIndex[2]);

      index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;

      /* Find subgrid number of this particle, add to count, and move
	 star to "fake" grid. */

      subgrid = nint(BaryonField[NumberOfBaryonFields][index])-1;

      if (subgrid < -1 || subgrid > NumberOfSubgrids-1) {
	ENZO_VFAIL("particle subgrid (%"ISYM"/%"ISYM") out of range\n", subgrid,
		NumberOfSubgrids)
      }

      MoveStar = PopStar(cstar);  // also advances to NextStar

      if (subgrid >= 0) {

	MoveStar->CurrentGrid = ToGrids[subgrid];
	MoveStar->IncreaseLevel();
	MoveStar->GridID = ToGrids[subgrid]->ID;
	InsertStarAfter(ToGrids[subgrid]->Stars, MoveStar);

	this->NumberOfStars--;
	ToGrids[subgrid]->NumberOfStars++;
	StarsToMove[subgrid]++;

      } // end: if (subgrid >= 0)
      else
	InsertStarAfter(this->Stars, MoveStar);

    }  // end: loop over stars

  } // end: if (MyProcessorNumber)

  /* Communicate number of send stars to subgrids */

  if (AllLocal == FALSE)
#ifdef USE_MPI
    CommunicationAllReduceValues(StarsToMove, NumberOfSubgrids, MPI_MAX);
#endif

  /* Transfer stars to other processors (and clean up). */

  for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
    if ((MyProcessorNumber == ProcessorNumber ||
         MyProcessorNumber == ToGrids[subgrid]->ProcessorNumber) &&
	ProcessorNumber != ToGrids[subgrid]->ProcessorNumber)
      if (StarsToMove[subgrid] > 0) {
	if (this->CommunicationSendStars(ToGrids[subgrid], 
		  ToGrids[subgrid]->ProcessorNumber) == FAIL) {
	  ENZO_FAIL("Error in grid->CommunicationSendStars.\n");

	}
      }

  delete [] StarsToMove;

  return SUCCESS;
}
