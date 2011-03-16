/***********************************************************************
/
/  GRID CLASS (COPY STARS INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/  modified2:  John Wise
/  date:       July, 2009 (modified grid:CTP)
/
/  PURPOSE:
/
************************************************************************/
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
 
Star *PopStar(Star * &Node);
Star* StarBufferToList(StarBuffer *buffer, int n);
void InsertStarAfter(Star * &Node, Star * &NewNode);
void DeleteStarList(Star * &Node);
 
int grid::CommunicationTransferStars(grid* Grids[], int NumberOfGrids,
		 int ToGrid[6], int NumberToMove[6],
		 StarBuffer *StarData[6], int CopyDirection)
{
 
  /* Declarations. */
 
  int i, j, k, dim, grid;
  Star *cstar, *MoveStar, *StarList[6];
 
  /* ----------------------------------------------------------------- */
  /* Copy star out of grid. */
 
  if (CopyDirection == COPY_OUT) {

    /* If there are no stars to move, we're done. */
 
    if (NumberOfStars == 0)
      return SUCCESS;
 
    /* Count stars to move (mark already counted by setting mass < 0). */
 
    for (dim = 0; dim < GridRank; dim++)
      for (cstar = Stars; cstar; cstar = cstar->NextStar)
	if (cstar->Mass >= 0) {
	  if (cstar->pos[dim] < GridLeftEdge[dim]) {
	    NumberToMove[dim*2+0]++;
	    cstar->Mass *= -cstar->Mass;
	  }
	  if (cstar->pos[dim] > GridRightEdge[dim]) {
	    NumberToMove[dim*2+1]++;
	    cstar->Mass = -cstar->Mass;
	  }
	}
 
    /* Allocate space. */
 
    int TotalToMove = 0;
    for (i = 0; i < 6; i++) {
      TotalToMove += NumberToMove[i];
      StarList[i] = NULL;
      if (NumberToMove[i] > 0)
	StarData[i] = new StarBuffer[NumberToMove[i]];
      else
	StarData[i] = NULL;
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
 
    /* Move stars into linked lists for each face */

    for (dim = 0; dim < GridRank; dim++) {
      int n1 = 0, n2 = 0;
      NumberOfStars = 0;
      cstar = Stars;
      Stars = NULL;
      while (cstar != NULL) {

	MoveStar = PopStar(cstar);  // also advances to NextStar
 
	/* shift left. */
 
	if (MoveStar->pos[dim] < GridLeftEdge[dim]) {
	  MoveStar->Mass = -MoveStar->Mass;
	  InsertStarAfter(StarList[dim*2], MoveStar);
	}
 
	/* shift right. */
 
	else if (MoveStar->pos[dim] > GridRightEdge[dim]) {
	  MoveStar->Mass = -MoveStar->Mass;
	  InsertStarAfter(StarList[dim*2+1], MoveStar);
	}

	/* no shift */

	else {
	  InsertStarAfter(this->Stars, MoveStar);
	  NumberOfStars++;
	}
 
      } // ENDFOR stars
    } // end: loop over dims
 
    /* Convert linked lists into arrays and delete them */

    for (i = 0; i < 6; i++)
      if (NumberToMove[i] > 0) {
	StarData[i] = StarList[i]->StarListToBuffer(NumberToMove[i]);
	DeleteStarList(StarList[i]);
      }

  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy stars back into grid. */
 
  else {

    /* Count up total number. */
 
    int TotalNumberOfStars = NumberOfStars;
 
    for (i = 0; i < 6; i++)
      TotalNumberOfStars += NumberToMove[i];
 
    if (TotalNumberOfStars == 0 && NumberOfStars == 0)
      return SUCCESS;

    /* Copy stars from buffer into linked list */

    for (i = 0; i < 6; i++)
      if (NumberToMove[i] > 0) {
	MoveStar = StarBufferToList(StarData[i], NumberToMove[i]);
	MoveStar->CurrentGrid = this;
	MoveStar->GridID = this->ID;
	InsertStarAfter(this->Stars, MoveStar);
      } // ENDIF NumberToMove > 0
 
    /* Set new number of stars in this grid. */
 
    NumberOfStars = TotalNumberOfStars;
 
  } // end: if (COPY_IN)

 
  return SUCCESS;
}
