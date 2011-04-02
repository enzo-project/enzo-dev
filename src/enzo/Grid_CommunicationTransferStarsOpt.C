/***********************************************************************
/
/  GRID CLASS (COPY STARS INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/  modified2:  May, 2009 by John Wise: optimized version to transfer
/                particles in one sweep with collective calls.
/  modified3:  July, 2009 by John Wise: adapted for stars
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

Star *PopStar(Star * &Node);
Star* StarBufferToList(StarBuffer buffer);
void InsertStarAfter(Star * &Node, Star * &NewNode);
int search_lower_bound(int *arr, int value, int low, int high, 
		       int total);
 
int grid::CommunicationTransferStars(grid* Grids[], int NumberOfGrids,
				     int ThisGridNum, int TopGridDims[],
				     int *&NumberToMove, 
				     int StartIndex, int EndIndex, 
				     star_data *&List, int *Layout, 
				     int *GStartIndex[], int *GridMap, 
				     int CopyDirection)
{
 
  /* Declarations. */
 
  int i, j, k, dim, grid, proc, grid_num, width, bin, CenterIndex;
  int GridPosition[MAX_DIMENSION];
  FLOAT r[MAX_DIMENSION];
  int *ToGrid, *pbin;
  Star *cstar, *MoveStar;

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    GridPosition[dim] = 0;
 
  /* ----------------------------------------------------------------- */
  /* Copy stars out of grid. */
 
  if (CopyDirection == COPY_OUT) {

    /* If there are no stars to move, we're done. */
 
    if (NumberOfStars == 0)
      return SUCCESS;

//    if (MyProcessorNumber != ProcessorNumber)
//      return SUCCESS;
 
    /* Count the number of stars already moved */

    int PreviousTotalToMove = 0;
    for (i = 0; i < NumberOfProcessors; i++)
      PreviousTotalToMove += NumberToMove[i];

    /* Count stars to move.  Apply perioidic wrap to the stars. */
 
    ToGrid = new int[NumberOfStars];

    float DomainWidth[MAX_DIMENSION], DomainWidthInv[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
      DomainWidthInv[dim] = 1.0/DomainWidth[dim];
    }

    // Periodic boundaries
    for (dim = 0; dim < GridRank; dim++) 
      for (cstar = Stars; cstar; cstar = cstar->NextStar) {
	if (cstar->pos[dim] > DomainRightEdge[dim])
	  cstar->pos[dim] -= DomainWidth[dim];
	else if (cstar->pos[dim] < DomainLeftEdge[dim])
	  cstar->pos[dim] += DomainWidth[dim];
      }

    for (cstar = Stars, i = 0; cstar; cstar = cstar->NextStar, i++) {

      for (dim = 0; dim < GridRank; dim++) {

	if (Layout[dim] == 1) {
	  GridPosition[dim] = 0;
	} else {

	  CenterIndex = 
	  (int) (TopGridDims[dim] * 
		 (cstar->pos[dim] - DomainLeftEdge[dim]) *
		 DomainWidthInv[dim]);

	  GridPosition[dim] = 
	    search_lower_bound(GStartIndex[dim], CenterIndex, 0, Layout[dim],
			       Layout[dim]);
	  GridPosition[dim] = min(GridPosition[dim], Layout[dim]-1);

	} // ENDELSE Layout

      } // ENDFOR dim

      grid_num = GridPosition[0] + 
	Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
      grid = GridMap[grid_num];
      if (grid != ThisGridNum) {
	proc = Grids[grid]->ReturnProcessorNumber();
	NumberToMove[proc]++;
      }

      ToGrid[i] = grid;

    } // ENDFOR stars

    /* Allocate space. */
 
    int TotalToMove = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++)
      TotalToMove += NumberToMove[proc];

    if (TotalToMove > PreviousTotalToMove) {
 
      // Increase the size of the list to include the stars from this grid

      star_data *NewList = new star_data[TotalToMove];
      memcpy(NewList, List, PreviousTotalToMove * sizeof(star_data));
      delete [] List;
      List = NewList;
 
      /* Move stars into list */

      int n1 = PreviousTotalToMove;
      StarBuffer *TempBuffer;

      cstar = Stars;
      Stars = NULL;
      NumberOfStars = 0;

      i = 0;
      while (cstar != NULL) {

	MoveStar = PopStar(cstar);  // also advances to NextStar
	grid = ToGrid[i];

	if (grid != ThisGridNum) {
	  TempBuffer = MoveStar->StarListToBuffer(1);
	  List[n1].data = *TempBuffer;
	  List[n1].grid = grid;
	  List[n1].proc = MyProcessorNumber;
	  n1++;
	  delete MoveStar;
	}  // ENDIF different processor

	// Same processor -- no move
	else {
	  InsertStarAfter(Stars, MoveStar);
	  NumberOfStars++;
	}

	i++;

      } // ENDWHILE stars
      
    } // ENDIF TotalToMove > PreviousTotalToMove

    delete [] ToGrid;
 
  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy stars back into grid. */
 
  else {

    /* Count up total number. */
 
    int TotalNumberOfStars;
    int NumberOfNewStars = EndIndex - StartIndex;

    TotalNumberOfStars = NumberOfStars + NumberOfNewStars;
 
    /* Copy stars from buffer into linked list */

    if (NumberOfNewStars > 0)
      for (i = StartIndex; i < EndIndex; i++) {
	MoveStar = StarBufferToList(List[i].data);
	MoveStar->GridID = this->ID;
	MoveStar->CurrentGrid = this;
	InsertStarAfter(this->Stars, MoveStar);
      } // ENDFOR stars
 
    /* Set new number of stars in this grid. */
 
    NumberOfStars = TotalNumberOfStars;

  } // end: if (COPY_IN)
 
  return SUCCESS;
}
