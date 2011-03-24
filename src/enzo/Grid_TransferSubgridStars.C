/***********************************************************************
/
/  GRID CLASS (COPY SUBGRID STARS INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/  modified2:  May, 2009 by John Wise: modified to move subgrid particles
/  modified2:  July, 2009 by John Wise: modified to move stars
/
/  PURPOSE:
/
************************************************************************/
 
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
 
int grid::TransferSubgridStars(grid* Subgrids[], int NumberOfSubgrids, 
			       int* &NumberToMove, int StartIndex, 
			       int EndIndex, star_data* &List, 
			       bool KeepLocal, bool ParticlesAreLocal,
			       int CopyDirection, int IncludeGhostZones)
{
 
  /* Declarations. */

  int i, j, index, dim, n1, grid, proc;
  int i0, j0, k0;
  Star *cstar, *MoveStar;
  StarBuffer *TempBuffer;

  /* ----------------------------------------------------------------- */
  /* Copy stars out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If stars aren't distributed over several processors, exit
       if this isn't the host processor. */

    if (ParticlesAreLocal && MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* If there are no stars to move, we're done. */

    if (NumberOfStars == 0)
      return SUCCESS;

    /* Set boundaries (with and without ghost zones) */

    int StartIndex[] = {1,1,1}, EndIndex[] = {1,1,1};
    if (IncludeGhostZones)
      for (dim = 0; dim < GridRank; dim++) {
	StartIndex[dim] = 0;
	EndIndex[dim] = GridDimension[dim]-1;
      }
    else
      for (dim = 0; dim < GridRank; dim++) {
	StartIndex[dim] = GridStartIndex[dim];
	EndIndex[dim] = GridEndIndex[dim];
      }
 
    /* Count the number of stars already moved */

    int PreviousTotalToMove = 0;
    for (i = 0; i < NumberOfProcessors; i++)
      PreviousTotalToMove += NumberToMove[i];
 
    /* Count stars to move */

    
    int *subgrid = NULL;
    subgrid = new int[NumberOfStars];

    for (i = 0, cstar = Stars; cstar; i++, cstar = cstar->NextStar) {

      /* Compute index of star position. */
 
      i0 = int((cstar->pos[0] - CellLeftEdge[0][0])/CellWidth[0][0]);
      if (GridRank > 0)
       j0 = int((cstar->pos[1] - CellLeftEdge[1][0])/CellWidth[1][0]);
      if (GridRank > 1)
       k0 = int((cstar->pos[2] - CellLeftEdge[2][0])/CellWidth[2][0]);
 
      i0 = max(min(EndIndex[0], i0), StartIndex[0]);
      j0 = max(min(EndIndex[1], j0), StartIndex[1]);
      k0 = max(min(EndIndex[2], k0), StartIndex[2]);
 
      index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;
 
      /* Find and store subgrid number of this star, and add to
	 count. */
 
      subgrid[i] = nint(BaryonField[NumberOfBaryonFields][index])-1;
      if (subgrid[i] >= 0) {
	if (KeepLocal)
	  proc = MyProcessorNumber;
	else
	  proc = Subgrids[subgrid[i]]->ReturnProcessorNumber();
	NumberToMove[proc]++;
      }
      if (subgrid[i] < -1 || subgrid[i] > NumberOfSubgrids-1) {
	ENZO_VFAIL("star subgrid (%"ISYM"/%"ISYM") out of range\n", 
		subgrid[i], NumberOfSubgrids)
      }
      
    } // ENDFOR stars

    /* Allocate space. */
 
    int TotalToMove = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++)
      TotalToMove += NumberToMove[proc];
 
    if (TotalToMove > PreviousTotalToMove) {

      // Increase the size of the list to include the stars from
      // this grid

      star_data *NewList = new star_data[TotalToMove];
      memcpy(NewList, List, PreviousTotalToMove * sizeof(star_data));
      delete [] List;
      List = NewList;

      /* Move stars */

      n1 = PreviousTotalToMove;
      NumberOfStars = 0;
      cstar = Stars;
      Stars = NULL;
      i = 0;
      while (cstar != NULL) {

	MoveStar = PopStar(cstar); // also advances to NextStar

	if (subgrid[i] >= 0) {
	  TempBuffer = MoveStar->StarListToBuffer(1);
	  delete MoveStar;

	  List[n1].data = *TempBuffer;
	  List[n1].grid = subgrid[i];
	  List[n1].proc = (KeepLocal) ? MyProcessorNumber :
	    Subgrids[subgrid[i]]->ReturnProcessorNumber();
	  n1++;

	} // ENDIF move to subgrid

	else {
	  InsertStarAfter(Stars, MoveStar);
	  NumberOfStars++;
	}

	i++;

      } // ENDWHILE stars

    } // ENDIF stars to move

    delete [] subgrid;
 
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
	MoveStar->GridID = List[i].grid;
	MoveStar->CurrentGrid = this;
	// FALSE if going to a subgrid
	if (IncludeGhostZones == FALSE)
	  MoveStar->IncreaseLevel();
	InsertStarAfter(this->Stars, MoveStar);
      } // ENDFOR stars
 
    /* Set new number of stars in this grid. */
 
    NumberOfStars = TotalNumberOfStars;
  
  } // end: if (COPY_IN)

 
  return SUCCESS;
}
