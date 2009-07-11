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

void DeleteStarList(Star * &Node);
Star* StarBufferToList(StarBuffer buffer);
void InsertStarAfter(Star * &Node, Star * &NewNode);
 
int grid::CollectStars(int GridNum, int* &NumberToMove, 
		       int &StartIndex, int &EndIndex, 
		       star_data* &List, int CopyDirection)
{
 
  /* Declarations. */

  int i, j, dim, n1, grid, proc;
  Star *MoveStar;
  StarBuffer *TempBuffer;

  /* ----------------------------------------------------------------- */
  /* Copy star out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If there are no stars to move, we're done. */

    if (NumberOfStars == 0)
      return SUCCESS;

    /* If this is the correct processor, no copy-outs required. */

    if (MyProcessorNumber == ProcessorNumber)
      return SUCCESS;

    /* Add to the star count to move */

    // NumberOfStars is still the number of local stars, not the
    // actual total!
    NumberToMove[ProcessorNumber] += NumberOfStars;
 
    /* Move and delete stars */

    if (Stars == NULL)
      ENZO_FAIL("Star pointer cannot be NULL here.  NumberOfStars "
		"and star pointer are mismatched.");

    n1 = StartIndex;
    TempBuffer = Stars->StarListToBuffer(NumberOfStars);
    DeleteStarList(Stars);
    
    for (i = 0, n1 = StartIndex; i < NumberOfStars; i++, n1++) {
      List[n1].data = TempBuffer[i];
      List[n1].grid = GridNum;
      List[n1].proc = ProcessorNumber;
    } // ENDFOR stars

    StartIndex = n1;
    NumberOfStars = 0;

  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy stars back into grid. */
 
  else {

    if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* Count up total number. */
 
    int TotalNumberOfStars;
    int NumberOfNewStars = EndIndex - StartIndex;

    TotalNumberOfStars = NumberOfStars + NumberOfNewStars;

    if (NumberOfNewStars > 0)
      for (i = StartIndex; i < EndIndex; i++) {
	MoveStar = StarBufferToList(List[i].data);
	MoveStar->CurrentGrid = this;
	MoveStar->GridID = this->ID;
	InsertStarAfter(this->Stars, MoveStar);
      }
 
    /* Set new number of stars in this grid. */
 
    NumberOfStars = TotalNumberOfStars;
 
  } // end: if (COPY_IN)
 
  return SUCCESS;
}
