/***********************************************************************
/
/  GRID CLASS (MOVE ALL STARS FROM SPECIFIED GRID TO THIS GRID)
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE:
/
/    NOTE: Adapted from grid::MoveAllParticles()
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void DeleteStarList(Star * &Node);
void InsertStarAfter(Star * &Node, Star * &NewNode);
Star *PopStar(Star * &Node);
Star* StarBufferToList(StarBuffer *buffer, int n);

int grid::MoveAllStars(int NumberOfGrids, grid* FromGrid[], int TopGridDimension)
{

  StarBuffer *buffer;
  Star *NewStar, *cstar;

  if (NumberOfGrids < 1) {
    fprintf(stderr, "NumberOfGrids(%"ISYM") must be > 0.\n", NumberOfGrids);
    return FAIL;
  }

  /* Determine total number of stars. */

  int TotalNumberOfStars = NumberOfStars;
  int i, j, grid, dim;

  for (grid = 0; grid < NumberOfGrids; grid++)
    TotalNumberOfStars += FromGrid[grid]->NumberOfStars;
  if (TotalNumberOfStars == 0)
    return SUCCESS;

  /* Determine level of this grid */

  int ThisLevel = nint(-logf(TopGridDimension * CellWidth[0][0]) / logf(2.0f));

  /* Debugging info. */

//  if (debug) printf("MoveAllStars: %"ISYM" (before: ThisGrid = %"ISYM").\n",
//		    TotalNumberOfStars, NumberOfStars);

  /* Copy FromGrids' stars to new grid. */

  for (grid = 0; grid < NumberOfGrids; grid++) {

   /* If on the same processor, just copy. */

    if (MyProcessorNumber == ProcessorNumber &&
        MyProcessorNumber == FromGrid[grid]->ProcessorNumber) {

      if (FromGrid[grid]->Stars != NULL) {
	cstar = FromGrid[grid]->Stars;
	while (cstar != NULL) {
	  NewStar = PopStar(cstar); // also advances cstar pointer
	  NewStar->CurrentGrid = this;
	  NewStar->level = ThisLevel;
	  NewStar->GridID = this->ID;
	  InsertStarAfter(this->Stars, NewStar);
	  //cstar = cstar->NextStar;
	} // ENDWHILE stars
	FromGrid[grid]->Stars = NULL;
      }

    }

    /* Otherwise, communicate. */

    else {
      if (MyProcessorNumber == ProcessorNumber ||
          MyProcessorNumber == FromGrid[grid]->ProcessorNumber)
	if (FromGrid[grid]->CommunicationSendStars(this, ProcessorNumber) == FAIL) {
	  fprintf(stderr, "Error in grid->CommunicationSendStars.\n");
	  return FAIL;
        }

    } // ENDELSE same processor

  } // end: loop over grids.

  /* Set new number of particles in this grid. */

  NumberOfStars = TotalNumberOfStars;

  /* Set number of particles to zero.  No need to delete particles
     here, which is done in CommunicationSendStars.  For local copy,
     we copy the pointer to the new grid, and we shouldn't delete the
     linked list.  */

  for (grid = 0; grid < NumberOfGrids; grid++) {
    FromGrid[grid]->NumberOfStars = 0;
//    if (MyProcessorNumber == FromGrid[grid]->ProcessorNumber)
//      DeleteStarList(FromGrid[grid]->Stars);
  }

  return SUCCESS;
}
