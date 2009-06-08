/*------------------------------------------------------------------------
  WRITE STREAMING DATA (Binary)
  By John Wise
  
  Created : 09 Sep 2004
  Modified: 

  Purpose : To walk through the hierarchy and call WriteNewMovieData
            for each grid.

  History : 
------------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"

/****************************** Prototypes ******************************/
int WriteStreamData(HierarchyEntry *Grids[], int NumberOfGrids,
		    TopGridData *MetaData, int CycleCount, int EndStep = FALSE) {

  int iGrid;

  for (iGrid = 0; iGrid < NumberOfGrids; iGrid++) {

    /*** Initialize the UNDER_SUBGRID_FIELD for this grid. ***/

    if (Grids[iGrid]->NextGridNextLevel != NULL) {

      Grids[iGrid]->GridData->
	ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
      
      /* Zero the solution (on this grid) which is underneath any
	 subgrid (so we get only the high resolution solution from the
	 subgrid). */
    
      HierarchyEntry *Temp2 = Grids[iGrid]->NextGridNextLevel;
      
      while (Temp2 != NULL) {
	Grids[iGrid]->GridData->
	  ZeroSolutionUnderSubgrid(Temp2->GridData, ZERO_UNDER_SUBGRID_FIELD);
	Temp2 = Temp2->NextGridThisLevel;
      }

    } /* END: Zero cells where children reside */

    /* Write data (if necessary : check inside routine) for this grid */
    if (Grids[iGrid]->GridData->
	WriteNewMovieData(MetaData->NewMovieLeftEdge,
			  MetaData->NewMovieRightEdge,
			  MetaData->StopTime,
			  EndStep, CycleCount) == FAIL) {
      fprintf(stderr, "Error in WriteNewMovie.\n");
      return FAIL;
    }

  } /* END: grid loop */

  return SUCCESS;

}

