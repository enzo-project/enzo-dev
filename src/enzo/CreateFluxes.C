
//
// CreateFluxes
// David Collins, June 11 2009
// Generate the subgrid fluxes for the given level
//


 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"


int CreateFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		 int NumberOfGrids,int NumberOfSubgrids[]){


    /* For each grid, create the subgrid list. */


    LCAPERF_START("CreateFluxes"); // create subgrid list
    int grid1, counter,subgrid;
    HierarchyEntry *NextGrid;
    int RefinementFactors[MAX_DIMENSION];


    /* For each grid, compute the number of it's subgrids. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      NextGrid = Grids[grid1]->NextGridNextLevel;
      counter = 0;
      while (NextGrid != NULL) {
	NextGrid = NextGrid->NextGridThisLevel;
	if (++counter > MAX_NUMBER_OF_SUBGRIDS) {
	  ENZO_FAIL("More subgrids than MAX_NUMBER_OF_SUBGRIDS.\n");
	}
      }
      NumberOfSubgrids[grid1] = counter + 1;
    }


    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      /* Allocate the subgrid fluxes for this grid. */
 
      SubgridFluxesEstimate[grid1] = new fluxes *[NumberOfSubgrids[grid1]];
 
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid1]; subgrid++)
	SubgridFluxesEstimate[grid1][subgrid] = NULL;
 
      /* Collect the flux data and store it in the newly minted fluxes.
	 Or rather that's what we should do.  Instead, we create fluxes one
	 by one in this awkward array of pointers to pointers.  This should be
	 changed so that all the routines take arrays of flux rather than
	 arrays of pointers to flux.  Dumb. */
 
      counter = 0;

      // Only allocate fluxes for local grids: saves a *lot* of storage
      
      if (MyProcessorNumber ==

          Grids[grid1]->GridData->ReturnProcessorNumber()) {
	
	NextGrid = Grids[grid1]->NextGridNextLevel;
	while (NextGrid != NULL) {
	  SubgridFluxesEstimate[grid1][counter] = new fluxes;
	  Grids[grid1]->GridData->ComputeRefinementFactors
	    (NextGrid->GridData, RefinementFactors);
	  NextGrid->GridData->ReturnFluxDims
	    (*(SubgridFluxesEstimate[grid1][counter++]), RefinementFactors);
	  NextGrid = NextGrid->NextGridThisLevel;
	}
	
	/* Add the external boundary of this subgrid to the subgrid list. This
	   makes it easy to keep adding up the fluxes of this grid, but we must
	   keep in mind that the last subgrid should be ignored elsewhere. */
	
	SubgridFluxesEstimate[grid1][counter] = new fluxes;
	Grids[grid1]->GridData->ComputeRefinementFactors
	  (Grids[grid1]->GridData, RefinementFactors);
	Grids[grid1]->GridData->ReturnFluxDims
	  (*(SubgridFluxesEstimate[grid1][counter]), RefinementFactors);
	
      }
      
    } // end loop over grids (create Subgrid list)
    // dcc flux cut stop
    LCAPERF_STOP("CreateFluxes"); // create subgrid list
    return SUCCESS;
}
