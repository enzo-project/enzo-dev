 
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

void DeleteFluxes(fluxes *Fluxes);
int FinalizeFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		   int NumberOfGrids,int NumberOfSubgrids[]){
  
  int grid1,subgrid;
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
    
    // Only deallocate fluxes for local grids
    
    if (MyProcessorNumber ==
	Grids[grid1]->GridData->ReturnProcessorNumber()) {
      
      if (FluxCorrection)
	if (Grids[grid1]->GridData->AddToBoundaryFluxes
	    (SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1])
	    == FAIL) {
	  ENZO_FAIL("Error in grid->AddToBoundaryFluxes.\n");

	}
      
      /* Delete fluxes pointed to by SubgridFluxesEstimate[subgrid]. */
      
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid1]; subgrid++) {
	DeleteFluxes(SubgridFluxesEstimate[grid1][subgrid]);
	delete       SubgridFluxesEstimate[grid1][subgrid];
      }
      delete [] SubgridFluxesEstimate[grid1];
      
    }
    
  } // end of loop over grids
  return SUCCESS;
}
