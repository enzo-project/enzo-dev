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


int ComputeDomainBoundaryMassFlux(HierarchyEntry *Grids[], int level, int NumberOfGrids){

  if (!StoreDomainBoundaryMassFlux) return SUCCESS;

  if (level != 0) return SUCCESS; // only do for root grid

  int grid1;

  float allgrid_BoundaryMassFluxContainer[MAX_NUMBER_OF_BARYON_FIELDS];
  for (int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i ++) allgrid_BoundaryMassFluxContainer[i] = 0.0;

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++){

    if( Grids[grid1]->GridData->ComputeDomainBoundaryMassFlux(allgrid_BoundaryMassFluxContainer) == FAIL)
      ENZO_FAIL("Error in grid->ComputeDomainBoundaryMassFlux.\n");
  }

  /* now communicate */
  CommunicationSumValues(allgrid_BoundaryMassFluxContainer, MAX_NUMBER_OF_BARYON_FIELDS);

  /* now store total in root grid for output */
  if (MyProcessorNumber == ROOT_PROCESSOR){
    for( int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i ++){
        BoundaryMassFluxContainer[i] += max(allgrid_BoundaryMassFluxContainer[i], 0.0);
    }
  }

  return SUCCESS;
}
