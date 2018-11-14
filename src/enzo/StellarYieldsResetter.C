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
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#include "CommunicationUtilities.h"

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
                      HierarchyEntry **Grids[]);

int StellarYieldsResetter(LevelHierarchyEntry *LevelArray[], int ThisLevel,
                          TopGridData *MetaData){


  int level, i, grid1;
  HierarchyEntry **Grids;
  int NumberOfGrids;

  fprintf(stdout, "Resetting certain stellar abundance fields to zero\n");

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++){

    NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++){

      if (Grids[grid1]->GridData->StellarYieldsResetter(level) == FAIL){
        ENZO_FAIL("Error in grid::StellarYieldsResetter.\n");
      }
    }
    CommunicationBarrier();
  }


 ResetStellarAbundances = 0; // set to zero - make this manual ONLY

 return SUCCESS;
}
