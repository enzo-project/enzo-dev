/***********************************************************************
/
/  MAGNETIC FIELD RESETTER
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/  modified1:
/
/  PURPOSE: Simply resets the magnetic field in the entire simulation.
/           Will be useful when one wants to restart or resimulate 
/           the dump with nonzero magnetic field imposed.
/
************************************************************************/

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

int MagneticFieldResetter(LevelHierarchyEntry *LevelArray[], int ThisLevel,
			  TopGridData *MetaData)
{

  /* Return if this does not concern us */

  if (ResetMagneticField != TRUE || 
      (HydroMethod != MHD_RK && HydroMethod != HD_RK) ||
      !(MetaData->FirstTimestepAfterRestart)) 
    return SUCCESS;

  int level, i, grid1;
  HierarchyEntry **Grids;
  int NumberOfGrids;
  
  fprintf(stdout, "Resetting magnetic fields at the start.\n");

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
    
//  fprintf(stdout, "Resetting magnetic fields at level = %d ...\n", level);

    NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
    
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      
      if (Grids[grid1]->GridData->MagneticFieldResetter(level) == FAIL) {
	ENZO_FAIL("Error in grid::MagneticFieldResetter.\n");

      }
      
    }  // loop for grid1
    
    CommunicationBarrier(); 

  }  // loop for level

  /* Set resetter parameter zero; otherwise it will reset the field at next restart */

  ResetMagneticField = 0;

  return SUCCESS;

}







