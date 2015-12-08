/***********************************************************************
/
/  COMPUTE STOCHASTIC FORCING OF THE FLOW
/
/  written by: W. Schmidt
/  date:       September 2005
/  modified1: Sep, 2014: updated to support Enzo 2.4 // P. Grete
/
/  PURPOSE: for each grid, the inverse FT is called to compute 
/           the pyhsical force field from the forcing spectrum
/
************************************************************************/
#include "preincludes.h"
#include <iostream>
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
using namespace std;


int ComputeStochasticForcing(TopGridData *MetaData, 
                 HierarchyEntry *Grids[], int NumberOfGrids)
{
  int grid;

  /* Compute vector components of the force field in succession. */

  for (grid = 0; grid < NumberOfGrids; grid++) {
      if (debug) cout << "ComputeStochasticForcing: computing phase factors for grid " << grid << endl;
      Grids[grid]->GridData->Phases();
      if (debug) cout << "ComputeStochasticForcing: computing force field for grid " << grid << endl;
      for (int dim = 0; dim < MetaData->TopGridRank; dim++)
      if (Grids[grid]->GridData->FTStochasticForcing(dim) == FAIL) {
          fprintf(stderr, "Error in grid->FTStochasticForcing\n");
          return FAIL;
      }
  }

  return SUCCESS;
}
