/***********************************************************************
/
/  SET THE TIMESTEP FOR A LEVEL
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Matthew Turk, split off
/  date:       June 2009
/
/  PURPOSE:
/       Determine the timestep for this iteration of the loop.
/
************************************************************************/
 
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

 
float CommunicationMinValue(float Value);

int SetLevelTimeStep(HierarchyEntry *Grids[], int NumberOfGrids, int level,
		     float *dtThisLevelSoFar, float *dtThisLevel,
		     float dtLevelAbove)
{
  float dtGrid, dtActual, dtLimit;
  int grid1;

  LCAPERF_START("SetLevelTimeStep"); // SetTimeStep()

  if (level == 0) {
 
    /* For root level, use dtLevelAbove. */
 
    *dtThisLevel      = dtLevelAbove;
    *dtThisLevelSoFar = dtLevelAbove;
    dtActual          = dtLevelAbove;
 
  } else {
 
    /* Compute the mininum timestep for all grids. */
 
    *dtThisLevel = huge_number;
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      dtGrid      = Grids[grid1]->GridData->ComputeTimeStep();
      *dtThisLevel = min(*dtThisLevel, dtGrid);
    }
    *dtThisLevel = CommunicationMinValue(*dtThisLevel);

    dtActual = *dtThisLevel;

    /* Compute timestep without conduction and use to set the number 
       of iterations without rebuiding the hierarchy. */

    if (ConductionDynamicRebuildHierarchy && LevelSubCycleCount[level] == 0) {
      float dt_no_conduction, my_dt_no_conduction;
      int my_isotropic_conduction = IsotropicConduction;
      int my_anisotropic_conduction = AnisotropicConduction;
      IsotropicConduction = AnisotropicConduction = FALSE;

      dt_no_conduction = huge_number;
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
        my_dt_no_conduction = Grids[grid1]->GridData->ComputeTimeStep();
        dt_no_conduction = min(dt_no_conduction, my_dt_no_conduction);
      }
      dt_no_conduction = CommunicationMinValue(dt_no_conduction);
      int my_cycle_skip = (int) (dt_no_conduction / dtActual);
      RebuildHierarchyCycleSkip[level] = max(1, my_cycle_skip);
      if (debug) fprintf(stderr, "AdaptiveRebuildHierarchyCycleSkip[%"ISYM"] = %"ISYM"\n",
                         level, RebuildHierarchyCycleSkip[level]);

      IsotropicConduction = my_isotropic_conduction;
      AnisotropicConduction = my_anisotropic_conduction;
    }

#ifdef USE_DT_LIMIT

    //    dtLimit = LevelZeroDeltaT/(4.0)/POW(RefineBy,level);

    dtLimit = 0.5/(4.0)/POW(2.0,level);

    if ( dtActual < dtLimit )
      *dtThisLevel = dtLimit;

#endif
 
    /* Advance dtThisLevelSoFar (don't go over dtLevelAbove). */
 
    if (*dtThisLevelSoFar+*dtThisLevel*1.05 >= dtLevelAbove) {
      *dtThisLevel      = dtLevelAbove - *dtThisLevelSoFar;
      *dtThisLevelSoFar = dtLevelAbove;
    }
    else
      *dtThisLevelSoFar += *dtThisLevel;
 
  }

  if (debug) 
    printf("Level[%"ISYM"]: dt = %"GSYM"  %"GSYM" (%"GSYM"/%"GSYM")\n", 
	   level, *dtThisLevel, dtActual, *dtThisLevelSoFar, dtLevelAbove);
 
  /* Set all grid's timestep to this minimum dt. */
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->SetTimeStep(*dtThisLevel);
 

  LCAPERF_STOP("SetLevelTimeStep"); // SetTimeStep()
  return SUCCESS;
}
