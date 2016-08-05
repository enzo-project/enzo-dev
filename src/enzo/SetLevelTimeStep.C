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
 
    /* Calculate timestep without conduction and get conduction separately later. */

    int my_isotropic_conduction = IsotropicConduction;
    int my_anisotropic_conduction = AnisotropicConduction;
    int dynamic_hierarchy_rebuild = ConductionDynamicRebuildHierarchy &&
      level >= ConductionDynamicRebuildMinLevel &&
      dtRebuildHierarchy[level] <= 0.0;

    if (dynamic_hierarchy_rebuild) {
      IsotropicConduction = AnisotropicConduction = FALSE;      
    }

    /* Compute the mininum timestep for all grids. */
 
    *dtThisLevel = huge_number;
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      dtGrid      = Grids[grid1]->GridData->ComputeTimeStep();
      *dtThisLevel = min(*dtThisLevel, dtGrid);
    }
    *dtThisLevel = CommunicationMinValue(*dtThisLevel);

    /* Compute conduction timestep and use to set the number 
       of iterations without rebuiding the hierarchy. */

    if (dynamic_hierarchy_rebuild) {

      /* Return conduction parameters to original values. */
      IsotropicConduction = my_isotropic_conduction;
      AnisotropicConduction = my_anisotropic_conduction;

      float dt_conduction;
      float dt_cond_temp;
      dt_conduction = huge_number;
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
        if (Grids[grid1]->GridData->ComputeConductionTimeStep(dt_cond_temp) == FAIL) 
          ENZO_FAIL("Error in ComputeConductionTimeStep.\n");
	dt_conduction = min(dt_conduction,dt_cond_temp);
      }
      dt_conduction = CommunicationMinValue(dt_conduction);
      dt_conduction *= float(NumberOfGhostZones);  // for subcycling

      int my_cycle_skip = max(1, (int) (*dtThisLevel / dt_conduction));
      dtRebuildHierarchy[level] = *dtThisLevel;
      if (debug)
        fprintf(stderr, "Conduction dt[%"ISYM"] = %"GSYM", will rebuild hierarchy in about %"ISYM" cycles.\n",
                         level, dt_conduction, my_cycle_skip);

      /* Set actual timestep correctly. */
      *dtThisLevel = min(*dtThisLevel, dt_conduction);

    }

    dtActual = *dtThisLevel;

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
