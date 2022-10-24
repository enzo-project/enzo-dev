/***********************************************************************
/
/  CALL LIBYT AT FIXED TIMESTEPS
/
/  written by: Matthew Turk
/  date:       October, 2022
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

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
#include "CosmologyParameters.h"
#include "TopGridData.h"


int  GetUnits(float *DensityUnits, float *LengthUnits,
		       float *TemperatureUnits, float *TimeUnits,
		       float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int ExposeDataHierarchy(TopGridData *MetaData, HierarchyEntry *Grid, 
		       int &GridID, FLOAT WriteTime, int reset, int ParentID, int level);
void ExposeGridHierarchy(int NumberOfGrids);
void ExportParameterFile(TopGridData *MetaData, FLOAT CurrentTime, FLOAT OldTime, float dtFixed);
void CommunicationBarrier();

int CallInSitulibyt(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level, int from_topgrid)
{
#ifndef USE_LIBYT
    return SUCCESS;
#else

  yt_param_yt param_yt;


/*
    NumberOfPythonCalls++;
    if (from_topgrid) {
      NumberOfPythonTopGridCalls++;
      if (!(PythonTopGridSkip) ||
	  (NumberOfPythonTopGridCalls % PythonTopGridSkip) != 0) return SUCCESS;
    }
    else {
      if (LevelArray[level+1] != NULL) return SUCCESS;
      NumberOfPythonSubcycleCalls++;
      if (!(PythonSubcycleSkip) ||
	  (NumberOfPythonSubcycleCalls % PythonSubcycleSkip) != 0) return SUCCESS;
    }
*/

    FLOAT CurrentTime, OldTime;
    float dtFixed;
    int num_grids, start_index;
    num_grids = 0; start_index = 1;

    LevelHierarchyEntry *Temp2 = LevelArray[0];
    /* Count the grids */
    /* I think there is a better idiom for this somewhere
       but I couldn't find it, and I think this works     */
    for (int lc = 0; LevelArray[lc] != NULL; lc++){
        Temp2 = LevelArray[lc];
        while (Temp2 != NULL) {
            num_grids++; Temp2 = Temp2->NextGridThisLevel;
        }
    }
    ExposeGridHierarchy(num_grids);
    Temp2 = LevelArray[0];
    while (Temp2->NextGridThisLevel != NULL)
        Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
    CurrentTime = LevelArray[level]->GridData->ReturnTime();
    OldTime = LevelArray[level]->GridData->ReturnOldTime();
    dtFixed = LevelArray[level]->GridData->ReturnTimeStep();
    /*
    if (ExposeDataHierarchy(MetaData, Temp2->GridHierarchyEntry, start_index,
                CurrentTime, 1, 0, 0) == FAIL) {
        fprintf(stderr, "Error in ExposeDataHierarchy\n");
        return FAIL;
    }
    */

    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
          VelocityUnits = 1;

    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, CurrentTime);


    /* ExportParameterFile(MetaData, CurrentTime, OldTime, dtFixed); */
    param_yt.frontend = "enzo";
    param_yt.domain_left_edge[0] = (double) DomainLeftEdge[0];
    param_yt.domain_left_edge[1] = (double) DomainLeftEdge[1];
    param_yt.domain_left_edge[2] = (double) DomainLeftEdge[2];
    param_yt.domain_right_edge[0] = (double) DomainRightEdge[0];
    param_yt.domain_right_edge[1] = (double) DomainRightEdge[1];
    param_yt.domain_right_edge[2] = (double) DomainRightEdge[2];

    param_yt.current_time = CurrentTime;

    if (ComovingCoordinates) {
      FLOAT a, dadt, FinalRedshift, CurrentRedshift;
      CosmologyComputeExpansionFactor(MetaData->StopTime, &a, &dadt);

      FinalRedshift = (1 + InitialRedshift)/a - 1;

      /* Compute the current redshift (for information only). */

      CosmologyComputeExpansionFactor(CurrentTime, &a, &dadt);
      CurrentRedshift = (1 + InitialRedshift)/a - 1;

      param_yt.cosmological_simulation = 1;
      param_yt.current_redshift = CurrentRedshift;
      param_yt.omega_lambda = OmegaLambdaNow;
      param_yt.omega_matter = OmegaMatterNow;
      param_yt.hubble_constant = HubbleConstantNow;

      /* We will need to add a bunch of additional parameters later. */

    } else {
      param_yt.cosmological_simulation = 0;
    }

    param_yt.length_unit = LengthUnits;
    param_yt.mass_unit = DensityUnits * LengthUnits * LengthUnits * LengthUnits; /* this right? */
    param_yt.time_unit = TimeUnits;
    param_yt.magnetic_unit = 0.0; /* Not right */
    param_yt.periodicity[0] = 1; /* also not right */
    param_yt.periodicity[1] = 1;
    param_yt.periodicity[2] = 1;
    param_yt.dimensionality = MetaData->TopGridRank;
    param_yt.domain_dimensions[0] = MetaData->TopGridDims[0];
    param_yt.domain_dimensions[1] = MetaData->TopGridDims[1];
    param_yt.domain_dimensions[2] = MetaData->TopGridDims[2];
    param_yt.refine_by = RefineBy;
    /*param_yt.num_grids;
    param_yt.num_fields;
    param_yt.num_species;*/


    CommunicationBarrier();
    return SUCCESS;
#endif
}
