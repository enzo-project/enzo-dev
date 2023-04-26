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

#ifdef USE_LIBYT
#include "libyt.h"
#include "libyt_interactive_mode.h"
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
#include "CosmologyParameters.h"
#include "TopGridData.h"


int  GetUnits(float *DensityUnits, float *LengthUnits,
		       float *TemperatureUnits, float *TimeUnits,
		       float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int ExposeHierarchyToLibyt(TopGridData *MetaData, HierarchyEntry *Grid, int
        &GridID, int &LocalGridID, FLOAT WriteTime, int ParentID, int level,
        yt_grid *GridInfoArray);
void ExposeGridHierarchy(int NumberOfGrids);
void ExportParameterFile(TopGridData *MetaData, FLOAT CurrentTime, FLOAT OldTime, float dtFixed);
void CommunicationBarrier();

int CallInSitulibyt(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level, int from_topgrid)
{
#ifndef USE_LIBYT
    return SUCCESS;
#else

  yt_param_yt *params = (yt_param_yt*) param_yt;
  int i, j;


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
    int num_grids, num_local_grids, start_index;
    num_grids = num_local_grids = 0; start_index = 1;

    LevelHierarchyEntry *Temp2 = LevelArray[0];
    /* Count the grids */
    /* I think there is a better idiom for this somewhere
       but I couldn't find it, and I think this works     */
    for (int lc = 0; LevelArray[lc] != NULL; lc++){
        Temp2 = LevelArray[lc];
        while (Temp2 != NULL) {
            num_grids++;
            if (Temp2->GridData->ReturnProcessorNumber() == MyProcessorNumber)
                num_local_grids++;
            Temp2 = Temp2->NextGridThisLevel;
        }
    }

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
    params->frontend = "enzo";
    params->fig_basename = "Fig";
    params->domain_left_edge[0] = (double) DomainLeftEdge[0];
    params->domain_left_edge[1] = (double) DomainLeftEdge[1];
    params->domain_left_edge[2] = (double) DomainLeftEdge[2];
    params->domain_right_edge[0] = (double) DomainRightEdge[0];
    params->domain_right_edge[1] = (double) DomainRightEdge[1];
    params->domain_right_edge[2] = (double) DomainRightEdge[2];

    params->current_time = CurrentTime;

    if (ComovingCoordinates) {
      FLOAT a, dadt, FinalRedshift, CurrentRedshift;
      CosmologyComputeExpansionFactor(MetaData->StopTime, &a, &dadt);

      FinalRedshift = (1 + InitialRedshift)/a - 1;

      /* Compute the current redshift (for information only). */

      CosmologyComputeExpansionFactor(CurrentTime, &a, &dadt);
      CurrentRedshift = (1 + InitialRedshift)/a - 1;

      params->cosmological_simulation = 1;
      params->current_redshift = CurrentRedshift;
      params->omega_lambda = OmegaLambdaNow;
      params->omega_matter = OmegaMatterNow;
      params->hubble_constant = HubbleConstantNow;

      /* We will need to add a bunch of additional parameters later. */

    } else {
      params->cosmological_simulation = 0;
    }

    params->length_unit = LengthUnits;
    params->mass_unit = DensityUnits * LengthUnits * LengthUnits * LengthUnits; /* this right? */
    params->time_unit = TimeUnits;
    params->magnetic_unit = 0.0; /* Not right */
    params->periodicity[0] = 1; /* also not right */
    params->periodicity[1] = 1;
    params->periodicity[2] = 1;
    params->dimensionality = MetaData->TopGridRank;
    params->domain_dimensions[0] = MetaData->TopGridDims[0];
    params->domain_dimensions[1] = MetaData->TopGridDims[1];
    params->domain_dimensions[2] = MetaData->TopGridDims[2];
    params->refine_by = RefineBy;
    params->num_grids = num_grids;
    params->num_grids_local = num_local_grids;
    /* We do things by DataLabel */

    for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
        if (!DataLabel[i]) {
            params->num_fields = i;
            break;
        }
    }

    /* TODO: set num_par_types and par_type_list, so that libyt can properly initialize
     *       particle related stuff, etc yt_get_ParticlesPtr gets the pointer,
     *       and Grid_ConvertToLibyt.C have par_count_list initialized. */
    params->num_par_types = 0;

    if (yt_set_Parameters(params) != YT_SUCCESS){
        fprintf(stderr, "Error in libyt API yt_set_Parameters\n");
        return FAIL;
    }

    /* Set code-specific parameter
     * TODO: including this reach seg fault, but success to compile it.
     *       I think the seg fault comes from libyt and from yt_set_Parameters */
    char tempname[256];
    #include "InitializeLibytInterface_finderfunctions.inc"

    /* Here, we have a delicate operation to conduct.  We are setting up the fields
     * supplied to libyt.  The issue we need to be wary of is that we are setting them
     * up in the order they are in DataLabel, which *may* not be the same as in
     * the grids.  (We can't know for sure!)
     * */

    yt_field *field_list;
    yt_get_FieldsPtr( &field_list );

    int libyt_field_i = 0;
    for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
        /* This will be out of date when/if a new field is added to DataLables.
         *
         * We could potentially be much more careful about this, but looking at
         * the logic in grid::WriteGrid it is clear that the assumption, for IO
         * purposes, that DataLabel[f] maps to BaryonField[f], is implicit
         * throughout many relevant places.
         *
         * Note that this does not account for some other potential fields,
         * such as Temperature, that are not stored in DataLabel, which we will
         * address by hand.
         *
         * */
        if (!DataLabel[i]) break;
        /* This tells us that BaryonFields[i] maps to
         * libyt_field[libyt_field_lookup[i]]
         *
         * */

        field_list[libyt_field_i].field_name = DataLabel[i];
        field_list[libyt_field_i].field_type = "cell-centered";
        field_list[libyt_field_i].field_dtype = EYT_FLOAT;
        for (j = 0; j < 6; j++) {
            /*
             * It may be possible that in some cases, this global value is not
             * correct; however, it's pretty unlikely, and non-cell-centered
             * fields will be stored in different member fields anyway.
             *
             * */
            field_list[libyt_field_i].field_ghost_cell[j] = NumberOfGhostZones;
        }
        libyt_field_lookup[i] = libyt_field_i++;
    }

    /* Now we add on the following fields, as per grid::WriteGrid:
     *
     *  - Temperature
     *  - Dust_Temperature
     *  - Cooling_Time
     *
     *  Each of these is predicated on the global parameter associated with
     *  them.
     *
     * */

    /* We now have to do everything we do in CallPython.C, which amounts to
     *
     *  - ExposeGridHierarchy (not necessary anymore)
     *  - ExposeDataHierarchy (a recursive call)
     *
     * */

    /*
     * TODO: call yt_get_ParticlesPtr for num_par_types > 0
     */

    /* As I'm writing this in 2023, I cannot recall specifically why we need to
     * find the last entry in the linked list.  But, it shows up a number of
     * places, including in CallPython.C and OutputFromEvolveLevel.C.  I
     * suspect it's because the grids are added in the reverse order to the
     * LevelArray.
     *
     * I have left the original 'ugh:' so that it can be found by anyone
     * grepping for that.  I don't actually feel 'ugh' about it.  In fact,
     * in revisiting Enzo, I have found a number of places where in yt we
     * took 'the easy way' and the Enzo way is so much more elegant, even if
     * the impedance mismatch means I have to jump through a few more hoops.
     *
     * */

    yt_grid* GridInfoArray;
    yt_get_GridsPtr(&GridInfoArray);

    /* These are 1-indexed, so when we access elements we -1 them */
    int GlobalGridID = 1;
    int LocalGridID = 1;
    Temp2 = LevelArray[0];
    while (Temp2->NextGridThisLevel != NULL)
        Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
    if(ExposeHierarchyToLibyt(MetaData, Temp2->GridHierarchyEntry, 
                GlobalGridID, LocalGridID, CurrentTime, 0, 0, GridInfoArray) == FAIL) {
        fprintf(stderr, "Error in ExposeHierarchyToLibyt\n");
        return FAIL;
    }

    /* Commit all settings to libyt. */
    if (yt_commit() != YT_SUCCESS){
        fprintf(stderr, "Error in libyt API yt_commit\n");
        return FAIL;
    }

    /* TODO: Call another other inline Python Function using
     *       yt_run_Function and yt_run_FunctionArguments... */

    /* Call interactive mode. */
    if (yt_run_InteractiveMode("LIBYT_STOP") != YT_SUCCESS) {
        fprintf(stderr, "Error in libyt API yt_run_InteractiveMode\n");
        return FAIL;
    }

    /* Free resources allocated for libyt. */
    if (yt_free() != YT_SUCCESS) {
        fprintf(stderr, "Error in libyt API yt_free\n");
        return FAIL;
    }

    CommunicationBarrier();
    return SUCCESS;
#endif
}
