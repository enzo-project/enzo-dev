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
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

#ifdef USE_LIBYT

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

static yt_dtype MapHDF5TypeToYTType(hid_t hdf5type) {
    if (hdf5type == HDF5_INT) {
        return EYT_INT;
    } else if (hdf5type == HDF5_REAL) {
        return EYT_BFLOAT;
    } else if (hdf5type == HDF5_R8) {
        return YT_DOUBLE;
    } else if (hdf5type == HDF5_PREC) {
        return EYT_PFLOAT;
    } else {
        return YT_DTYPE_UNKNOWN;
    }
}

#endif

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
    params->velocity_unit = VelocityUnits;
    params->magnetic_unit = sqrt(4.0 * 3.141592653589793238462643383279502884L * DensityUnits) * VelocityUnits;
    params->periodicity[0] = MetaData->LeftFaceBoundaryCondition[0];
    params->periodicity[1] = MetaData->LeftFaceBoundaryCondition[1];
    params->periodicity[2] = MetaData->LeftFaceBoundaryCondition[2];
    params->dimensionality = MetaData->TopGridRank;
    params->domain_dimensions[0] = MetaData->TopGridDims[0];
    params->domain_dimensions[1] = MetaData->TopGridDims[1];
    params->domain_dimensions[2] = MetaData->TopGridDims[2];
    params->refine_by = RefineBy;
    params->index_offset = 1;
    params->num_grids = num_grids;
    params->num_grids_local = num_local_grids;
    /* We do things by DataLabel */

    for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
        if (!DataLabel[i]) {
            params->num_fields = i;
            break;
        }
    }

    // Add fields for Temperature/Cooling_Time derived field from enzo
    params->num_fields += 2;

    // Add active particle ptypes
    params->num_par_types = 1 + EnabledActiveParticlesCount; // DarkMatter and Other ActiveParticle (ex: SmartStar)
    yt_par_type *par_type_list = new yt_par_type [params->num_par_types];

    par_type_list[0].par_type = "DarkMatter";
    par_type_list[0].num_attr = 3 + 3 + 1 + 1 + 1 + NumberOfParticleAttributes;

    // the attributes name should be alive within the entire libyt in situ analysis,
    // because libyt does not make a copy of the names.
    std::vector<std::vector<std::string>> active_particles_attributes;
    std::vector<std::vector<hid_t>> active_particles_attributes_hdf5type;

    for (int i = 0; i < EnabledActiveParticlesCount; i++) {
        ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
        active_particles_attributes.emplace_back(ActiveParticleTypeToEvaluate->GetParticleAttributeNames());
        active_particles_attributes_hdf5type.emplace_back(ActiveParticleTypeToEvaluate->GetParticleAttributesHDF5DataType());

        par_type_list[1 + i].par_type = ActiveParticleTypeToEvaluate->particle_name.c_str();
        par_type_list[1 + i].num_attr = active_particles_attributes[i].size();
    }

    params->par_type_list = par_type_list;

    if (yt_set_Parameters(params) != YT_SUCCESS){
        fprintf(stderr, "Error in libyt API yt_set_Parameters\n");
        return FAIL;
    }

    delete [] par_type_list;

    yt_particle *particle_list;
    yt_get_ParticlesPtr(&particle_list);

    // TODO: make sure enzo's particle is always DarkMatter
    particle_list[0].par_type = "DarkMatter";

    // We have the attributes: 3 positions, 3 velocities, "mass", ID and Type
    // and extras.
    particle_list[0].num_attr = 3 + 3 + 1 + 1 + 1 + NumberOfParticleAttributes;
    const char *attr_name[] = {"particle_position_x",
                               "particle_position_y",
                               "particle_position_z",
                               "particle_velocity_x",
                               "particle_velocity_y",
                               "particle_velocity_z",
                               "particle_mass",
                               "particle_index",
                               "particle_type",
#ifdef WINDS
                               "creation_time",
                               "dynamical_time",
                               "metallicity_fraction",
                               "particle_jet_x", 
                               "particle_jet_y",
                               "particle_jet_z",
                               "typeia_fraction"
#else
                               "creation_time",
                               "dynamical_time",
                               "metallicity_fraction",
                               "typeia_fraction"
#endif
                              };
    for (int v = 0; v < particle_list[0].num_attr; v++) {
        particle_list[0].attr_list[v].attr_name = attr_name[v];
        particle_list[0].attr_list[v].attr_dtype = (v < 3) ? EYT_PFLOAT : EYT_BFLOAT;
    }
    // Now go back and reset it for index and type
    particle_list[0].attr_list[7].attr_dtype = EYT_PINT;
    particle_list[0].attr_list[8].attr_dtype = EYT_INT;
    particle_list[0].coor_x = attr_name[0];
    particle_list[0].coor_y = attr_name[1];
    particle_list[0].coor_z = attr_name[2];

    // we need to map hdf5type to yt datatype,
    // since we want to avoid storing yt data type in ParticleAttributeHandler class and causing
    for (int i = 0; i < EnabledActiveParticlesCount; i++) {
        for (int v = 0; v < active_particles_attributes[i].size(); v++) {
            particle_list[1 + i].attr_list[v].attr_name = active_particles_attributes[i][v].c_str();
            particle_list[1 + i].attr_list[v].attr_dtype = MapHDF5TypeToYTType(active_particles_attributes_hdf5type[i][v]);
        }

        particle_list[1 + i].coor_x = "particle_position_x";
        particle_list[1 + i].coor_y = "particle_position_y";
        particle_list[1 + i].coor_z = "particle_position_z";
    }

    /* Set code-specific parameter */
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
        field_list[libyt_field_i].field_dtype = EYT_BFLOAT;
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
     *  - Dust_Temperature (ignore for now)
     *  - Cooling_Time
     *
     *  Each of these is predicated on the global parameter associated with
     *  them.
     *
     * */

    field_list[libyt_field_i].field_name = "Temperature";
    field_list[libyt_field_i].field_type = "cell-centered";
    field_list[libyt_field_i].field_dtype = EYT_BFLOAT;
    for (j = 0; j < 6; j++) {
        field_list[libyt_field_i].field_ghost_cell[j] = NumberOfGhostZones;
    }
    libyt_field_i = libyt_field_i + 1;

    field_list[libyt_field_i].field_name = "Cooling_Time";
    field_list[libyt_field_i].field_type = "cell-centered";
    field_list[libyt_field_i].field_unit = "code_time";
    field_list[libyt_field_i].field_dtype = EYT_BFLOAT;
    for (j = 0; j < 6; j++) {
        field_list[libyt_field_i].field_ghost_cell[j] = NumberOfGhostZones;
    }

    /* We now have to do everything we do in CallPython.C, which amounts to
     *
     *  - ExposeGridHierarchy (not necessary anymore)
     *  - ExposeDataHierarchy (a recursive call)
     *
     * */

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

    /* Run yt_run_Function and yt_run_FunctionArguments */
	// if (yt_run_Function("yt_inline") != YT_SUCCESS) {
	// 	   fprintf(stderr, "Error while running yt_run_Function and call yt_inline\n");
	// 	   return FAIL;
	// }

    // if (yt_run_FunctionArguments("yt_inline_args", 1, "\'density\'") != YT_SUCCESS) {
 	//     fprintf(stderr, "Error while running yt_run_FunctionArguments and call yt_inline_args\n");
	//     return FAIL;
    // }

#ifdef USE_LIBYT_INTERACTIVE
    /* Call interactive Python prompt. */
    if (yt_run_InteractiveMode("LIBYT_STOP") != YT_SUCCESS) {
        fprintf(stderr, "Error in libyt API yt_run_InteractiveMode\n");
        fprintf(stderr, "One reason might be compiling libyt without -DINTERACTIVE_MODE=ON, "
                        "which does not support yt_run_InteractiveMode.\n");
    }
#endif

#ifdef USE_LIBYT_RELOAD
    /* Reloading script */
    if (yt_run_ReloadScript("LIBYT_STOP", "RELOAD", "reload.py") != YT_SUCCESS) {
        fprintf(stderr, "Error in libyt API yt_run_ReloadScript\n");
        fprintf(stderr, "One reason might be compiling libyt without -DINTERACTIVE_MODE=ON, "
                        "which does not support yt_run_ReloadScript.\n");
    }
#endif

#ifdef USE_LIBYT_JUPYTER
    /* Launch libyt Jupyter kernel */
    if (yt_run_JupyterKernel("LIBYT_STOP", false) != YT_SUCCESS) {
         fprintf(stderr, "Error in libyt API yt_run_JupyterKernel\n");
         fprintf(stderr, "One reason might be compiling libyt without -DJUPYTER_KERNEL=ON, "
                         "which does not support yt_run_JupyterKernel.\n");
    }
#endif

    /* Free resources allocated for libyt. */
    if (yt_free() != YT_SUCCESS) {
        fprintf(stderr, "Error in libyt API yt_free\n");
        return FAIL;
    }

    for (std::size_t i = 0; i < libyt_generated_data.size(); i++) {
        delete [] libyt_generated_data[i];
    }
    libyt_generated_data.clear();

    CommunicationBarrier();
    return SUCCESS;
#endif
}
