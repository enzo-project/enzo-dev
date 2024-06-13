#ifdef USE_LIBYT
/***********************************************************************
/
/  GRID CLASS (GIVE OUR INFO TO LIBYT)
/
/  written by: Matthew Turk
/  date:       April, 2023
/
/  PURPOSE:
/
/  RETURNS:
/    ENZO_SUCCESS or FAIL
/
************************************************************************/

#include <string>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "libyt.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

int GetUnits(float *DensityUnits, float *LengthUnits, float *TemperatureUnits, float *TimeUnits, float *VelocityUnits, FLOAT Time);

/* Note that previously we allowed WriteTime to be supplied and we do not now */
void grid::ConvertToLibyt(int LocalGridID, int GlobalGridID, int ParentID, int level, yt_grid &GridInfo)
{

    this->DebugCheck("Converting to libyt arrays");

    /* Declarations */

    /* libyt only wants us to set any array values if the grid is local.  Note
     * that this is opposite what we do in other places, where operations on
     * the hierarchy (that don't involve fields) are often replicated on all
     * processes. 
     *
     * *But* there may be times when we actually do want to do this.  So,
     * we don't do the check in here -- instead, we rely on it being checked
     * externally.
     *
     * In here, we'll happily do whatever we're told.
     *
     * */
  
    for (int i = 0; i < MAX_DIMENSION; i++) {
        GridInfo.grid_dimensions[i] = (this->GridEndIndex[i]) - (this->GridStartIndex[i]) + 1; // this is active dimension
        GridInfo.left_edge[i] = this->GridLeftEdge[i];
        GridInfo.right_edge[i] = this->GridRightEdge[i];
    }
    GridInfo.id = GlobalGridID;
    GridInfo.parent_id = ParentID;
    GridInfo.level = level;

    for (int field = 0; field < NumberOfBaryonFields; field++) {
        /* These are pointers, and *not* copies. Ownership is retained here. *
         * Note that we are only setting with the *current* BaryonField values,
         * and ignoring OldBaryonField.
         *
         * The other thing to note here is that we need to be extra careful
         * about fields being in different order in different grids.  So what we
         * do is use our global lookup table to map the field_type enum to the
         * libyt field entry.
         * 
         * This means we do not need to do any species field identification.
         *
         * (This may be far too cautious in almost all cases.  I mean,
         * grid::WriteGrid assumes that DataLabel corresponds to the field
         * index, and that is also what we assume internally for libyt
         * initialization.)
         *
         * */
        int libyt_field = libyt_field_lookup[field];
        if (libyt_field == -1) continue; // TODO: will this ever be -1?
        GridInfo.field_data[libyt_field].data_ptr = BaryonField[field];
    }

    long grid_size = GridDimension[0] * GridDimension[1] * GridDimension[2];
    float *temp_field = new float[grid_size];
    float *cooling_time_field = new float[grid_size];

    float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time);

    if (this->ComputeTemperatureField(temp_field) == SUCCESS) {
        libyt_generated_data.push_back(temp_field);
        int temp_field_i = ((yt_param_yt*)param_yt)->num_fields - 2;
        GridInfo.field_data[temp_field_i].data_ptr = temp_field;
    }

    if (this->ComputeCoolingTime(cooling_time_field) == SUCCESS) {
        for (long i = 0; i < grid_size; i++) { cooling_time_field[i] = fabs(cooling_time_field[i]) * TimeUnits; }
        libyt_generated_data.push_back(cooling_time_field);
        int cooling_time_field_i = ((yt_param_yt*)param_yt)->num_fields - 1;
        GridInfo.field_data[cooling_time_field_i].data_ptr = cooling_time_field;
    }

    /* par_count_list can take multiple particle types */
    if (this->ReturnNumberOfParticles() > 0) {
        GridInfo.par_count_list[0] = this->ReturnNumberOfParticles();
        for (int field = 0; field < 3; field++) {
            /* Set particle positions and velocities */
            GridInfo.particle_data[0][field].data_ptr = this->ParticlePosition[field];
            GridInfo.particle_data[0][field + 3].data_ptr = this->ParticleVelocity[field];
        }
        GridInfo.particle_data[0][6].data_ptr = this->ParticleMass;
        GridInfo.particle_data[0][7].data_ptr = this->ParticleNumber;
        GridInfo.particle_data[0][8].data_ptr = this->ParticleType;
        for (int field = 0; field < NumberOfParticleAttributes; field++) {
            GridInfo.particle_data[0][9 + field].data_ptr = this->ParticleAttribute[field];
        }
    }
    else {
        GridInfo.par_count_list[0] = 0;
    }

    /* fill in active particle count */
    std::vector<int> active_particle_count(EnabledActiveParticlesCount, 0);
    for (int i = 0; i < NumberOfActiveParticles; i++) {
        active_particle_count[(this->ActiveParticles)[i]->GetEnabledParticleID()]++;
    }
    for (int i = 0; i < EnabledActiveParticlesCount; i++) {
        GridInfo.par_count_list[1 + i] = active_particle_count[i];
    }

    /* fill in buffer only if there is active particle, which is sum of active particle count > 0
     * ignore AccretionRate/AccretionRateTime/Accreted_angmom for now since they are not one dimensional */
    if (NumberOfActiveParticles > 0) {
        for (int i = 0; i < EnabledActiveParticlesCount; i++) {
            ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];

            // the returned buffer ignores multi-dimensional array for now, it set them as nullptr.
            std::vector<void*> par_attr_buffer = ActiveParticleTypeToEvaluate->GetParticleAttributes(
                    this->ActiveParticles, i, NumberOfActiveParticles, active_particle_count[i],
                    ActiveParticleTypeToEvaluate->particle_name
                    );

            for (int a = 0; a < par_attr_buffer.size(); a++) {
                GridInfo.particle_data[1 + i][a].data_ptr = par_attr_buffer[a];
            }

            // free these pre-allocated buffer after libyt in situ analysis is done.
            libyt_generated_data.insert(libyt_generated_data.end(), par_attr_buffer.begin(), par_attr_buffer.end());
        }
    }
}
#endif
