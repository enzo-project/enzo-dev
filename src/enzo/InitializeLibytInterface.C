#ifdef USE_LIBYT
/***********************************************************************
/
/  INITIALIZE LIBYT INTERFACE
/
/  written by: Matthew Turk
/  date:       April, 2023
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "libyt.h"
#include "libyt_interactive_mode.h"

// We do this before any Enzo includes

#include <stdlib.h>
#include <stdio.h>

extern void* param_yt;
extern void* param_libyt;

int InitializeLibytByItself(long long argc, char *argv[])
{

    param_libyt = (void*) malloc(sizeof(yt_param_libyt));
    param_yt = (void*) malloc(sizeof(yt_param_yt));
    yt_param_libyt *params = (yt_param_libyt*) param_libyt;
    params->verbose = YT_VERBOSE_INFO;
    params->script = "inline";
    params->check_data = false;
    yt_initialize(argc, argv, params);
    fprintf(stderr, "Finished calling initialize!\n");
    return 0;
}

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

int InitializeLibytInterface(int argc, char *argv[])
{

    InitializeLibytByItself(argc, argv);
    char tempname[256];
    int i, j;

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
        if (!DataLabel[i]) continue;
        /* This tells us that BaryonFields[i] maps to
         * libyt_field[libyt_field_lookup[i]]
         *
         * */

        field_list[libyt_field_i].field_name = DataLabel[i];
        field_list[libyt_field_i].field_type = "cell-centered";
        field_list[libyt_field_i].field_dtype = ( sizeof(float) == 32 ) ? YT_FLOAT : YT_DOUBLE;
        for (j = 0; j < 6; j++) {
            /* This next bit is tricky.  The issue we have is that fields, in general,
             * will likely have the same number of ghost cells, regardless of
             * which grid they are on.  *But*, the way Enzo always does this is
             * to use the start and end index values of the grids to compute
             * the active dimensions.
             *
             * However, we don't have that information here -- we don't have
             * any grids, and we're just setting up the field info.
             *
             * */
            field_list[libyt_field_i].field_ghost_cell[j] = 3;
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

    

/* We call this parameter setting function here, but we *also* call it every
 * time, so that any updated parameters are caught. This is just to set the
 * stage. */

#include "InitializeLibytInterface_finderfunctions.inc"

  if (yt_run_InteractiveMode("LIBYT_STOP") != YT_SUCCESS) {
      return 1;
  }
  return 0;
}

int FinalizeLibytInterface()
{
  yt_finalize();
  if(debug)fprintf(stdout, "Completed Python interpreter finalization.\n");
  free(param_libyt);
  return SUCCESS;
}

void ExportParameterFileToLibyt(TopGridData *MetaData, FLOAT CurrentTime, FLOAT OldTime, 
			 float dtFixed)
{
  /* We need: */

  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
    VelocityUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, CurrentTime);

  if (ComovingCoordinates) {
    FLOAT a, dadt, FinalRedshift, CurrentRedshift;
    CosmologyComputeExpansionFactor(MetaData->StopTime, &a, &dadt);

    FinalRedshift = (1 + InitialRedshift)/a - 1;

    /* Compute the current redshift (for information only). */

    CosmologyComputeExpansionFactor(CurrentTime, &a, &dadt);
    CurrentRedshift = (1 + InitialRedshift)/a - 1;

    yt_set_UserParameterDouble("CosmologyCurrentRedshift", 1, &CurrentRedshift);
    yt_set_UserParameterDouble("CosmologyComovingBoxSize", 1, &ComovingBoxSize);
    yt_set_UserParameterDouble("CosmologyOmegaMatterNow", 1, &OmegaMatterNow);
    yt_set_UserParameterDouble("CosmologyOmegaLambdaNow", 1, &OmegaLambdaNow);
    yt_set_UserParameterDouble("CosmologyHubbleConstantNow", 1, &HubbleConstantNow);
    yt_set_UserParameterDouble("CosmologyInitialRedshift", 1, &InitialRedshift);
  }

  yt_set_UserParameterDouble("DensityUnits", 1, &DensityUnits);
  yt_set_UserParameterDouble("LengthUnits", 1, &LengthUnits);
  yt_set_UserParameterDouble("TemperatureUnits", 1, &TemperatureUnits);
  yt_set_UserParameterDouble("TimeUnits", 1, &TimeUnits);
  yt_set_UserParameterLongLong("HydroMethod", 1, &HydroMethod);
  yt_set_UserParameterLongLong("DualEnergyFormalism", 1, &DualEnergyFormalism);
  yt_set_UserParameterDouble("InitialTime", 1, &CurrentTime);
  yt_set_UserParameterDouble("StopTime", 1, &MetaData->StopTime);
  yt_set_UserParameterDouble("OldTime", 1, &OldTime);
  yt_set_UserParameterDouble("dtFixed", 1, &dtFixed);
  yt_set_UserParameterLongLong("ComovingCoordinates", 1, &ComovingCoordinates);

  /* Some of the things we do in InitializePythonInterface we do in the finder functions here:
   *
   *  - HydroMethod
   *  - RefineBy
   *  - DomainLeftEdge
   *  - DomainRightEdge
   *
   * Additionally, we do not need to specify the various conversion factors on
   * a field-by-field basis, as we provide them and supply the fields in code
   * units.
   */

  yt_set_UserParameterLongLong("TopGridRank", 1, &MetaData->TopGridRank);
  yt_set_UserParameterLongLong("NumberOfLibytCalls", 1, &NumberOfLibytCalls);
  yt_set_UserParameterLongLong("NumberOfLibytSubcycleCalls", 1, &NumberOfLibytSubcycleCalls);
  yt_set_UserParameterLongLong("NumberOfLibytTopGridCalls", 1, &NumberOfLibytTopGridCalls);
  yt_set_UserParameterDouble("CurrentMaximumDensity", 1, &CurrentMaximumDensity);
  yt_set_UserParameterDouble("AngularVelocity", 1, &AngularVelocity);
  yt_set_UserParameterDouble("VelocityGradient", 1, &VelocityGradient);
  yt_set_UserParameterDouble("dtDataDump", 1, &MetaData->dtDataDump);
  yt_set_UserParameterDouble("TimeLastDataDump", 1, &MetaData->TimeLastDataDump);
  yt_set_UserParameterLongLong("WroteData", 1, &MetaData->WroteData);

  yt_set_UserParameterLongLong("TopGridDimensions", MAX_DIMENSION, MetaData->TopGridDims);
  yt_set_UserParameterLongLong("LeftFaceBoundaryCondition", MAX_DIMENSION, MetaData->LeftFaceBoundaryCondition);

  return;
}


#endif
