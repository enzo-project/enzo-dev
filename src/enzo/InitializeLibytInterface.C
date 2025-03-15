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

// We do this before any Enzo includes

#include <stdlib.h>
#include <stdio.h>

#ifdef CONFIG_BFLOAT_4
#define YT_SET_USERPARAM_NONINTEGER(X, Y, Z) yt_set_UserParameterFloat(X, Y, Z)
#else
#define YT_SET_USERPARAM_NONINTEGER(X, Y, Z) yt_set_UserParameterDouble(X, Y, Z)
#endif

#ifdef SMALL_INTS
#define YT_SET_USERPARAM_INTEGER(X, Y, Z) yt_set_UserParameterInt(X, Y, Z)
#else
#define YT_SET_USERPARAM_INTEGER(X, Y, Z) yt_set_UserParameterLongLong(X, Y, Z)
#endif


extern char libyt_script_name[512];
extern void* param_yt;
extern void* param_libyt;

int InitializeLibytByItself(long long argc, char *argv[])
{

    param_libyt = (void*) malloc(sizeof(yt_param_libyt));
    param_yt = (void*) malloc(sizeof(yt_param_yt));
    yt_param_libyt *params = (yt_param_libyt*) param_libyt;
    params->verbose = YT_VERBOSE_INFO;
    params->script = libyt_script_name;  // inline.py will be loaded
    params->counter = 0;                 // Fig basename will start from count 0
    params->check_data = false;          // do not check hierarchy
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

    YT_SET_USERPARAM_NONINTEGER("CosmologyCurrentRedshift", 1, &CurrentRedshift);
    YT_SET_USERPARAM_NONINTEGER("CosmologyComovingBoxSize", 1, &ComovingBoxSize);
    YT_SET_USERPARAM_NONINTEGER("CosmologyOmegaMatterNow", 1, &OmegaMatterNow);
    YT_SET_USERPARAM_NONINTEGER("CosmologyOmegaLambdaNow", 1, &OmegaLambdaNow);
    YT_SET_USERPARAM_NONINTEGER("CosmologyHubbleConstantNow", 1, &HubbleConstantNow);
    YT_SET_USERPARAM_NONINTEGER("CosmologyInitialRedshift", 1, &InitialRedshift);
  }

  YT_SET_USERPARAM_NONINTEGER("DensityUnits", 1, &DensityUnits);
  YT_SET_USERPARAM_NONINTEGER("LengthUnits", 1, &LengthUnits);
  YT_SET_USERPARAM_NONINTEGER("TemperatureUnits", 1, &TemperatureUnits);
  YT_SET_USERPARAM_NONINTEGER("TimeUnits", 1, &TimeUnits);
  YT_SET_USERPARAM_INTEGER("HydroMethod", 1, &HydroMethod);
  YT_SET_USERPARAM_INTEGER("DualEnergyFormalism", 1, &DualEnergyFormalism);
  YT_SET_USERPARAM_NONINTEGER("InitialTime", 1, &CurrentTime);
  YT_SET_USERPARAM_NONINTEGER("StopTime", 1, &MetaData->StopTime);
  YT_SET_USERPARAM_NONINTEGER("OldTime", 1, &OldTime);
  YT_SET_USERPARAM_NONINTEGER("dtFixed", 1, &dtFixed);
  YT_SET_USERPARAM_INTEGER("ComovingCoordinates", 1, &ComovingCoordinates);

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

  YT_SET_USERPARAM_INTEGER("TopGridRank", 1, &MetaData->TopGridRank);
  YT_SET_USERPARAM_INTEGER("NumberOfLibytCalls", 1, &NumberOfLibytCalls);
  YT_SET_USERPARAM_INTEGER("NumberOfLibytSubcycleCalls", 1, &NumberOfLibytSubcycleCalls);
  YT_SET_USERPARAM_INTEGER("NumberOfLibytTopGridCalls", 1, &NumberOfLibytTopGridCalls);
  YT_SET_USERPARAM_NONINTEGER("CurrentMaximumDensity", 1, &CurrentMaximumDensity);
  YT_SET_USERPARAM_NONINTEGER("AngularVelocity", 1, &AngularVelocity);
  YT_SET_USERPARAM_NONINTEGER("VelocityGradient", 1, &VelocityGradient);
  YT_SET_USERPARAM_NONINTEGER("dtDataDump", 1, &MetaData->dtDataDump);
  YT_SET_USERPARAM_NONINTEGER("TimeLastDataDump", 1, &MetaData->TimeLastDataDump);
  YT_SET_USERPARAM_INTEGER("WroteData", 1, &MetaData->WroteData);

  YT_SET_USERPARAM_INTEGER("TopGridDimensions", MAX_DIMENSION, MetaData->TopGridDims);
  YT_SET_USERPARAM_INTEGER("LeftFaceBoundaryCondition", MAX_DIMENSION, MetaData->LeftFaceBoundaryCondition);

  return;
}

#undef YT_SET_USERPARAM_NONINTEGER
#undef YT_SET_USERPARAM_INTEGER

#endif
