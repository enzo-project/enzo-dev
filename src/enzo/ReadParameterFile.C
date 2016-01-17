/***********************************************************************
/
/  READ A PARAMETER FILE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       June 25, 2006
/  modified2:  Robert Harkness
/  date:       February 29th, 2008
/  modified3:  Robert Harkness
/  date:       May, 2008
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine reads the parameter file in the argument and sets parameters
//   based on it.

#include "preincludes.h" 
#include <stdlib.h>
#include <unistd.h>
#include <vector>

#ifdef CONFIG_USE_LIBCONFIG
#include <libconfig.h++>
#endif
 
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "hydro_rk/EOS.h" 
#include "CosmologyParameters.h"

/* This variable is declared here and only used in Grid_ReadGrid. */
 


/* function prototypes */

void my_exit(int status); 
int ReadListOfFloats(FILE *fptr, int N, float floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
int CosmologyReadParameters(FILE *fptr, FLOAT *StopTime, FLOAT *InitTime);
int ReadUnits(FILE *fptr);
int InitializeCosmicRayData();
int InitializeRateData(FLOAT Time);
int InitializeEquilibriumCoolData(FLOAT Time);
int InitializeGadgetEquilibriumCoolData(FLOAT Time);
int InitializeRadiationFieldData(FLOAT Time);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int ReadEvolveRefineFile(void);
 
int CheckShearingBoundaryConsistency(TopGridData &MetaData); 
void get_uuid(char *buffer);

int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt)
{
#ifndef CONFIG_USE_LIBCONFIG
  /* declarations */

  
  char line[MAX_LINE_LENGTH];
  int i, dim, ret, int_dummy;
  float TempFloat, float_dummy;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  int comment_count = 0;
 
  /* read until out of lines */

  rewind(fptr);
  while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL) 
      && (comment_count < 2)) {

    ret = 0;
 
    /* read MetaData parameters */
 
    ret += sscanf(line, "InitialCycleNumber = %"ISYM, &MetaData.CycleNumber);
    ret += sscanf(line, "InitialTime        = %"PSYM, &MetaData.Time);
    ret += sscanf(line, "InitialCPUTime     = %lf", &MetaData.CPUTime);
    ret += sscanf(line, "Initialdt          = %"FSYM, Initialdt);
 
    ret += sscanf(line, "CheckpointRestart = %"ISYM, &CheckpointRestart);
    ret += sscanf(line, "StopTime    = %"PSYM, &MetaData.StopTime);
    ret += sscanf(line, "StopCycle   = %"ISYM, &MetaData.StopCycle);
    ret += sscanf(line, "StopSteps   = %"ISYM, &MetaData.StopSteps);
    ret += sscanf(line, "StopCPUTime = %"FSYM, &MetaData.StopCPUTime);
    ret += sscanf(line, "ResubmitOn  = %"ISYM, &MetaData.ResubmitOn);
    if (sscanf(line, "ResubmitCommand = %s", dummy) == 1) 
      MetaData.ResubmitCommand = dummy;

    ret += sscanf(line, "MaximumTopGridTimeStep = %"FSYM,
		  &MetaData.MaximumTopGridTimeStep);

    ret += sscanf(line, "TimeLastRestartDump = %"FSYM,
		  &MetaData.TimeLastRestartDump);
    ret += sscanf(line, "dtRestartDump       = %"FSYM, &MetaData.dtRestartDump);
    ret += sscanf(line, "TimeLastDataDump    = %"PSYM,
		  &MetaData.TimeLastDataDump);
    ret += sscanf(line, "dtDataDump          = %"PSYM, &MetaData.dtDataDump);
    ret += sscanf(line, "TimeLastHistoryDump = %"PSYM,
		  &MetaData.TimeLastHistoryDump);
    ret += sscanf(line, "dtHistoryDump       = %"PSYM, &MetaData.dtHistoryDump);
 
    ret += sscanf(line, "TracerParticleOn  = %"ISYM, &TracerParticleOn);
    ret += sscanf(line, "TracerParticleOutputVelocity  = %"ISYM, &TracerParticleOutputVelocity);
    ret += sscanf(line, "WriteGhostZones = %"ISYM, &WriteGhostZones);
    ret += sscanf(line, "ReadGhostZones = %"ISYM, &ReadGhostZones);
    ret += sscanf(line, "OutputParticleTypeGrouping = %"ISYM,
                        &OutputParticleTypeGrouping);
    ret += sscanf(line, "TimeLastTracerParticleDump = %"PSYM,
                  &MetaData.TimeLastTracerParticleDump);
    ret += sscanf(line, "dtTracerParticleDump       = %"PSYM,
                  &MetaData.dtTracerParticleDump);
    ret += sscanf(line, "TimeLastInterpolatedDataDump    = %"PSYM,
		  &MetaData.TimeLastInterpolatedDataDump);
    ret += sscanf(line, "dtInterpolatedDataDump          = %"PSYM, 
		  &MetaData.dtInterpolatedDataDump);
 
    ret += sscanf(line, "NewMovieLeftEdge  = %"PSYM" %"PSYM" %"PSYM, 
		  MetaData.NewMovieLeftEdge,
		  MetaData.NewMovieLeftEdge+1, 
		  MetaData.NewMovieLeftEdge+2);
    ret += sscanf(line, "NewMovieRightEdge = %"PSYM" %"PSYM" %"PSYM, 
		  MetaData.NewMovieRightEdge, 
		  MetaData.NewMovieRightEdge+1,
		  MetaData.NewMovieRightEdge+2);

    ret += sscanf(line, "CycleLastRestartDump = %"ISYM,
		  &MetaData.CycleLastRestartDump);
    ret += sscanf(line, "CycleSkipRestartDump = %"ISYM,
		  &MetaData.CycleSkipRestartDump);
    ret += sscanf(line, "CycleLastDataDump    = %"ISYM,
		  &MetaData.CycleLastDataDump);
    ret += sscanf(line, "CycleSkipDataDump    = %"ISYM,
		  &MetaData.CycleSkipDataDump);
    ret += sscanf(line, "CycleLastHistoryDump = %"ISYM,
		  &MetaData.CycleLastHistoryDump);
    ret += sscanf(line, "CycleSkipHistoryDump = %"ISYM,
		  &MetaData.CycleSkipHistoryDump);
    ret += sscanf(line, "CycleSkipGlobalDataDump = %"ISYM, //AK
                  &MetaData.CycleSkipGlobalDataDump);
    ret += sscanf(line, "OutputFirstTimeAtLevel = %"ISYM,
		  &MetaData.OutputFirstTimeAtLevel);
    ret += sscanf(line, "StopFirstTimeAtLevel = %"ISYM,
		  &MetaData.StopFirstTimeAtLevel);
    ret += sscanf(line, "NumberOfOutputsBeforeExit = %"ISYM,
		  &MetaData.NumberOfOutputsBeforeExit);

    /* Maximum density directed output */
    ret += sscanf(line, "OutputOnDensity = %"ISYM,
           &OutputOnDensity);
    ret += sscanf(line, "StartDensityOutputs = %"FSYM,
           &StartDensityOutputs);
    ret += sscanf(line, "CurrentDensityOutput = %"FSYM,
           &CurrentDensityOutput);
    ret += sscanf(line, "IncrementDensityOutput = %"FSYM,
           &IncrementDensityOutput);
    ret += sscanf(line, "StopFirstTimeAtDensity = %"FSYM,
           &StopFirstTimeAtDensity);
    ret += sscanf(line, "StopFirstTimeAtMetalEnrichedDensity = %"FSYM,
           &StopFirstTimeAtMetalEnrichedDensity);
    ret += sscanf(line, "EnrichedMetalFraction = %"FSYM,
           &EnrichedMetalFraction);

    /* Subcycle directed output */
    ret += sscanf(line, "SubcycleSkipDataDump = %"ISYM, 
                  &MetaData.SubcycleSkipDataDump);
    ret += sscanf(line, "SubcycleLastDataDump = %"ISYM, 
                  &MetaData.SubcycleLastDataDump);
    ret += sscanf(line, "SubcycleNumber = %"ISYM, 
                  &MetaData.SubcycleNumber);

    ret += sscanf(line,"FileDirectedOutput = %"ISYM,
		  &FileDirectedOutput);

    ret += sscanf(line,"HierarchyFileInputFormat = %"ISYM,
		  &HierarchyFileInputFormat);    
    ret += sscanf(line,"HierarchyFileOutputFormat = %"ISYM,
		  &HierarchyFileOutputFormat);    

    ret += sscanf(line, "RestartDumpNumber = %"ISYM, &MetaData.RestartDumpNumber);
    ret += sscanf(line, "DataDumpNumber    = %"ISYM, &MetaData.DataDumpNumber);
    ret += sscanf(line, "HistoryDumpNumber = %"ISYM, &MetaData.HistoryDumpNumber);
    ret += sscanf(line, "TracerParticleDumpNumber = %"ISYM, &MetaData.TracerParticleDumpNumber);
 
    if (sscanf(line, "RestartDumpName      = %s", dummy) == 1)
      MetaData.RestartDumpName = dummy;
    if (sscanf(line, "DataDumpName         = %s", dummy) == 1)
      MetaData.DataDumpName = dummy;
    if (sscanf(line, "HistoryDumpName      = %s", dummy) == 1)
      MetaData.HistoryDumpName = dummy;
    if (sscanf(line, "TracerParticleDumpName = %s", dummy) == 1)
      MetaData.TracerParticleDumpName = dummy;
    if (sscanf(line, "RedshiftDumpName     = %s", dummy) == 1)
      MetaData.RedshiftDumpName = dummy;
 
    if (sscanf(line, "RestartDumpDir      = %s", dummy) == 1)
      MetaData.RestartDumpDir = dummy;
    if (sscanf(line, "DataDumpDir         = %s", dummy) == 1)
      MetaData.DataDumpDir = dummy;
    if (sscanf(line, "HistoryDumpDir      = %s", dummy) == 1)
      MetaData.HistoryDumpDir = dummy;
    if (sscanf(line, "TracerParticleDumpDir = %s", dummy) == 1)
      MetaData.TracerParticleDumpDir = dummy;
    if (sscanf(line, "RedshiftDumpDir     = %s", dummy) == 1)
      MetaData.RedshiftDumpDir = dummy;
 
    if (sscanf(line, "LocalDir            = %s", dummy) == 1)
      MetaData.LocalDir = dummy;
    if (sscanf(line, "GlobalDir           = %s", dummy) == 1)
      MetaData.GlobalDir = dummy;
 
    if (sscanf(line, "CubeDump[%"ISYM"] = %s", &dim, dummy) == 2) {
      ret++; CubeDumps[dim] = dummy;
      if (dim >= MAX_CUBE_DUMPS) {
        ENZO_VFAIL("CubeDump %"ISYM" > maximum allowed.\n", dim)
      }
    }

    ret += sscanf(line, "NumberOfGhostZones = %"ISYM, &NumberOfGhostZones);
    ret += sscanf(line, "LoadBalancing = %"ISYM, &LoadBalancing);
    ret += sscanf(line, "ResetLoadBalancing = %"ISYM, &ResetLoadBalancing);
    ret += sscanf(line, "LoadBalancingCycleSkip = %"ISYM, &LoadBalancingCycleSkip);
    ret += sscanf(line, "LoadBalancingMinLevel = %"ISYM, &LoadBalancingMinLevel);
    ret += sscanf(line, "LoadBalancingMaxLevel = %"ISYM, &LoadBalancingMaxLevel);
 
    ret += sscanf(line, "ConductionDynamicRebuildHierarchy = %"ISYM, 
                  &ConductionDynamicRebuildHierarchy);
    ret += sscanf(line, "ConductionDynamicRebuildMinLevel = %"ISYM, 
                  &ConductionDynamicRebuildMinLevel);
    if (sscanf(line, "RebuildHierarchyCycleSkip[%"ISYM"] =", &int_dummy) == 1) {
      if (int_dummy > MAX_DEPTH_OF_HIERARCHY) {
	ENZO_VFAIL("Cannot set RebuildHierarchyCycleSkip[%"ISYM"], max hierarchy depth = %"ISYM".\n", int_dummy, MAX_DEPTH_OF_HIERARCHY);
      }
      ret += sscanf(line, "RebuildHierarchyCycleSkip[%"ISYM"] = %"ISYM,
		    &int_dummy, &RebuildHierarchyCycleSkip[int_dummy]);
    }

    if (sscanf(line, "TimeActionType[%"ISYM"] = %"ISYM, &dim, &int_dummy) == 2) {
      ret++;
      if (dim >= MAX_TIME_ACTIONS-1) {
	ENZO_VFAIL("Time action %"ISYM" > maximum allowed.\n", dim)
      }
      TimeActionType[dim] = int_dummy;
    }
    if (sscanf(line, "TimeActionRedshift[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionRedshift[%"ISYM"] = %"PSYM, &dim,
		    TimeActionRedshift+dim);
    if (sscanf(line, "TimeActionTime[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionTime[%"ISYM"] = %"PSYM, &dim,
		    TimeActionTime+dim);
    if (sscanf(line, "TimeActionParameter[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line, "TimeActionParameter[%"ISYM"] = %"FSYM, &dim,
		    TimeActionParameter+dim);
 
    ret += sscanf(line, "StaticHierarchy = %"ISYM, &MetaData.StaticHierarchy);
 
    ret += sscanf(line, "TopGridRank       = %"ISYM, &MetaData.TopGridRank);
    ret += sscanf(line, "TopGridDimensions = %"ISYM" %"ISYM" %"ISYM, MetaData.TopGridDims,
		  MetaData.TopGridDims+1, MetaData.TopGridDims+2);
 
    ret += sscanf(line, "TopGridGravityBoundary = %"ISYM,
		  &MetaData.GravityBoundary);
 
#ifdef TRANSFER
    if (sscanf(line, "RadHydroParamfile = %s", dummy) == 1)
      MetaData.RadHydroParameterFname = dummy;
#endif
    ret += sscanf(line, "ImplicitProblem = %"ISYM, &ImplicitProblem);
    ret += sscanf(line, "RadiativeTransferFLD   = %"ISYM, &RadiativeTransferFLD);
#ifdef EMISSIVITY
    ret += sscanf(line, "StarMakerEmissivityField = %"ISYM, 
		  &StarMakerEmissivityField);
    ret += sscanf(line, "uv_param = %"FSYM, &uv_param);
#endif

    ret += sscanf(line, "ParticleBoundaryType   = %"ISYM,
		  &MetaData.ParticleBoundaryType);
    ret += sscanf(line, "NumberOfParticles      = %"PISYM,
		  &MetaData.NumberOfParticles);
 
    ret += sscanf(line, "CourantSafetyNumber    = %"FSYM,
		  &MetaData.CourantSafetyNumber);
    ret += sscanf(line, "PPMFlatteningParameter = %"ISYM,
		  &MetaData.PPMFlatteningParameter);
    ret += sscanf(line, "PPMDiffusionParameter  = %"ISYM,
		  &MetaData.PPMDiffusionParameter);
    ret += sscanf(line, "PPMSteepeningParameter = %"ISYM,
		  &MetaData.PPMSteepeningParameter);
 
    /* read global Parameters */
 
    ret += sscanf(line, "ProblemType            = %"ISYM, &ProblemType);
#ifdef NEW_PROBLEM_TYPES
    if (sscanf(line, "ProblemTypeName = %s", dummy) == 1) {
      ProblemTypeName = dummy;
      ProblemType = -978;
      ret = 1;
    }
#endif
    ret += sscanf(line, "HydroMethod            = %"ISYM, &HydroMethod);

    ret += sscanf(line, "huge_number            = %"FSYM, &huge_number);
    ret += sscanf(line, "tiny_number            = %"FSYM, &tiny_number);
    ret += sscanf(line, "Gamma                  = %"FSYM, &Gamma);
    ret += sscanf(line, "PressureFree           = %"ISYM, &PressureFree);
    ret += sscanf(line, "RefineBy               = %"ISYM, &RefineBy);
    ret += sscanf(line, "MaximumRefinementLevel = %"ISYM,
		  &MaximumRefinementLevel);
    ret += sscanf(line, "MaximumGravityRefinementLevel = %"ISYM,
		  &MaximumGravityRefinementLevel);
    ret += sscanf(line, "MaximumParticleRefinementLevel = %"ISYM,
		  &MaximumParticleRefinementLevel);
    ret += sscanf(line, "CellFlaggingMethod     = %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM,
	     CellFlaggingMethod+0, CellFlaggingMethod+1, CellFlaggingMethod+2,
	     CellFlaggingMethod+3, CellFlaggingMethod+4, CellFlaggingMethod+5,
	     CellFlaggingMethod+6);
    ret += sscanf(line, "FluxCorrection         = %"ISYM, &FluxCorrection);
    ret += sscanf(line, "UseCoolingTimestep     = %"ISYM, &UseCoolingTimestep);
    ret += sscanf(line, "CoolingTimestepSafetyFactor = %"FSYM, &CoolingTimestepSafetyFactor);
    ret += sscanf(line, "InterpolationMethod    = %"ISYM, &InterpolationMethod);
    ret += sscanf(line, "ConservativeInterpolation = %"ISYM,
		  &ConservativeInterpolation);
    ret += sscanf(line, "MinimumEfficiency      = %"FSYM, &MinimumEfficiency);
    ret += sscanf(line, "SubgridSizeAutoAdjust  = %"ISYM, &SubgridSizeAutoAdjust);
    ret += sscanf(line, "OptimalSubgridsPerProcessor = %"ISYM, 
		  &OptimalSubgridsPerProcessor);
    ret += sscanf(line, "MinimumSubgridEdge     = %"ISYM, &MinimumSubgridEdge);
    ret += sscanf(line, "MaximumSubgridSize     = %"ISYM, &MaximumSubgridSize);
    ret += sscanf(line, "CriticalGridRatio      = %"FSYM, &CriticalGridRatio);
    ret += sscanf(line, "NumberOfBufferZones    = %"ISYM, &NumberOfBufferZones);
    ret += sscanf(line, "FastSiblingLocatorEntireDomain = %"ISYM, &FastSiblingLocatorEntireDomain);
    ret += sscanf(line, "MustRefineRegionMinRefinementLevel = %"ISYM,
		  &MustRefineRegionMinRefinementLevel);
    ret += sscanf(line, "MetallicityRefinementMinLevel = %"ISYM,
		  &MetallicityRefinementMinLevel);
    ret += sscanf(line, "MetallicityRefinementMinMetallicity = %"FSYM, 
		  &MetallicityRefinementMinMetallicity);
    ret += sscanf(line, "MetallicityRefinementMinDensity = %"FSYM, 
		  &MetallicityRefinementMinDensity);

    ret += sscanf(line, "DomainLeftEdge        = %"PSYM" %"PSYM" %"PSYM, DomainLeftEdge,
		  DomainLeftEdge+1, DomainLeftEdge+2);
    ret += sscanf(line, "DomainRightEdge       = %"PSYM" %"PSYM" %"PSYM, DomainRightEdge,
		  DomainRightEdge+1, DomainRightEdge+2);
    ret += sscanf(line, "GridVelocity          = %"FSYM" %"FSYM" %"FSYM, GridVelocity,
		  GridVelocity+1, GridVelocity+2);
    ret += sscanf(line, "RefineRegionAutoAdjust = %"ISYM, &RefineRegionAutoAdjust);
    ret += sscanf(line, "RefineRegionLeftEdge  = %"PSYM" %"PSYM" %"PSYM,
		  RefineRegionLeftEdge, RefineRegionLeftEdge+1,
		  RefineRegionLeftEdge+2);
    ret += sscanf(line, "RefineRegionRightEdge = %"PSYM" %"PSYM" %"PSYM,
		  RefineRegionRightEdge, RefineRegionRightEdge+1,
		  RefineRegionRightEdge+2);
    ret += sscanf(line, "MustRefineRegionLeftEdge  = %"PSYM" %"PSYM" %"PSYM,
		  MustRefineRegionLeftEdge, MustRefineRegionLeftEdge+1,
		  MustRefineRegionLeftEdge+2);
    ret += sscanf(line, "MustRefineRegionRightEdge  = %"PSYM" %"PSYM" %"PSYM,
		  MustRefineRegionRightEdge, MustRefineRegionRightEdge+1,
		  MustRefineRegionRightEdge+2);
  
    /* Parameters for the MultiRefineRegion mechanics */
    
    ret += sscanf(line, "MultiRefineRegionMaximumOuterLevel  = %"ISYM, &MultiRefineRegionMaximumOuterLevel);
    ret += sscanf(line, "MultiRefineRegionMinimumOuterLevel  = %"ISYM, &MultiRefineRegionMinimumOuterLevel);
    if (sscanf(line, "MultiRefineRegionMaximumLevel[%"ISYM"] = %"ISYM, &dim, &int_dummy) == 2) 
      {
	if (dim > MAX_STATIC_REGIONS-1) 
	  ENZO_VFAIL("MultiRefineRegion number %"ISYM" (MAX_STATIC_REGIONS) > MAX allowed\n", dim);
	ret++;
	MultiRefineRegionMaximumLevel[dim] = int_dummy;
      }
    if (sscanf(line, "MultiRefineRegionGeometry[%"ISYM"] = %"ISYM, &dim, &int_dummy) == 2){
      ret++;
      MultiRefineRegionGeometry[dim] = int_dummy;
    }
    if (sscanf(line, "MultiRefineRegionMinimumLevel[%"ISYM"] = %"ISYM, &dim, &int_dummy) == 2){
      ret++;
      MultiRefineRegionMinimumLevel[dim] = int_dummy;
    }
    if (sscanf(line, "MultiRefineRegionRadius[%"ISYM"] = %"PSYM, &dim, &float_dummy) == 2){
      ret++;
      MultiRefineRegionRadius[dim] = float_dummy;
    }
    if (sscanf(line, "MultiRefineRegionWidth[%"ISYM"] = %"PSYM, &dim, &float_dummy) == 2){
      ret++;
      MultiRefineRegionWidth[dim] = float_dummy;
    }
    if (sscanf(line, "MultiRefineRegionStaggeredRefinement[%"ISYM"] = %"PSYM, &dim, &float_dummy) == 2){
      ret++;
      MultiRefineRegionStaggeredRefinement[dim] = float_dummy;
    }
    if (sscanf(line, "MultiRefineRegionCenter[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "MultiRefineRegionCenter[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, MultiRefineRegionCenter[dim],
		    MultiRefineRegionCenter[dim]+1,
		    MultiRefineRegionCenter[dim]+2);
    if (sscanf(line, "MultiRefineRegionOrientation[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "MultiRefineRegionOrientation[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, MultiRefineRegionOrientation[dim],
		    MultiRefineRegionOrientation[dim]+1,
		    MultiRefineRegionOrientation[dim]+2);
    if (sscanf(line, "MultiRefineRegionLeftEdge[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "MultiRefineRegionLeftEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, MultiRefineRegionLeftEdge[dim],
		    MultiRefineRegionLeftEdge[dim]+1,
		    MultiRefineRegionLeftEdge[dim]+2);
    if (sscanf(line, "MultiRefineRegionRightEdge[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "MultiRefineRegionRightEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, MultiRefineRegionRightEdge[dim],
		    MultiRefineRegionRightEdge[dim]+1,
		    MultiRefineRegionRightEdge[dim]+2);

    /* Read evolving RefineRegion */

    ret += sscanf(line, "RefineRegionTimeType = %"ISYM, &RefineRegionTimeType);
    if (sscanf(line, "RefineRegionFile = %s", dummy) == 1) {
      RefineRegionFile = dummy;
      ret++;
    }

    if (sscanf(line, "DatabaseLocation = %s", dummy) == 1) {
      DatabaseLocation = dummy;
      ret++;
    }

    if (sscanf(line, "DataLabel[%"ISYM"] = %s\n", &dim, dummy) == 2)
      DataLabel[dim] = dummy;
    if (sscanf(line, "DataUnits[%"ISYM"] = %s\n", &dim, dummy) == 2)
      DataUnits[dim] = dummy;
 
    ret += sscanf(line, "UniformGravity          = %"ISYM, &UniformGravity);
    ret += sscanf(line, "UniformGravityDirection = %"ISYM,
		  &UniformGravityDirection);
    ret += sscanf(line, "UniformGravityConstant  = %"FSYM,
		  &UniformGravityConstant);
 
    ret += sscanf(line, "PointSourceGravity         = %"ISYM,&PointSourceGravity);
    ret += sscanf(line, "PointSourceGravityPosition = %"PSYM" %"PSYM" %"PSYM,
		  PointSourceGravityPosition, PointSourceGravityPosition+1,
		  PointSourceGravityPosition+2);
    ret += sscanf(line, "PointSourceGravityConstant = %"FSYM,
		  &PointSourceGravityConstant);
    ret += sscanf(line, "PointSourceGravityCoreRadius = %"FSYM,
		  &PointSourceGravityCoreRadius);
 
    ret += sscanf(line, "DiskGravity                        = %"ISYM,&DiskGravity);
    ret += sscanf(line, "DiskGravityPosition                = %"PSYM" %"PSYM" %"PSYM,
      DiskGravityPosition, DiskGravityPosition+1, DiskGravityPosition+2);
    ret += sscanf(line, "DiskGravityAngularMomentum         = %"PSYM" %"PSYM" %"PSYM,
      DiskGravityAngularMomentum,DiskGravityAngularMomentum+1,
      DiskGravityAngularMomentum+2);
    ret += sscanf(line, "DiskGravityStellarDiskMass         = %"FSYM,&DiskGravityStellarDiskMass);
    ret += sscanf(line, "DiskGravityStellarDiskScaleHeightR = %"FSYM,&DiskGravityStellarDiskScaleHeightR);
    ret += sscanf(line, "DiskGravityStellarDiskScaleHeightz = %"FSYM,&DiskGravityStellarDiskScaleHeightz);
    ret += sscanf(line, "DiskGravityStellarBulgeMass        = %"FSYM,&DiskGravityStellarBulgeMass);
    ret += sscanf(line, "DiskGravityStellarBulgeR           = %"FSYM,&DiskGravityStellarBulgeR);
    ret += sscanf(line, "DiskGravityDarkMatterR             = %"FSYM,&DiskGravityDarkMatterR);
    ret += sscanf(line, "DiskGravityDarkMatterDensity       = %"FSYM,&DiskGravityDarkMatterDensity);

    ret += sscanf(line, "ExternalGravity         = %"ISYM,&ExternalGravity);
    ret += sscanf(line, "ExternalGravityConstant = %"FSYM, &ExternalGravityConstant);
    ret += sscanf(line, "ExternalGravityRadius   = %"FSYM,&ExternalGravityRadius);
    ret += sscanf(line, "ExternalGravityDensity  = %"FSYM,&ExternalGravityDensity);
    ret += sscanf(line, "ExternalGravityPosition = %"PSYM" %"PSYM" %"PSYM,
		  ExternalGravityPosition, ExternalGravityPosition+1,
		  ExternalGravityPosition+2);
    ret += sscanf(line, "ExternalGravityOrientation = %"FSYM" %"FSYM" %"FSYM, 
		  ExternalGravityOrientation, ExternalGravityOrientation+1, 
		  ExternalGravityOrientation+2);

    ret += sscanf(line, "SelfGravity           = %"ISYM, &SelfGravity);
    ret += sscanf(line, "SelfGravityGasOff     = %"ISYM, &SelfGravityGasOff);
    ret += sscanf(line, "AccretionKernal       = %"ISYM, &AccretionKernal);
    ret += sscanf(line, "GravitationalConstant = %"FSYM, &GravitationalConstant);
    ret += sscanf(line, "S2ParticleSize        = %"FSYM, &S2ParticleSize);
    ret += sscanf(line, "GravityResolution     = %"FSYM, &GravityResolution);
    ret += sscanf(line, "ComputePotential      = %"ISYM, &ComputePotential);
    ret += sscanf(line, "PotentialIterations   = %"ISYM, &PotentialIterations);
    ret += sscanf(line, "WritePotential        = %"ISYM, &WritePotential);
    ret += sscanf(line, "ParticleSubgridDepositMode  = %"ISYM, &ParticleSubgridDepositMode);
    ret += sscanf(line, "WriteAcceleration      = %"ISYM, &WriteAcceleration);
    ret += sscanf(line, "BaryonSelfGravityApproximation = %"ISYM,
		  &BaryonSelfGravityApproximation);
 
    ret += sscanf(line, "GreensFunctionMaxNumber   = %"ISYM,
		  &GreensFunctionMaxNumber);
    ret += sscanf(line, "GreensFunctionMaxSize     = %"ISYM,
		  &GreensFunctionMaxSize);
 
    ret += sscanf(line, "DualEnergyFormalism     = %"ISYM, &DualEnergyFormalism);
    ret += sscanf(line, "DualEnergyFormalismEta1 = %"FSYM,
		  &DualEnergyFormalismEta1);
    ret += sscanf(line, "DualEnergyFormalismEta2 = %"FSYM,
		  &DualEnergyFormalismEta2);
    ret += sscanf(line, "ParticleCourantSafetyNumber = %"FSYM,
		  &ParticleCourantSafetyNumber);
    ret += sscanf(line, "RootGridCourantSafetyNumber = %"FSYM,
		  &RootGridCourantSafetyNumber);
    ret += sscanf(line, "RandomForcing = %"ISYM, &RandomForcing); //AK
    ret += sscanf(line, "RandomForcingEdot = %"FSYM, &RandomForcingEdot); //AK
    ret += sscanf(line, "RandomForcingMachNumber = %"FSYM, //AK
                  &RandomForcingMachNumber);

    ret += sscanf(line, "DrivenFlowProfile = %"ISYM, &DrivenFlowProfile);
    ret += sscanf(line, "DrivenFlowWeight = %"FSYM, &DrivenFlowWeight);
    ret += sscanf(line, "DrivenFlowAlpha = %"ISYM" %"ISYM" %"ISYM,
                  DrivenFlowAlpha, DrivenFlowAlpha+1, DrivenFlowAlpha+2);
    ret += sscanf(line, "DrivenFlowSeed = %"ISYM, &DrivenFlowSeed);
    ret += sscanf(line, "DrivenFlowBandWidth = %"FSYM"%"FSYM"%"FSYM,
                  DrivenFlowBandWidth, DrivenFlowBandWidth+1, DrivenFlowBandWidth+2);
    ret += sscanf(line, "DrivenFlowVelocity = %"FSYM"%"FSYM"%"FSYM,
                  DrivenFlowVelocity, DrivenFlowVelocity+1, DrivenFlowVelocity+2);
    ret += sscanf(line, "DrivenFlowAutoCorrl = %"FSYM"%"FSYM"%"FSYM,
                     DrivenFlowAutoCorrl, DrivenFlowAutoCorrl+1, DrivenFlowAutoCorrl+2);

#ifdef USE_GRACKLE
    /* Grackle chemistry parameters */
    ret += sscanf(line, "use_grackle = %d", &grackle_data.use_grackle);
    ret += sscanf(line, "with_radiative_cooling = %d",
                  &grackle_data.with_radiative_cooling);
    if (sscanf(line, "grackle_data_file = %s", dummy) == 1) {
      grackle_data.grackle_data_file = dummy;
      ret++;
    }
    ret += sscanf(line, "UVbackground = %d", &grackle_data.UVbackground);
    ret += sscanf(line, "Compton_xray_heating = %d", 
                  &grackle_data.Compton_xray_heating);
    ret += sscanf(line, "LWbackground_intensity = %lf", 
                  &grackle_data.LWbackground_intensity);
    ret += sscanf(line, "LWbackground_sawtooth_suppression = %d",
                  &grackle_data.LWbackground_sawtooth_suppression);
    /********************************/
#endif
    ret += sscanf(line, "RadiativeCooling = %"ISYM, &RadiativeCooling);
    ret += sscanf(line, "RadiativeCoolingModel = %"ISYM, &RadiativeCoolingModel);
    ret += sscanf(line, "GadgetEquilibriumCooling = %"ISYM, &GadgetEquilibriumCooling);
    ret += sscanf(line, "MultiSpecies = %"ISYM, &MultiSpecies);
    ret += sscanf(line, "CIECooling = %"ISYM, &CIECooling);
    ret += sscanf(line, "H2OpticalDepthApproximation = %"ISYM, &H2OpticalDepthApproximation);
    ret += sscanf(line, "ThreeBodyRate = %"ISYM, &ThreeBodyRate);
    ret += sscanf(line, "H2FormationOnDust = %"ISYM, &H2FormationOnDust);
    if (sscanf(line, "CloudyCoolingGridFile = %s", dummy) == 1) {
      CloudyCoolingData.CloudyCoolingGridFile = dummy;
      ret++;
    }
    ret += sscanf(line, "IncludeCloudyHeating = %"ISYM, &CloudyCoolingData.IncludeCloudyHeating);
    ret += sscanf(line, "CMBTemperatureFloor = %"ISYM, &CloudyCoolingData.CMBTemperatureFloor);
    ret += sscanf(line, "CloudyElectronFractionFactor = %"FSYM,&CloudyCoolingData.CloudyElectronFractionFactor);
    ret += sscanf(line, "MetalCooling = %"ISYM"", &MetalCooling);
    if (sscanf(line, "MetalCoolingTable = %s", dummy) == 1) {
      MetalCoolingTable = dummy;
      ret++;
    }

    ret += sscanf(line, "CRModel = %"ISYM, &CRModel); 
    ret += sscanf(line, "CRDiffusion = %"ISYM, &CRDiffusion);
    ret += sscanf(line, "CRkappa = %"FSYM, &CRkappa);
    ret += sscanf(line, "CRCourantSafetyNumber = %"FSYM, &CRCourantSafetyNumber);
    ret += sscanf(line, "CRFeedback = %"FSYM, &CRFeedback);
    ret += sscanf(line, "CRdensFloor = %"FSYM, &CRdensFloor);
    ret += sscanf(line, "CRmaxSoundSpeed = %"FSYM, &CRmaxSoundSpeed);
    ret += sscanf(line, "CRgamma = %"FSYM, &CRgamma);
    ret += sscanf(line, "CosmologySimulationUniformCR = %"FSYM, &CosmologySimulationUniformCR); // FIXME

    ret += sscanf(line, "ShockMethod = %"ISYM, &ShockMethod);
    ret += sscanf(line, "ShockTemperatureFloor = %"FSYM, &ShockTemperatureFloor);
    ret += sscanf(line, "StorePreShockFields = %"ISYM, &StorePreShockFields);
    ret += sscanf(line, "FindShocksOnlyOnOutput = %"ISYM, &FindShocksOnlyOnOutput);

    ret += sscanf(line, "RadiationFieldType = %"ISYM, &RadiationFieldType);
    ret += sscanf(line, "RadiationFieldRedshift = %"FSYM, &RadiationFieldRedshift);
    ret += sscanf(line, "TabulatedLWBackground = %"ISYM, &TabulatedLWBackground);
    ret += sscanf(line, "AdjustUVBackground = %"ISYM, &AdjustUVBackground);
    ret += sscanf(line, "AdjustUVBackgroundHighRedshift = %"ISYM, &AdjustUVBackgroundHighRedshift);
    ret += sscanf(line, "SetUVBAmplitude = %"FSYM, &SetUVBAmplitude);
    ret += sscanf(line, "SetHeIIHeatingScale = %"FSYM, &SetHeIIHeatingScale);
    ret += sscanf(line, "RadiationFieldLevelRecompute = %"ISYM, &RadiationFieldLevelRecompute);    
    ret += sscanf(line, "RadiationShield = %"ISYM, &RadiationData.RadiationShield);
    ret += sscanf(line, "RadiationSpectrumNormalization = %"FSYM, &CoolData.f3);
    ret += sscanf(line, "RadiationSpectrumSlope = %"FSYM, &CoolData.alpha0);
    ret += sscanf(line, "CoolDataf0to3 = %"FSYM, &CoolData.f0to3);
    ret += sscanf(line, "RadiationRedshiftOn = %"FSYM, &CoolData.RadiationRedshiftOn);
    ret += sscanf(line, "RadiationRedshiftOff = %"FSYM, &CoolData.RadiationRedshiftOff);
    ret += sscanf(line, "RadiationRedshiftFullOn = %"FSYM, &CoolData.RadiationRedshiftFullOn);
    ret += sscanf(line, "RadiationRedshiftDropOff = %"FSYM, &CoolData.RadiationRedshiftDropOff);
    ret += sscanf(line, "HydrogenFractionByMass = %"FSYM, &CoolData.HydrogenFractionByMass);
    ret += sscanf(line, "DeuteriumToHydrogenRatio = %"FSYM, &CoolData.DeuteriumToHydrogenRatio);
    ret += sscanf(line, "SolarMetalFractionByMass = %"FSYM, &CoolData.SolarMetalFractionByMass);
    ret += sscanf(line, "NumberOfTemperatureBins = %"ISYM, &CoolData.NumberOfTemperatureBins);
    ret += sscanf(line, "CoolDataIh2co = %"ISYM, &CoolData.ih2co);
    ret += sscanf(line, "CoolDataIpiht = %"ISYM, &CoolData.ipiht);
    ret += sscanf(line, "TemperatureStart = %"FSYM, &CoolData.TemperatureStart);
    ret += sscanf(line, "TemperatureEnd = %"FSYM, &CoolData.TemperatureEnd);
    ret += sscanf(line, "CoolDataCompXray = %"FSYM, &CoolData.comp_xray);
    ret += sscanf(line, "CoolDataTempXray = %"FSYM, &CoolData.temp_xray);
    ret += sscanf(line, "RateDataCaseBRecombination = %"ISYM, &RateData.CaseBRecombination);
    ret += sscanf(line, "NumberOfDustTemperatureBins = %"ISYM, &RateData.NumberOfDustTemperatureBins);
    ret += sscanf(line, "DustTemperatureStart = %"FSYM, &RateData.DustTemperatureStart);
    ret += sscanf(line, "DustTemperatureEnd = %"FSYM, &RateData.DustTemperatureEnd);
    ret += sscanf(line, "PhotoelectricHeating  = %"ISYM, &PhotoelectricHeating);
    ret += sscanf(line, "PhotoelectricHeatingRate = %"FSYM, &PhotoelectricHeatingRate);

    ret += sscanf(line, "OutputCoolingTime = %"ISYM, &OutputCoolingTime);
    ret += sscanf(line, "OutputTemperature = %"ISYM, &OutputTemperature);
    ret += sscanf(line, "OutputDustTemperature = %"ISYM, &OutputDustTemperature);

    ret += sscanf(line, "OutputSmoothedDarkMatter = %"ISYM, 
		  &OutputSmoothedDarkMatter);
    ret += sscanf(line, "SmoothedDarkMatterNeighbors = %"ISYM, 
		  &SmoothedDarkMatterNeighbors);
    ret += sscanf(line, "OutputGriddedStarParticle = %"ISYM, 
		  &OutputGriddedStarParticle);

    ret += sscanf(line, "ZEUSQuadraticArtificialViscosity = %"FSYM,
		  &ZEUSQuadraticArtificialViscosity);
    ret += sscanf(line, "ZEUSLinearArtificialViscosity = %"FSYM,
		  &ZEUSLinearArtificialViscosity);
 
    ret += sscanf(line, "UseMinimumPressureSupport = %"ISYM,
		  &UseMinimumPressureSupport);
    ret += sscanf(line, "MinimumPressureSupportParameter = %"FSYM,
		  &MinimumPressureSupportParameter);
    ret += sscanf(line, "RefineByJeansLengthSafetyFactor = %"FSYM,
		  &RefineByJeansLengthSafetyFactor);
    ret += sscanf(line, "JeansRefinementColdTemperature = %"FSYM,
		  &JeansRefinementColdTemperature);
    ret += sscanf(line, "RefineByResistiveLengthSafetyFactor = %" FSYM,
		  &RefineByResistiveLengthSafetyFactor);
    ret += sscanf(line, "MustRefineParticlesRefineToLevel = %"ISYM,
                  &MustRefineParticlesRefineToLevel);
    ret += sscanf(line, "MustRefineParticlesCreateParticles = %"ISYM,
                  &MustRefineParticlesCreateParticles);
    ret += sscanf(line, "MustRefineParticlesLeftEdge  = %"PSYM" %"PSYM" %"PSYM,
                  MustRefineParticlesLeftEdge, MustRefineParticlesLeftEdge+1, 
                  MustRefineParticlesLeftEdge+2);
    ret += sscanf(line, "MustRefineParticlesRightEdge = %"PSYM" %"PSYM" %"PSYM,
                  MustRefineParticlesRightEdge, MustRefineParticlesRightEdge+1,
                  MustRefineParticlesRightEdge+2);
    ret += sscanf(line, "MustRefineParticlesRefineToLevelAutoAdjust = %"ISYM,
                  &MustRefineParticlesRefineToLevelAutoAdjust);
    ret += sscanf(line, "MustRefineParticlesMinimumMass = %"FSYM,
                  &MustRefineParticlesMinimumMass);
    ret += sscanf(line, "ParticleTypeInFile = %"ISYM,
                  &ParticleTypeInFile);

    if (sscanf(line, "AvoidRefineRegionLevel[%"ISYM"] = %"ISYM,&dim,&int_dummy) == 2){
      if (dim > MAX_STATIC_REGIONS-1) {
        ENZO_VFAIL("AvoidRegion number %"ISYM" (MAX_STATIC_REGIONS) > MAX allowed\n", dim)
      }
      ret++;
      AvoidRefineRegionLevel[dim] = int_dummy;
    }
    if (sscanf(line, "AvoidRefineRegionLeftEdge[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "AvoidRefineRegionLeftEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, AvoidRefineRegionLeftEdge[dim],
		    AvoidRefineRegionLeftEdge[dim]+1,
		    AvoidRefineRegionLeftEdge[dim]+2);
    if (sscanf(line, "AvoidRefineRegionRightEdge[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "AvoidRefineRegionRightEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, AvoidRefineRegionRightEdge[dim],
		    AvoidRefineRegionRightEdge[dim]+1,
		    AvoidRefineRegionRightEdge[dim]+2);
 
    if (sscanf(line, "StaticRefineRegionLevel[%"ISYM"] = %"ISYM,&dim,&int_dummy) == 2){
      if (dim > MAX_STATIC_REGIONS-1) {
        ENZO_VFAIL("StaticRegion number %"ISYM" > MAX allowed\n", dim)
      }
      ret++;
      StaticRefineRegionLevel[dim] = int_dummy;
    }
    if (sscanf(line, "StaticRefineRegionLeftEdge[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "StaticRefineRegionLeftEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, StaticRefineRegionLeftEdge[dim],
		    StaticRefineRegionLeftEdge[dim]+1,
		    StaticRefineRegionLeftEdge[dim]+2);
    if (sscanf(line, "StaticRefineRegionRightEdge[%"ISYM"] = ", &dim) == 1)
      ret += sscanf(line,
		    "StaticRefineRegionRightEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &dim, StaticRefineRegionRightEdge[dim],
		    StaticRefineRegionRightEdge[dim]+1,
		    StaticRefineRegionRightEdge[dim]+2);
 
    ret += sscanf(line, "ParallelRootGridIO = %"ISYM, &ParallelRootGridIO);
 
    ret += sscanf(line, "ParallelParticleIO = %"ISYM, &ParallelParticleIO);
 
    ret += sscanf(line, "Unigrid = %"ISYM, &Unigrid);
    ret += sscanf(line, "UnigridTranspose = %"ISYM, &UnigridTranspose);
    ret += sscanf(line, "NumberOfRootGridTilesPerDimensionPerProcessor = %"ISYM, &NumberOfRootGridTilesPerDimensionPerProcessor);
    ret += sscanf(line, "UserDefinedRootGridLayout = %"ISYM" %"ISYM" %"ISYM, &UserDefinedRootGridLayout[0],
                  &UserDefinedRootGridLayout[1], &UserDefinedRootGridLayout[2]);

    ret += sscanf(line, "PartitionNestedGrids = %"ISYM, &PartitionNestedGrids);
 
    ret += sscanf(line, "ExtractFieldsOnly = %"ISYM, &ExtractFieldsOnly);
 
    ret += sscanf(line, "CubeDumpEnabled = %"ISYM, &CubeDumpEnabled);
 
    ret += sscanf(line, "Debug1 = %"ISYM, &debug1);

    ret += sscanf(line, "Debug2 = %"ISYM, &debug2);

    ret += sscanf(line, "MemoryLimit = %lld", &MemoryLimit);

#ifdef STAGE_INPUT
    ret += sscanf(line, "StageInput = %"ISYM, &StageInput);
#endif

#ifdef OOC_BOUNDARY

    ret += sscanf(line, "ExternalBoundaryIO = %"ISYM, &ExternalBoundaryIO);

    ret += sscanf(line, "ExternalBoundaryTypeIO = %"ISYM, &ExternalBoundaryTypeIO);

    ret += sscanf(line, "ExternalBoundaryValueIO = %"ISYM, &ExternalBoundaryValueIO);

    ret += sscanf(line, "SimpleConstantBoundary = %"ISYM, &SimpleConstantBoundary);

#endif

    ret += sscanf(line, "SlopeFlaggingFields = "
		  " %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM,
		  SlopeFlaggingFields+0, 
		  SlopeFlaggingFields+1,
		  SlopeFlaggingFields+2, 
		  SlopeFlaggingFields+3,
		  SlopeFlaggingFields+4,
		  SlopeFlaggingFields+5,
		  SlopeFlaggingFields+6);
    
    ret += sscanf(line, "MinimumSlopeForRefinement = " 	  
		  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
	
		  MinimumSlopeForRefinement+0,
		  MinimumSlopeForRefinement+1,
		  MinimumSlopeForRefinement+2,
		  MinimumSlopeForRefinement+3,
		  MinimumSlopeForRefinement+4,  
		  MinimumSlopeForRefinement+5,
		  MinimumSlopeForRefinement+6);

     ret += sscanf(line, "SecondDerivativeFlaggingFields = "
		  " %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM,
		  SecondDerivativeFlaggingFields+0, 
		  SecondDerivativeFlaggingFields+1,
		  SecondDerivativeFlaggingFields+2, 
		  SecondDerivativeFlaggingFields+3,
		  SecondDerivativeFlaggingFields+4,
		  SecondDerivativeFlaggingFields+5,
		  SecondDerivativeFlaggingFields+6);
    
    ret += sscanf(line, "MinimumSecondDerivativeForRefinement = " 	  
		  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
	
		  MinimumSecondDerivativeForRefinement+0,
		  MinimumSecondDerivativeForRefinement+1,
		  MinimumSecondDerivativeForRefinement+2,
		  MinimumSecondDerivativeForRefinement+3,
		  MinimumSecondDerivativeForRefinement+4,  
		  MinimumSecondDerivativeForRefinement+5,
		  MinimumSecondDerivativeForRefinement+6);

    ret += sscanf(line, "SecondDerivativeEpsilon  = %"FSYM,
		  &SecondDerivativeEpsilon);

    ret += sscanf(line, "MinimumOverDensityForRefinement  = "
		  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
		  MinimumOverDensityForRefinement+0, 
		  MinimumOverDensityForRefinement+1,
		  MinimumOverDensityForRefinement+2, 
		  MinimumOverDensityForRefinement+3,
		  MinimumOverDensityForRefinement+4, 
		  MinimumOverDensityForRefinement+5,
		  MinimumOverDensityForRefinement+6);

    ret += sscanf(line, "MinimumMassForRefinement  = "
		  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
		  MinimumMassForRefinement+0, 
		  MinimumMassForRefinement+1,
		  MinimumMassForRefinement+2, 
		  MinimumMassForRefinement+3,
		  MinimumMassForRefinement+4, 
		  MinimumMassForRefinement+5,
		  MinimumMassForRefinement+6);

    ret += sscanf(line, "MinimumMassForRefinementLevelExponent = "
		  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
		  MinimumMassForRefinementLevelExponent+0,
		  MinimumMassForRefinementLevelExponent+1,
		  MinimumMassForRefinementLevelExponent+2,
		  MinimumMassForRefinementLevelExponent+3,
		  MinimumMassForRefinementLevelExponent+4,
		  MinimumMassForRefinementLevelExponent+5,
		  MinimumMassForRefinementLevelExponent+6);

    ret += sscanf(line, "MinimumPressureJumpForRefinement = %"FSYM,
		  &MinimumPressureJumpForRefinement);
    ret += sscanf(line, "OldShearMethod = %"ISYM,
		  &OldShearMethod);
    ret += sscanf(line, "MinimumShearForRefinement = %"FSYM,
		  &MinimumShearForRefinement);
    ret += sscanf(line, "MinimumEnergyRatioForRefinement = %"FSYM,
		  &MinimumEnergyRatioForRefinement);
    ret += sscanf(line, "ShockwaveRefinementMinMach = %"FSYM,
                 &ShockwaveRefinementMinMach);
    ret += sscanf(line, "ShockwaveRefinementMinVelocity = %"FSYM,
                 &ShockwaveRefinementMinVelocity);
    ret += sscanf(line, "ShockwaveRefinementMaxLevel = %"ISYM,
                 &ShockwaveRefinementMaxLevel);
    ret += sscanf(line, "ComovingCoordinates = %"ISYM,&ComovingCoordinates);
    ret += sscanf(line, "StarParticleCreation = %"ISYM, &StarParticleCreation);
    ret += sscanf(line, "BigStarFormation = %"ISYM, &BigStarFormation);
    ret += sscanf(line, "BigStarFormationDone = %"ISYM, &BigStarFormationDone);
    ret += sscanf(line, "BigStarSeparation = %"FSYM, &BigStarSeparation);
    ret += sscanf(line, "SimpleQ = %lf", &SimpleQ);
    ret += sscanf(line, "SimpleRampTime = %"FSYM, &SimpleRampTime);
    ret += sscanf(line, "StarFormationOncePerRootGridTimeStep = %"ISYM, &StarFormationOncePerRootGridTimeStep);
    ret += sscanf(line, "StarParticleFeedback = %"ISYM, &StarParticleFeedback);
    ret += sscanf(line, "StarParticleRadiativeFeedback = %"ISYM, &StarParticleRadiativeFeedback);
    ret += sscanf(line, "NumberOfParticleAttributes = %"ISYM,
		  &NumberOfParticleAttributes);

    /* read data which defines the boundary conditions */
 
    ret += sscanf(line, "LeftFaceBoundaryCondition  = %"ISYM" %"ISYM" %"ISYM,
		  MetaData.LeftFaceBoundaryCondition,
		  MetaData.LeftFaceBoundaryCondition+1,
		  MetaData.LeftFaceBoundaryCondition+2);
    
     ret += sscanf(line, "RightFaceBoundaryCondition = %"ISYM" %"ISYM" %"ISYM,
 		  MetaData.RightFaceBoundaryCondition,
 		  MetaData.RightFaceBoundaryCondition+1,
 		  MetaData.RightFaceBoundaryCondition+2);
  
    if (sscanf(line, "BoundaryConditionName         = %s", dummy) == 1)
      MetaData.BoundaryConditionName = dummy;

    if (sscanf(line, "MetaDataIdentifier = %s", dummy) == 1) {
      MetaData.MetaDataIdentifier = dummy;
      ret++;
    }
    if (sscanf(line, "MetaDataSimulationUUID = %s", dummy) == 1) {
      MetaData.SimulationUUID = dummy;
      ret++;
    }
    if (sscanf(line, "MetaDataDatasetUUID = %s", dummy) == 1) {
      MetaData.RestartDatasetUUID = dummy;
      ret++;
    }
    if (sscanf(line, "MetaDataInitialConditionsUUID = %s", dummy) == 1) {
      MetaData.InitialConditionsUUID = dummy;
      ret++;
    }
 
    /* Check version number. */
 
    if (sscanf(line, "VersionNumber = %"FSYM, &TempFloat) == 1) {
      ret++;
      if (fabs(TempFloat - VERSION) >= 1.0e-3 &&
	  MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Warning: Incorrect version number.\n");
    }

    /* Read Galaxy Simulation Wind Boundary Variabels */

     ret += sscanf(line, "GalaxySimulationRPSWind = %"ISYM,&GalaxySimulationRPSWind);
     ret += sscanf(line, "GalaxySimulationRPSWindShockSpeed = %"FSYM,&GalaxySimulationRPSWindShockSpeed);
     ret += sscanf(line, "GalaxySimulationRPSWindDelay = %"FSYM,&GalaxySimulationRPSWindDelay);
     ret += sscanf(line, "GalaxySimulationRPSWindDensity = %"FSYM,&GalaxySimulationRPSWindDensity);
     ret += sscanf(line, "GalaxySimulationRPSWindTotalEnergy = %"FSYM,&GalaxySimulationRPSWindTotalEnergy);
     ret += sscanf(line, "GalaxySimulationRPSWindPressure = %"FSYM,&GalaxySimulationRPSWindPressure);
     ret += sscanf(line, "GalaxySimulationRPSWindVelocity = %"PSYM" %"PSYM" %"PSYM,
      GalaxySimulationRPSWindVelocity, GalaxySimulationRPSWindVelocity+1, GalaxySimulationRPSWindVelocity+2);
     ret += sscanf(line, "GalaxySimulationPreWindDensity = %"FSYM,&GalaxySimulationPreWindDensity);
     ret += sscanf(line, "GalaxySimulationPreWindTotalEnergy = %"FSYM,&GalaxySimulationPreWindTotalEnergy);
     ret += sscanf(line, "GalaxySimulationPreWindVelocity = %"PSYM" %"PSYM" %"PSYM,
         GalaxySimulationPreWindVelocity,GalaxySimulationPreWindVelocity+1,GalaxySimulationPreWindVelocity+2);
 
    /* Read star particle parameters. */

    ret += sscanf(line, "StarMakerTypeIaSNe = %"ISYM,
		  &StarMakerTypeIaSNe);
    ret += sscanf(line, "StarMakerTypeIISNeMetalField = %"ISYM,
		  &StarMakerTypeIISNeMetalField);
    ret += sscanf(line, "StarMakerPlanetaryNebulae = %"ISYM,
		  &StarMakerPlanetaryNebulae);
    ret += sscanf(line, "StarMakerOverDensityThreshold = %"FSYM,
		  &StarMakerOverDensityThreshold);
    ret += sscanf(line, "StarMakerUseOverDensityThreshold = %"ISYM,
          &StarMakerUseOverDensityThreshold);
    ret += sscanf(line, "StarMakerMaximumFractionCell = %"FSYM,
          &StarMakerMaximumFractionCell);
    ret += sscanf(line, "StarMakerSHDensityThreshold = %"FSYM,
		  &StarMakerSHDensityThreshold);
    ret += sscanf(line, "StarMakerTimeIndependentFormation = %"ISYM,
		  &StarMakerTimeIndependentFormation);
    ret += sscanf(line, "StarMakerMassEfficiency = %"FSYM,
		  &StarMakerMassEfficiency);
    ret += sscanf(line, "StarMakerMinimumMass = %"FSYM, &StarMakerMinimumMass);
    ret += sscanf(line, "StarMakerMinimumDynamicalTime = %"FSYM,
                  &StarMakerMinimumDynamicalTime);
    ret += sscanf(line, "StarMassEjectionFraction = %"FSYM,
		  &StarMassEjectionFraction);
    ret += sscanf(line, "StarMetalYield = %"FSYM, &StarMetalYield);
    ret += sscanf(line, "StarEnergyToThermalFeedback = %"FSYM,
		  &StarEnergyToThermalFeedback);
    ret += sscanf(line, "StarEnergyToStellarUV = %"FSYM, &StarEnergyToStellarUV);
    ret += sscanf(line, "StarEnergyToQuasarUV = %"FSYM, &StarEnergyToQuasarUV);
    ret += sscanf(line, "StarFeedbackDistRadius = %"ISYM, &StarFeedbackDistRadius);
    ret += sscanf(line, "StarFeedbackDistCellStep = %"ISYM, &StarFeedbackDistCellStep);

    ret += sscanf(line, "StarClusterUseMetalField = %"ISYM, 
		  &StarClusterUseMetalField);
    ret += sscanf(line, "StarClusterMinDynamicalTime = %"FSYM, 
		  &StarClusterMinDynamicalTime);
    ret += sscanf(line, "StarClusterHeliumIonization = %"ISYM, 
		  &StarClusterHeliumIonization);
    ret += sscanf(line, "StarClusterUnresolvedModel = %"ISYM, 
		  &StarClusterUnresolvedModel);
    ret += sscanf(line, "StarClusterIonizingLuminosity = %lf", 
		  &StarClusterIonizingLuminosity);
    ret += sscanf(line, "StarClusterSNEnergy = %lf", &StarClusterSNEnergy);
    ret += sscanf(line, "StarClusterSNRadius = %"FSYM, &StarClusterSNRadius);
    ret += sscanf(line, "StarClusterFormEfficiency = %"FSYM, 
		  &StarClusterFormEfficiency);
    ret += sscanf(line, "StarClusterMinimumMass = %"FSYM, 
		  &StarClusterMinimumMass);
    ret += sscanf(line, "StarClusterCombineRadius = %"FSYM,
		  &StarClusterCombineRadius);
    ret += sscanf(line, "StarClusterRegionLeftEdge = %"FSYM" %"FSYM" %"FSYM,
		  StarClusterRegionLeftEdge, StarClusterRegionLeftEdge+1, 
		  StarClusterRegionLeftEdge+2);
    ret += sscanf(line, "StarClusterRegionRightEdge = %"FSYM" %"FSYM" %"FSYM,
		  StarClusterRegionRightEdge, StarClusterRegionRightEdge+1, 
		  StarClusterRegionRightEdge+2);

    ret += sscanf(line, "PopIIIStarMass = %"FSYM, &PopIIIStarMass);
    ret += sscanf(line, "PopIIIInitialMassFunction = %"ISYM, 
		  &PopIIIInitialMassFunction);
    ret += sscanf(line, "PopIIIInitialMassFunctionSeed = %"ISYM, 
		  &PopIIIInitialMassFunctionSeed);
    ret += sscanf(line, "PopIIIInitialMassFunctionCalls = %"ISYM, 
		  &PopIIIInitialMassFunctionCalls);
    ret += sscanf(line, "PopIIIMassRange = %"FSYM" %"FSYM,
		  &PopIIILowerMassCutoff, &PopIIIUpperMassCutoff);
    ret += sscanf(line, "PopIIIInitialMassFunctionSlope = %"FSYM, 
		  &PopIIIInitialMassFunctionSlope);
    ret += sscanf(line, "PopIIIBlackHoles = %"ISYM, &PopIIIBlackHoles);
    ret += sscanf(line, "PopIIIBHLuminosityEfficiency = %"FSYM, 
		  &PopIIIBHLuminosityEfficiency);
    ret += sscanf(line, "PopIIIOverDensityThreshold = %"FSYM,
		  &PopIIIOverDensityThreshold);
    ret += sscanf(line, "PopIIIH2CriticalFraction = %"FSYM,
		  &PopIIIH2CriticalFraction);
    ret += sscanf(line, "PopIIIMetalCriticalFraction = %"FSYM,
		  &PopIIIMetalCriticalFraction);
    ret += sscanf(line, "PopIIISupernovaRadius = %"FSYM, &PopIIISupernovaRadius);
    ret += sscanf(line, "PopIIISupernovaUseColour = %"ISYM, 
		  &PopIIISupernovaUseColour);
    ret += sscanf(line, "PopIIISupernovaMustRefine = %"ISYM,
		  &PopIIISupernovaMustRefine);
    ret += sscanf(line, "PopIIISupernovaMustRefineResolution = %"ISYM,
		  &PopIIISupernovaMustRefineResolution);
    ret += sscanf(line, "PopIIIHeliumIonization = %"ISYM, 
		  &PopIIIHeliumIonization);

    ret += sscanf(line, "PopIIIColorDensityThreshold = %"FSYM,
		  &PopIIIColorDensityThreshold);
    ret += sscanf(line, "PopIIIColorMass = %"FSYM,
		  &PopIIIColorMass);
    ret += sscanf(line, "PopIIIUseHypernova = %"ISYM,
		  &PopIIIUseHypernova);
    ret += sscanf(line, "PopIIISupernovaExplosions = %"ISYM,
		  &PopIIISupernovaExplosions);
    ret += sscanf(line, "PopIIIOutputOnFeedback = %"ISYM,
		  &PopIIIOutputOnFeedback);

    ret += sscanf(line, "MBHAccretion = %"ISYM, &MBHAccretion);
    ret += sscanf(line, "MBHAccretionRadius = %"FSYM, &MBHAccretionRadius);
    ret += sscanf(line, "MBHAccretingMassRatio = %"FSYM, &MBHAccretingMassRatio);
    ret += sscanf(line, "MBHAccretionFixedTemperature = %"FSYM, &MBHAccretionFixedTemperature);
    ret += sscanf(line, "MBHAccretionFixedRate = %"FSYM, &MBHAccretionFixedRate);
    ret += sscanf(line, "MBHTurnOffStarFormation = %"ISYM, &MBHTurnOffStarFormation);
    ret += sscanf(line, "MBHCombineRadius = %"FSYM, &MBHCombineRadius);
    ret += sscanf(line, "MBHMinDynamicalTime = %"FSYM, &MBHMinDynamicalTime);
    ret += sscanf(line, "MBHMinimumMass = %"FSYM, &MBHMinimumMass);

    ret += sscanf(line, "MBHFeedback = %"ISYM, &MBHFeedback);
    ret += sscanf(line, "MBHFeedbackRadiativeEfficiency = %"FSYM, &MBHFeedbackRadiativeEfficiency);
    ret += sscanf(line, "MBHFeedbackEnergyCoupling = %"FSYM, &MBHFeedbackEnergyCoupling);
    ret += sscanf(line, "MBHFeedbackMassEjectionFraction = %"FSYM, &MBHFeedbackMassEjectionFraction);
    ret += sscanf(line, "MBHFeedbackMetalYield = %"FSYM, &MBHFeedbackMetalYield);
    ret += sscanf(line, "MBHFeedbackThermalRadius = %"FSYM, &MBHFeedbackThermalRadius);
    ret += sscanf(line, "MBHFeedbackJetsThresholdMass = %"FSYM, &MBHFeedbackJetsThresholdMass);

    ret += sscanf(line, "MBHParticleIO = %"ISYM,
		  &MBHParticleIO);
    if (sscanf(line, "MBHParticleIOFilename = %s", dummy) == 1)
      MBHParticleIOFilename = dummy;
    if (sscanf(line, "MBHInsertLocationFilename = %s", dummy) == 1)
      MBHInsertLocationFilename = dummy;

    ret += sscanf(line, "H2StarMakerEfficiency = %"FSYM,
		  &H2StarMakerEfficiency);
    ret += sscanf(line, "H2StarMakerNumberDensityThreshold = %"FSYM,
		  &H2StarMakerNumberDensityThreshold);
    ret += sscanf(line, "H2StarMakerMinimumMass = %"FSYM,
		  &H2StarMakerMinimumMass);
    ret += sscanf(line, "H2StarMakerMinimumH2FractionForStarFormation = %"FSYM,
		  &H2StarMakerMinimumH2FractionForStarFormation);
    ret += sscanf(line, "H2StarMakerStochastic = %"ISYM,
		  &H2StarMakerStochastic);
    ret += sscanf(line, "H2StarMakerUseSobolevColumn = %"ISYM,
		  &H2StarMakerUseSobolevColumn);
    ret += sscanf(line, "H2StarMakerSigmaOverR = %"FSYM,
		  &H2StarMakerSigmaOverR);
    ret += sscanf(line, "H2StarMakerAssumeColdWarmPressureBalance = %"ISYM,
		  &H2StarMakerAssumeColdWarmPressureBalance);
    ret += sscanf(line, "H2StarMakerH2DissociationFlux_MW = %"FSYM,
		  &H2StarMakerH2DissociationFlux_MW);
    ret += sscanf(line, "H2StarMakerH2FloorInColdGas = %"FSYM,
		  &H2StarMakerH2FloorInColdGas);
    ret += sscanf(line, "H2StarMakerColdGasTemperature = %"FSYM,
		  &H2StarMakerColdGasTemperature);

    /* Read Movie Dump parameters */

    ret += sscanf(line, "MovieSkipTimestep = %"ISYM, &MovieSkipTimestep);
    ret += sscanf(line, "Movie3DVolumes = %"ISYM, &Movie3DVolumes);
    ret += sscanf(line, "MovieVertexCentered = %"ISYM, &MovieVertexCentered);
    ret += sscanf(line, "NewMovieParticleOn = %"ISYM, &NewMovieParticleOn);
    ret += sscanf(line, "MovieDataField = %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM,
		  MovieDataField+0, MovieDataField+1, MovieDataField+2,
		  MovieDataField+3, MovieDataField+4, MovieDataField+5);
    ret += sscanf(line, "NewMovieDumpNumber = %"ISYM, &NewMovieDumpNumber);
    if (sscanf(line, "NewMovieName = %s", dummy) == 1)
      NewMovieName = dummy;
    ret += sscanf(line, "MovieTimestepCounter = %"ISYM, &MetaData.MovieTimestepCounter);

    ret += sscanf(line, "MultiMetals = %"ISYM, &MultiMetals);
    ret += sscanf(line, "IsotropicConduction = %"ISYM, &IsotropicConduction);
    ret += sscanf(line, "AnisotropicConduction = %"ISYM, &AnisotropicConduction);
    ret += sscanf(line, "IsotropicConductionSpitzerFraction = %"FSYM, &IsotropicConductionSpitzerFraction);
    ret += sscanf(line, "AnisotropicConductionSpitzerFraction = %"FSYM, &AnisotropicConductionSpitzerFraction);
    ret += sscanf(line, "ConductionCourantSafetyNumber = %"FSYM, &ConductionCourantSafetyNumber);
    ret += sscanf(line, "SpeedOfLightTimeStepLimit = %"ISYM, &SpeedOfLightTimeStepLimit);

    ret += sscanf(line, "RadiativeTransfer = %"ISYM, &RadiativeTransfer);
    ret += sscanf(line, "RadiationXRaySecondaryIon = %"ISYM, &RadiationXRaySecondaryIon);
    ret += sscanf(line, "RadiationXRayComptonHeating = %"ISYM, &RadiationXRayComptonHeating);

    /* Shearing Box Boundary parameters */

    ret += sscanf(line, "AngularVelocity = %"FSYM, &AngularVelocity);
    ret += sscanf(line, "VelocityGradient = %"FSYM, &VelocityGradient);
    ret += sscanf(line, "ShearingVelocityDirection = %"ISYM, &ShearingVelocityDirection);
    ret += sscanf(line, "ShearingBoxProblemType = %"ISYM, &ShearingBoxProblemType);  
    

#ifdef STAGE_INPUT
    sscanf(line, "LocalPath = %s\n", LocalPath);
    sscanf(line, "GlobalPath = %s\n", GlobalPath);
#endif

    /* Embedded Python */
    ret += sscanf(line, "PythonTopGridSkip = %"ISYM, &PythonTopGridSkip);
    ret += sscanf(line, "PythonSubcycleSkip = %"ISYM, &PythonSubcycleSkip);
    ret += sscanf(line, "PythonReloadScript = %"ISYM, &PythonReloadScript);
#ifdef USE_PYTHON
    ret += sscanf(line, "NumberOfPythonCalls = %"ISYM, &NumberOfPythonCalls);
    ret += sscanf(line, "NumberOfPythonTopGridCalls = %"ISYM, &NumberOfPythonTopGridCalls);
    ret += sscanf(line, "NumberOfPythonSubcycleCalls = %"ISYM, &NumberOfPythonSubcycleCalls);
#endif

    /* EnzoTiming Parameters */
    ret += sscanf(line, "TimingCycleSkip = %"ISYM, &TimingCycleSkip);

    /* Inline halo finder */

    ret += sscanf(line, "InlineHaloFinder = %"ISYM, &InlineHaloFinder);
    ret += sscanf(line, "HaloFinderSubfind = %"ISYM, &HaloFinderSubfind);
    ret += sscanf(line, "HaloFinderOutputParticleList = %"ISYM, 
		  &HaloFinderOutputParticleList);
    ret += sscanf(line, "HaloFinderRunAfterOutput = %"ISYM, 
		  &HaloFinderRunAfterOutput);
    ret += sscanf(line, "HaloFinderLinkingLength = %"FSYM, 
		  &HaloFinderLinkingLength);
    ret += sscanf(line, "HaloFinderMinimumSize = %"ISYM, &HaloFinderMinimumSize);
    ret += sscanf(line, "HaloFinderCycleSkip = %"ISYM, &HaloFinderCycleSkip);
    ret += sscanf(line, "HaloFinderTimestep = %"FSYM, &HaloFinderTimestep);
    ret += sscanf(line, "HaloFinderLastTime = %"PSYM, &HaloFinderLastTime);

    /* This Block for Stanford Hydro */

    ret += sscanf(line, "UseHydro               = %"ISYM, &UseHydro);


    /* Sink particles (for present day star formation) & winds */
    ret += sscanf(line, "SinkMergeDistance     = %"FSYM, &SinkMergeDistance); 
    ret += sscanf(line, "SinkMergeMass         = %"FSYM, &SinkMergeMass);
    ret += sscanf(line, "StellarWindFeedback   = %"ISYM, &StellarWindFeedback);
    ret += sscanf(line, "StellarWindTurnOnMass = %"FSYM, &StellarWindTurnOnMass);
    ret += sscanf(line, "MSStellarWindTurnOnMass = %"FSYM, &MSStellarWindTurnOnMass);

    ret += sscanf(line, "VelAnyl = %"ISYM, &VelAnyl);
    ret += sscanf(line, "BAnyl = %"ISYM, &BAnyl);
    ret += sscanf(line, "WriteExternalAccel = %"ISYM, &WriteExternalAccel);


    /* Read MHD Paramters */
    ret += sscanf(line, "UseDivergenceCleaning = %"ISYM"", &UseDivergenceCleaning);
    ret += sscanf(line, "DivergenceCleaningBoundaryBuffer = %"ISYM"", &DivergenceCleaningBoundaryBuffer);
    ret += sscanf(line, "DivergenceCleaningThreshold = %"FSYM, &DivergenceCleaningThreshold);
    ret += sscanf(line, "PoissonApproximationThreshold = %"FSYM, &PoissonApproximationThreshold);
    ret += sscanf(line, "PoissonBoundaryType = %"ISYM"", &PoissonBoundaryType);
   

    ret += sscanf(line, "AngularVelocity = %"FSYM, &AngularVelocity);
    ret += sscanf(line, "VelocityGradient = %"FSYM, &VelocityGradient);
    ret += sscanf(line, "UseDrivingField = %"ISYM"", &UseDrivingField);
    ret += sscanf(line, "DrivingEfficiency = %"FSYM, &DrivingEfficiency);

    ret += sscanf(line, "StringKick = %"FSYM, &StringKick);
    ret += sscanf(line, "StringKickDimension = %"ISYM, &StringKickDimension);
    ret += sscanf(line, "UsePhysicalUnit = %"ISYM"", &UsePhysicalUnit);
    ret += sscanf(line, "Theta_Limiter = %"FSYM, &Theta_Limiter);
    ret += sscanf(line, "RKOrder = %"ISYM"", &RKOrder);
    ret += sscanf(line, "UseFloor = %"ISYM"", &UseFloor);
    ret += sscanf(line, "UseViscosity = %"ISYM"", &UseViscosity);
    ret += sscanf(line, "ViscosityCoefficient = %"FSYM, &ViscosityCoefficient);  
    ret += sscanf(line, "UseAmbipolarDiffusion = %"ISYM"", &UseAmbipolarDiffusion);
    ret += sscanf(line, "UseResistivity = %"ISYM"", &UseResistivity);
    ret += sscanf(line, "SmallRho = %"FSYM, &SmallRho);
    ret += sscanf(line, "SmallP = %"FSYM, &SmallP);
    ret += sscanf(line, "SmallT = %"FSYM, &SmallT);
    ret += sscanf(line, "MaximumAlvenSpeed = %"FSYM, &MaximumAlvenSpeed);
    ret += sscanf(line, "Coordinate = %"ISYM, &Coordinate);
    ret += sscanf(line, "RiemannSolver = %"ISYM, &RiemannSolver);
    ret += sscanf(line, "RiemannSolverFallback = %"ISYM, &RiemannSolverFallback);
    ret += sscanf(line, "ConservativeReconstruction = %"ISYM, &ConservativeReconstruction);
    ret += sscanf(line, "PositiveReconstruction = %"ISYM, &PositiveReconstruction);
    ret += sscanf(line, "ReconstructionMethod = %"ISYM, &ReconstructionMethod);

    ret += sscanf(line, "EOSType = %"ISYM, &EOSType);
    ret += sscanf(line, "EOSSoundSpeed = %"FSYM, &EOSSoundSpeed);
    ret += sscanf(line, "EOSCriticalDensity = %"FSYM, &EOSCriticalDensity);
    ret += sscanf(line, "EOSGamma = %"FSYM, &EOSGamma);
    ret += sscanf(line, "UseConstantAcceleration = %"ISYM, &UseConstantAcceleration);
    ret += sscanf(line, "ConstantAcceleration = %"GSYM" %"GSYM" %"GSYM, &ConstantAcceleration[0],
		  &ConstantAcceleration[1], &ConstantAcceleration[2]);
    ret += sscanf(line, "Mu = %"FSYM, &Mu);
    ret += sscanf(line, "DivBDampingLength = %"FSYM, &DivBDampingLength);
    ret += sscanf(line, "CoolingCutOffDensity1 = %"GSYM, &CoolingCutOffDensity1);
    ret += sscanf(line, "CoolingCutOffDensity2 = %"GSYM, &CoolingCutOffDensity2);
    ret += sscanf(line, "CoolingCutOffTemperature = %"GSYM, &CoolingCutOffTemperature);
    ret += sscanf(line, "CoolingPowerCutOffDensity1 = %"GSYM, &CoolingPowerCutOffDensity1);
    ret += sscanf(line, "CoolingPowerCutOffDensity2 = %"GSYM, &CoolingPowerCutOffDensity2);
    ret += sscanf(line, "UseCUDA = %"ISYM,&UseCUDA);
    ret += sscanf(line, "ClusterSMBHFeedback = %"ISYM, &ClusterSMBHFeedback);
    ret += sscanf(line, "ClusterSMBHJetMdot = %"FSYM, &ClusterSMBHJetMdot);
    ret += sscanf(line, "ClusterSMBHJetVelocity = %"FSYM, &ClusterSMBHJetVelocity);
    ret += sscanf(line, "ClusterSMBHJetRadius = %"FSYM, &ClusterSMBHJetRadius);
    ret += sscanf(line, "ClusterSMBHJetLaunchOffset = %"FSYM, &ClusterSMBHJetLaunchOffset);
    ret += sscanf(line, "ClusterSMBHStartTime = %"FSYM, &ClusterSMBHStartTime);
    ret += sscanf(line, "ClusterSMBHTramp = %"FSYM, &ClusterSMBHTramp);
    ret += sscanf(line, "ClusterSMBHJetOpenAngleRadius = %"FSYM, &ClusterSMBHJetOpenAngleRadius);
    ret += sscanf(line, "ClusterSMBHFastJetRadius = %"FSYM, &ClusterSMBHFastJetRadius);
    ret += sscanf(line, "ClusterSMBHFastJetVelocity = %"FSYM, &ClusterSMBHFastJetVelocity);
    ret += sscanf(line, "ClusterSMBHJetEdot = %"FSYM, &ClusterSMBHJetEdot);
    ret += sscanf(line, "ClusterSMBHKineticFraction = %"FSYM, &ClusterSMBHKineticFraction);
    ret += sscanf(line, "ClusterSMBHJetAngleTheta = %"FSYM, &ClusterSMBHJetAngleTheta);
    ret += sscanf(line, "ClusterSMBHJetAnglePhi = %"FSYM, &ClusterSMBHJetAnglePhi);
    ret += sscanf(line, "ClusterSMBHJetPrecessionPeriod = %"FSYM, &ClusterSMBHJetPrecessionPeriod);
    ret += sscanf(line, "ClusterSMBHCalculateGasMass = %"ISYM, &ClusterSMBHCalculateGasMass);
    ret += sscanf(line, "ClusterSMBHFeedbackSwitch = %"ISYM, &ClusterSMBHFeedbackSwitch);
    ret += sscanf(line, "ClusterSMBHEnoughColdGas = %"FSYM, &ClusterSMBHEnoughColdGas);
    ret += sscanf(line, "ClusterSMBHAccretionTime = %"FSYM, &ClusterSMBHAccretionTime);
    ret += sscanf(line, "ClusterSMBHJetDim = %"ISYM, &ClusterSMBHJetDim);
    ret += sscanf(line, "ClusterSMBHAccretionEpsilon = %"FSYM, &ClusterSMBHAccretionEpsilon);

    ret += sscanf(line, "ExtraOutputs = %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"", ExtraOutputs,
		  ExtraOutputs +1,ExtraOutputs +2,ExtraOutputs +3,
		  ExtraOutputs +4,ExtraOutputs +5,ExtraOutputs +6,
		  ExtraOutputs +7,ExtraOutputs +8,ExtraOutputs +9);

    //MHDCT variables
    ret += sscanf(line, "MHDCTPowellSource             = %"ISYM, &MHDCTPowellSource);
    ret += sscanf(line, "MHDCTDualEnergyMethod             = %"ISYM, &MHDCTDualEnergyMethod);
    ret += sscanf(line, "MHDCTSlopeLimiter             = %"ISYM, &MHDCTSlopeLimiter);
    ret += sscanf(line, "WriteBoundary          = %"ISYM, &WriteBoundary);
    ret += sscanf(line,"TracerParticlesAddToRestart = %"ISYM,&TracerParticlesAddToRestart);
    ret += sscanf(line,"RefineByJeansLengthUnits = %"ISYM,&RefineByJeansLengthUnits);

    ret += sscanf(line,"CT_AthenaDissipation = %"FSYM,&CT_AthenaDissipation);
    ret += sscanf(line,"MHD_WriteElectric = %"ISYM,&MHD_WriteElectric);

    ret += sscanf(line,"tiny_pressure = %"FSYM,&tiny_pressure);
    ret += sscanf(line,"MHD_CT_Method = %"ISYM,&MHD_CT_Method);
		  
    ret += sscanf(line,"NumberOfGhostZones = %"ISYM,&NumberOfGhostZones);
    ret += sscanf(line,"MHD_ProjectB = %"ISYM,&MHD_ProjectB);
    ret += sscanf(line,"MHD_ProjectE = %"ISYM,&MHD_ProjectE);
    ret += sscanf(line,"EquationOfState = %"ISYM,&EquationOfState);
    if(sscanf(line, "MHDLabel[%"ISYM"] = %s\n", &dim, dummy) == 2)
      MHDLabel[dim] = dummy;
    if(sscanf(line, "MHDUnits[%"ISYM"] = %s\n", &dim, dummy) == 2)
      MHDUnits[dim] = dummy;
    if(sscanf(line, "MHDcLabel[%"ISYM"] = %s\n", &dim, dummy) == 2){
        ENZO_FAIL("Looks like you're restarting an OLD MHDCT run. \n Run src/CenteredBremover.py on your dataset.\n");
    }
    if(sscanf(line, "MHDeLabel[%"ISYM"] = %s\n", &dim, dummy) ==2)
      MHDeLabel[dim] = dummy;
    if(sscanf(line, "MHDeUnits[%"ISYM"] = %s\n", &dim, dummy) == 2)
      MHDeUnits[dim] = dummy;

    ret += sscanf(line, "CorrectParentBoundaryFlux             = %"ISYM, &CorrectParentBoundaryFlux);
    ret += sscanf(line, "MoveParticlesBetweenSiblings = %"ISYM,
		  &MoveParticlesBetweenSiblings);
    ret += sscanf(line, "ParticleSplitterIterations = %"ISYM,
		  &ParticleSplitterIterations);
    ret += sscanf(line, "ParticleSplitterRandomSeed = %"ISYM,
		  &ParticleSplitterRandomSeed);
    ret += sscanf(line, "ParticleSplitterChildrenParticleSeparation = %"FSYM,
		  &ParticleSplitterChildrenParticleSeparation);
    ret += sscanf(line, "ResetMagneticField = %"ISYM,
		  &ResetMagneticField);
    ret += sscanf(line, "ResetMagneticFieldAmplitude  =  %"GSYM" %"GSYM" %"GSYM, 
		  ResetMagneticFieldAmplitude,
		  ResetMagneticFieldAmplitude+1,
		  ResetMagneticFieldAmplitude+2);

    ret += sscanf(line, "UseGasDrag = %"ISYM, &UseGasDrag);
    ret += sscanf(line, "GasDragCoefficient = %"GSYM, &GasDragCoefficient);


    /* If the dummy char space was used, then make another. */
 
    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      dummy[0] = 0;
      ret++;
    }
 
    /* check to see if the line belongs to one of the test problems */
 
    if (strstr(line, "ShockTube")           ) ret++;
    if (strstr(line, "WavePool" )           ) ret++;
    if (strstr(line, "ShockPool")           ) ret++;
    if (strstr(line, "DoubleMach")          ) ret++;
    if (strstr(line, "Implosion")           ) ret++;
    if (strstr(line, "SedovBlast")          ) ret++;
    if (strstr(line, "Units")               ) ret++;
    if (strstr(line, "RadiatingShock")      ) ret++;
    if (strstr(line, "RotatingCylinder")    ) ret++;
    if (strstr(line, "RotatingDisk")    ) ret++;
    if (strstr(line, "RotatingSphere")    ) ret++;
    if (strstr(line, "StratifiedMediumExplosion")) ret++;
    if (strstr(line, "TestOrbit")    ) ret++;
    if (strstr(line, "KelvinHelmholtz")     ) ret++;
    if (strstr(line, "KH")                  ) ret++;
    if (strstr(line, "Noh")                 ) ret++;
    if (strstr(line, "TestProblem")         ) ret++;
    if (strstr(line, "ZeldovichPancake")    ) ret++;
    if (strstr(line, "PressurelessCollapse")) ret++;
    if (strstr(line, "AdiabaticExpansion" ) ) ret++;
    if (strstr(line, "CosmologySimulation") ) ret++;
    if (strstr(line, "TestGravity"        ) ) ret++;
    if (strstr(line, "SphericalInfall"    ) ) ret++;
    if (strstr(line, "TestGravitySphere"  ) ) ret++;
    if (strstr(line, "Cluster"            ) ) ret++;
    if (strstr(line, "CollapseTest"       ) ) ret++;
    if (strstr(line, "Cosmology"          ) ) ret++;
    if (strstr(line, "SupernovaRestart"   ) ) ret++;
    if (strstr(line, "TracerParticleCreation")) ret++;
    if (strstr(line, "TurbulenceSimulation")) ret++;
    if (strstr(line, "ProtostellarCollapse")) ret++;
    if (strstr(line, "GalaxySimulation") 
			&& !strstr(line,"RPSWind") && !strstr(line,"PreWind") ) ret++;
    if (strstr(line, "AgoraRestart")) ret++;
    if (strstr(line, "ConductionTest")) ret++;
    if (strstr(line, "ConductionBubble")) ret++;
    if (strstr(line, "ConductionCloud")) ret++;
    if (strstr(line, "CoolingTest")) ret++;
    if (strstr(line, "OneZoneFreefall")) ret++;
    if (strstr(line, "ShearingBox")) ret++;
    if (strstr(line, "PoissonSolverTest")) ret++;
    /* 7.22.10 - CBH: Added 5 following lines to avoid runtime warnings from 
    extra params previously added to code (but not read_params) by others.*/
    if (strstr(line, "Cloudy")              ) ret++;
    if (strstr(line, "IsothermalSoundSpeed")) ret++;
    if (strstr(line, "dtPhoton")            ) ret++;
    if (strstr(line, "CurrentTimeIdentifier")) ret++;
    if (strstr(line, "MetaDataRestart")     ) ret++;
    if (strstr(line, "MustRefine") ) ret++;
    if (strstr(line, "AccretionKernal")     ) ret++;
    if (strstr(line, "PopIII")              ) ret++;
#ifdef TRANSFER
    if (strstr(line, "Radiative")           ) ret++;
    if (strstr(line, "PhotonTest")          ) ret++;
#endif
    if (strstr(line, "MHDDRF")              ) ret++;
    if (strstr(line, "DrivenFlowMach")      ) ret++;
    if (strstr(line, "DrivenFlowMagField")  ) ret++;
    if (strstr(line, "DrivenFlowDensity")      ) ret++;
    if (strstr(line, "DrivenFlowPressure")      ) ret++;

    if (strstr(line, "\"\"\"")              ) comment_count++;

    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#')
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s", line);
 
  }

  // HierarchyFile IO sanity check

  // Note that although I only do not allow HierarchyFileInputFormat=2
  // (both ASCII and HDF5 input), it is supported internally for
  // debugging purpose.
  if ((HierarchyFileInputFormat < 0) || (HierarchyFileInputFormat > 1))
    ENZO_FAIL("Invalid HierarchyFileInputFormat. Must be 0 (HDF5) or 1 (ASCII).")

  if ((HierarchyFileOutputFormat < 0) || (HierarchyFileOutputFormat > 2))
    ENZO_FAIL("Invalid HierarchyFileOutputFormat. Must be 0 (HDF5), 1 (ASCII), or 2 (both).")
  
  // While we're examining the hierarchy, check that the MultiRefinedRegion doesn't demand more refinement that we've got                                                                                     
  for (int ireg = 0; ireg < MAX_STATIC_REGIONS; ireg++)
    if (MultiRefineRegionGeometry[ireg] >= 0)
      if (MultiRefineRegionMaximumLevel[ireg] > MaximumRefinementLevel)
	ENZO_VFAIL("MultiRefineRegionMaximumLevel[%"ISYM"] = %"ISYM"  > MaximumRefinementLevel\n", ireg, MultiRefineRegionMaximumLevel[ireg]);



  /* clean up */
 
  delete [] dummy;
  rewind(fptr);


/*  If stochastic forcing is used, initialize the object 
 *  but only if it is not a fresh simulation*/
  if (DrivenFlowProfile && MetaData.Time != 0.) {
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
                 DrivenFlowDomainLength[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

    Forcing.Init(MetaData.TopGridRank,
                 DrivenFlowProfile,
                 DrivenFlowAlpha,
                 DrivenFlowDomainLength,
                 DrivenFlowBandWidth,
                 DrivenFlowVelocity,
                 DrivenFlowAutoCorrl,
                 DrivenFlowWeight,
                 DrivenFlowSeed);
  }


  /* Now we know which hydro solver we're using, we can assign the
     default Riemann solver and flux reconstruction methods.  These
     parameters aren't used for PPM_LagrangeRemap and Zeus. */

  if (HydroMethod == PPM_DirectEuler) {
    if (RiemannSolver == INT_UNDEFINED) 
      RiemannSolver = TwoShock;
    if (ReconstructionMethod == INT_UNDEFINED)
      ReconstructionMethod = PPM;
    if (ReconstructionMethod == PLM) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("ReconstructionMethod = PLM.\n"
	       "These are the defaults for the MUSCL (hydro_rk) solvers,\n"
	       "but don't exist for the FORTRAN solvers (HydroMethod = 0).  "
	       "-- To override this, do not set ReconstructionMethod or set it to 1\n");
      ENZO_FAIL("Stopped in ReadParameterFile.\n");
      //      RiemannSolver = TwoShock;
      //      ReconstructionMethod = PPM;
    }
    if (RiemannSolver == -HLL) RiemannSolver = HLL;
  } else if (HydroMethod == HD_RK || HydroMethod == MHD_RK) {
    if (RiemannSolver == INT_UNDEFINED) 
      RiemannSolver = HLL;
    if (ReconstructionMethod == INT_UNDEFINED)
      ReconstructionMethod = PLM;
  }

  else if (HydroMethod == MHD_Li )
    if (RiemannSolver == INT_UNDEFINED) 
        RiemannSolver = HLLD;
    if (ReconstructionMethod == INT_UNDEFINED)
        ReconstructionMethod = PLM;

  if (HydroMethod==MHD_RK) UseMHD = 1;
  if (HydroMethod==MHD_Li) {UseMHDCT = 1; UseMHD = 1;}
  if (HydroMethod==MHD_Li ||HydroMethod==MHD_RK || HydroMethod==HD_RK ){
      MaxVelocityIndex = 3;
  }else{
      MaxVelocityIndex = MetaData.TopGridRank ;
  }
  if (UseMHDCT) CorrectParentBoundaryFlux = TRUE;

    if (DualEnergyFormalism == FALSE)
        MHDCTDualEnergyMethod = 0;
    else
      if ( MHDCTDualEnergyMethod == INT_UNDEFINED || MHDCTDualEnergyMethod == 0)
        MHDCTDualEnergyMethod = 2;

  //  OutputTemperature = ((ProblemType == 7) || (ProblemType == 11));

  /* Even if this is not cosmology, due to a check for nested grid cosmology
     in ProtoSubgrid_AcceptableGrid.C, we'll set the default for this here. */
  CosmologySimulationNumberOfInitialGrids = 1;

  if (HydroMethod != MHD_RK && UseMHDCT != 1)
    BAnyl = 0; // set this to zero no matter what unless we have a magnetic field to analyze.

  if ((HydroMethod != MHD_RK) && (UseGasDrag != 0))
    {
      if(MyProcessorNumber == ROOT_PROCESSOR ) {
	fprintf(stderr, "WARNING:  UseGasDrag != 0 yet HM=MHD_RK \n");
	fprintf(stderr, "WARNING:  setting UseGasDrag = 0. I.e no Gas drag included. \n");
	UseGasDrag = 0; 
      }
    }


  /* Count static nested grids since this isn't written in the
     parameter file */

  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    if (StaticRefineRegionLevel[i] != INT_UNDEFINED)
      CosmologySimulationNumberOfInitialGrids++;
 
  /* If we have turned on Comoving coordinates, read cosmology parameters. */
 
  if (ComovingCoordinates) {

    // Always output temperature in cosmology runs
    OutputTemperature = TRUE;

    if (CosmologyReadParameters(fptr, &MetaData.StopTime, &MetaData.Time)
	== FAIL) {
      ENZO_FAIL("Error in ReadCosmologyParameters.\n");
    }
    rewind(fptr);
  }
  else {
    if (ReadUnits(fptr) == FAIL){
      ENZO_FAIL("Error in ReadUnits. ");
    }
    rewind(fptr);
  }

  // make sure that MHD is turned on if we're trying to use anisotropic conduction.
  // if not, alert user.
  if(AnisotropicConduction==TRUE && UseMHD==0){
    ENZO_FAIL("AnisotropicConduction can only be used if MHD is turned on!\n");
  }  
  if(AnisotropicConduction==TRUE && MetaData.TopGridRank < 2){
    ENZO_FAIL("AnisotropicConduction can only be used if TopGridRank is >= 2!\n");
  }

  if(EquationOfState == 1 && HydroMethod != MHD_Li){
    ENZO_FAIL("If EquationOfState = 1, you must be using MHD-CT!\n");
  }


  /*
    if (EOSType == 3) // an isothermal equation of state implies the adiabatic index = 1 
    Gamma = 1; 
  */

  /* convert MustRefineParticlesMinimumMass from a mass into a density, 
     ASSUMING CUBE simulation space */
    MustRefineParticlesMinimumMass /= POW(1/(float(MetaData.TopGridDims[0])
				       *POW(float(RefineBy), float(MustRefineParticlesRefineToLevel))),3);


    /* Use Physical units stuff */

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, 
    TimeUnits = 1.0, VelocityUnits = 1.0, PressureUnits = 1.0;
  double MassUnits = 1.0;
  if (UsePhysicalUnit) {
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, 
	     &MassUnits, MetaData.Time);
    PressureUnits = DensityUnits*pow(VelocityUnits,2);
    /*IMOPORTANT: If change anything here must change both equivilant parts in WriteParameterFile.C as well */

    /* Change input physical parameters into code units */
    MustRefineParticlesMinimumMass /= MassUnits; 
    StarMakerOverDensityThreshold /= DensityUnits;
    //  StarEnergyFeedbackRate = StarEnergyFeedbackRate/pow(LengthUnits,2)*pow(TimeUnits,3);
    
    if (SinkMergeDistance > 1.0)
      SinkMergeDistance /= LengthUnits;
    //printf(" \n SinkMergeDistance = %"FSYM"\n \n", SinkMergeDistance);
    SmallRho /= DensityUnits;
    SmallP /= PressureUnits;
    SmallT /= TemperatureUnits;
    MaximumAlvenSpeed /= VelocityUnits;
    EOSSoundSpeed /=  VelocityUnits;
    float h, cs, dpdrho, dpde;
    EOS(SmallP, SmallRho, SmallEint, h, cs, dpdrho, dpde, EOSType, 1);
    if (debug && (HydroMethod == HD_RK || HydroMethod == MHD_RK))
      printf("smallrho=%g, smallp=%g, smalleint=%g, DensityUnits = %g, PressureUnits=%g, MaximumAlvenSpeed=%g\n",
	     SmallRho, SmallP, SmallEint,DensityUnits, PressureUnits, MaximumAlvenSpeed);
    for (int i = 0; i < MAX_FLAGGING_METHODS; i++) 
      if (MinimumMassForRefinement[i] != FLOAT_UNDEFINED) {
	MinimumMassForRefinement[i] /= MassUnits;
      }
    if (GravitationalConstant > 12.49 && GravitationalConstant < 12.61) {
      GravitationalConstant = 4.0 * 3.1415926 * 6.6726e-8 * DensityUnits * pow(TimeUnits,2);
      printf("Gravitational Constant recalculated from 4pi to 4piG in code units\n");
    }

  }

  /* For !restart, this only ruins the units because MinimumOverDensityForRefinement is already 
     set in SetDefaultGlobalValues and not FLOAT_UNDEFINED.
     For restart, MinimumOverDensityForRefinement is not even needs to be read because only 
     MinimumMassForRefinement is used for CellFlagging.  
     So, why did we have to do this in the first place?  - Ji-hoon Kim in Apr.2010
     (The counterpart in WriteParameterFile is also commented out) */   //#####

  /*
  if (!ComovingCoordinates && UsePhysicalUnit) 
    for (int i = 0; i < MAX_FLAGGING_METHODS; i++) 
      if (MinimumOverDensityForRefinement[i] != FLOAT_UNDEFINED) {
	MinimumOverDensityForRefinement[i] /= DensityUnits;
      }
  */
  
  /* If RefineRegionTimeType is 0 or 1, read in the input file. */
  if ((RefineRegionTimeType == 0) || (RefineRegionTimeType == 1)) {
      if (ReadEvolveRefineFile() == FAIL) {
        ENZO_FAIL("Error in ReadEvolveRefineFile.");
      }
  }

#ifdef USE_GRACKLE
  /* If using Grackle chemistry and cooling library, override all other 
     cooling machinery and do a translation of some of the parameters. */
  if (grackle_data.use_grackle == TRUE) {
    // grackle_data.use_grackle already set
    // grackle_data.with_radiative_cooling already set
    // grackle_data.grackle_data_file already set
    // grackle_data.UVbackground already set
    // grackle_data.Compton_xray_heating already set
    // grackle_data.LWbackground_intensity already set
    // grackle_data.LWbackground_sawtooth_suppression already set
    grackle_data.Gamma                          = (double) Gamma;
    grackle_data.primordial_chemistry           = (Eint32) MultiSpecies;
    grackle_data.metal_cooling                  = (Eint32) MetalCooling;
    grackle_data.h2_on_dust                     = (Eint32) H2FormationOnDust;
    grackle_data.cmb_temperature_floor          = (Eint32) CloudyCoolingData.CMBTemperatureFloor;
    grackle_data.three_body_rate                = (Eint32) ThreeBodyRate;
    grackle_data.cie_cooling                    = (Eint32) CIECooling;
    grackle_data.h2_optical_depth_approximation = (Eint32) H2OpticalDepthApproximation;
    grackle_data.photoelectric_heating          = (Eint32) PhotoelectricHeating;
    grackle_data.photoelectric_heating_rate     = (double) PhotoelectricHeatingRate;
    grackle_data.NumberOfTemperatureBins        = (Eint32) CoolData.NumberOfTemperatureBins;
    grackle_data.CaseBRecombination             = (Eint32) RateData.CaseBRecombination;
    grackle_data.TemperatureStart               = (double) CoolData.TemperatureStart;
    grackle_data.TemperatureEnd                 = (double) CoolData.TemperatureEnd;
    grackle_data.NumberOfDustTemperatureBins    = (Eint32) RateData.NumberOfDustTemperatureBins;
    grackle_data.DustTemperatureStart           = (double) RateData.DustTemperatureStart;
    grackle_data.DustTemperatureEnd             = (double) RateData.DustTemperatureEnd;
    grackle_data.HydrogenFractionByMass         = (double) CoolData.HydrogenFractionByMass;
    grackle_data.DeuteriumToHydrogenRatio       = (double) CoolData.DeuteriumToHydrogenRatio;
    grackle_data.SolarMetalFractionByMass       = (double) CoolData.SolarMetalFractionByMass;

    // Initialize units structure.
    FLOAT a_value, dadt;
    a_value = 1.0;
    code_units grackle_units;
    grackle_units.a_units = 1.0;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                 &TimeUnits, &VelocityUnits, MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in GetUnits.\n");
    }
    if (ComovingCoordinates) {
      if (CosmologyComputeExpansionFactor(MetaData.Time, &a_value, 
                                          &dadt) == FAIL) {
        ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
      }
      grackle_units.a_units            = (double) (1.0 / (1.0 + InitialRedshift));
    }
    grackle_units.comoving_coordinates = (Eint32) ComovingCoordinates;
    grackle_units.density_units        = (double) DensityUnits;
    grackle_units.length_units         = (double) LengthUnits;
    grackle_units.time_units           = (double) TimeUnits;
    grackle_units.velocity_units       = (double) VelocityUnits;

    // Initialize chemistry structure.
    if (initialize_chemistry_data(&grackle_units,
                                  (double) a_value) == FAIL) {
      ENZO_FAIL("Error in Grackle initialize_chemistry_data.\n");
    }
  }  // if (grackle_data.use_grackle == TRUE)

  else {
#endif // USE_GRACKE

    /* If GadgetEquilibriumCooling == TRUE, we don't want MultiSpecies
       or RadiationFieldType to be on - both are taken care of in
       the Gadget cooling routine.  Therefore, we turn them off!
       Also, initialize the Gadget equilibrium cooling data. */

    if(GadgetEquilibriumCooling == TRUE){

      if(MyProcessorNumber == ROOT_PROCESSOR ) {
        fprintf(stderr, "WARNING:  GadgetEquilibriumCooling = 1.  Forcing\n");
        fprintf(stderr, "WARNING:  RadiationFieldType = 0, MultiSpecies = 0, and\n");
        fprintf(stderr, "WARNING:  RadiativeCooling = 1.\n");
      }

      RadiationFieldType = 0;
      MultiSpecies       = 0;
      RadiativeCooling   = 1;

      // initialize Gadget equilibrium cooling
      if (InitializeGadgetEquilibriumCoolData(MetaData.Time) == FAIL) {
        ENZO_FAIL("Error in InitializeGadgetEquilibriumCoolData.");
      } 
    }

    /* If set, initialize the RadiativeCooling and RateEquations data. */

    if (MultiSpecies > 0) {
      if (InitializeRateData(MetaData.Time) == FAIL) {
        ENZO_FAIL("Error in InitializeRateData.");
      }
    }
 
    if (MultiSpecies             == 0 && 
        MetalCooling             == 0 &&
        GadgetEquilibriumCooling == 0 &&
        RadiativeCooling          > 0) {
      if (InitializeEquilibriumCoolData(MetaData.Time) == FAIL) {
        ENZO_FAIL("Error in InitializeEquilibriumCoolData.");
      }
    }

    /* If using the internal radiation field, initialize it. */
 
    if (RadiationFieldType == 11) 
      RadiationData.RadiationShield = TRUE; 
    else if (RadiationFieldType == 10)
      RadiationData.RadiationShield = FALSE; 

    if ((RadiationFieldType >= 10 && RadiationFieldType <= 11) ||
        RadiationData.RadiationShield == TRUE)
      if (InitializeRadiationFieldData(MetaData.Time) == FAIL) {
	ENZO_FAIL("Error in InitializeRadiationFieldData.");
      }
 
#ifdef USE_GRACKLE
  } // else (if Grackle == TRUE)
#endif

  /* If using MBHFeedback = 2 to 5 (Star->FeedbackFlag = MBH_JETS), 
     you need MBHParticleIO for angular momentum */

  if (MBHFeedback >= 2 && MBHFeedback <= 5) 
    MBHParticleIO = TRUE;

  /* Turn off DualEnergyFormalism for zeus hydro (and a few other things). */
 
  if (HydroMethod == Zeus_Hydro) {
    ConservativeInterpolation = FALSE;
    DualEnergyFormalism       = FALSE;
    //    FluxCorrection            = FALSE;
  }

  if (DualEnergyFormalism > 0 && EOSType > 0)
    ENZO_FAIL("DualEnergyFormalism should be off for EOSType > 0");

  /* Set some star feedback parameters. */

  if ((STARFEED_METHOD(NORMAL_STAR) || STARFEED_METHOD(UNIGRID_STAR)) && 
      (StarFeedbackDistRadius > 0)) {

    // Calculate number of cells in the shape over which to distribute feedback.
    StarFeedbackDistRadius = min(StarFeedbackDistRadius,
				 StarFeedbackDistCellStep);
    int i, j, k, cell_step;

    StarFeedbackDistTotalCells = 0;
    for (k = -StarFeedbackDistRadius;k <= StarFeedbackDistRadius;k++) {
      for (j = -StarFeedbackDistRadius;j <= StarFeedbackDistRadius;j++) {
	for (i = -StarFeedbackDistRadius;i <= StarFeedbackDistRadius;i++) {
	  cell_step = fabs(k) + fabs(j) + fabs(i);
	  if (cell_step <= StarFeedbackDistCellStep) {
	    StarFeedbackDistTotalCells++;
	  }
       }
      }
    }
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr,"Total cells for star feedback smoothing: %"ISYM".\n",
	      StarFeedbackDistTotalCells);
    }
  }
 
  /* For rk_hydro, we need to set some variables */

  if (DualEnergyFormalism) {
    NEQ_HYDRO = 6;
    NEQ_MHD   = 10;
    ieint = 5;
    iBx = 6;
    iBy = 7;
    iBz = 8;
    iPhi = 9;
    iEint = 5;
  }

  // Don't include free electron field
  switch (MultiSpecies) {
  case 0:  NSpecies = 0;  break;
  case 1:  NSpecies = 5;  break;
  case 2:  NSpecies = 8;  break;
  case 3:  NSpecies = 11; break;
  default: NSpecies = 0;  break;
  }

  // Determine color fields (NColor) later inside a grid object.
  // ...
#ifdef UNUSED
  if (MaximumGravityRefinementLevel == INT_UNDEFINED)
    MaximumGravityRefinementLevel = (RadiativeCooling && SelfGravity
				     && HydroMethod == Zeus_Hydro) ?
       max(MaximumRefinementLevel-2, 5) : MaximumRefinementLevel;
#else
  if (MaximumGravityRefinementLevel == INT_UNDEFINED)
    MaximumGravityRefinementLevel = MaximumRefinementLevel;
#endif

    ret += sscanf(line, "IsothermalSoundSpeed = %"GSYM, &IsothermalSoundSpeed);
    ret += sscanf(line, "RefineByJeansLengthUnits = %"ISYM, &RefineByJeansLengthUnits);
 
  MaximumGravityRefinementLevel =
    min(MaximumGravityRefinementLevel, MaximumRefinementLevel);
 
  /* If MultiSpecies < 2, we can't simulate Pop III star formation */

  if (MultiSpecies < 2 && STARMAKE_METHOD(POP3_STAR)) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Cannot form Pop III stars without H2 cooling!\n"
	      "Turning Pop III star formation OFF.\n");
    StarParticleCreation -= 1 << POP3_STAR;
  }

  /* Use the value in MaximumParticleRefinementLevel to set the smoothing
     radius for the particles, to be used to Grid_DepositPositions. */
 
  if (MaximumParticleRefinementLevel >= 0)
    DepositPositionsParticleSmoothRadius =
      (DomainRightEdge[0] - DomainLeftEdge[0])/
      (float(MetaData.TopGridDims[0])*
       POW(float(RefineBy), float(MaximumParticleRefinementLevel)));
  else
    DepositPositionsParticleSmoothRadius = 0;
 
  

//  PPMDiffusion causes an out-of-bounds condition as currently written
//  The following is an over-ride to force PPMDiffusion OFF. This has
//  been fixed in this latest version (AK).

//  NOTE: The fix keeps the code from crashing, but is not a proper 
//  implementation of PPM diffusion.  The reason why is that Enzo typically
//  uses 3 ghost zones, and the correct PPM diffusion implementation requires
//  4 parameters.  SO, you should not use this parameter for, e.g., cosmology
//  runs unless you know what you're doing.  (BWO)

  if (MetaData.PPMDiffusionParameter != 0 && ProblemType != 60 // Turbulence
                                          && ProblemType != 4  // Double Mach Reflection test
                                          && ProblemType != 6  // Implosion test
                                          && ProblemType != 7  // SedovBlast test
                                          && ProblemType != 8  // KH test
                                          && ProblemType != 9  // Noh test
                                          && ProblemType != 11 // Radiating shock test
                                          ) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("WARNING! Setting MetaData.PPMDiffusionParameter = 0\n");
    MetaData.PPMDiffusionParameter = 0;
  }
 
  if (PartitionNestedGrids == 1) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("WARNING! PartitionNestedGrids = 1 forces Parallel IO = 1\n");
    ParallelRootGridIO = 1;
    ParallelParticleIO = 1;
  }

  if(ProblemType==70 && UseHydro==1){
    printf("ReadParameterFile: ProblemType=70.  Disabling hydrodynamics!\n");
    UseHydro=FALSE;
  }



  if ((MetaData.GravityBoundary != TopGridPeriodic) &&
      (UnigridTranspose)) {
    /* it turns out that Robert Harkness' unigrid transpose stuff is incompatible with the top
       grid isolated gravity boundary conditions.  I'm not 100 percent sure why this is - in the 
       meantime, just double-check to make sure that if one tries to use the isolated boundary
       conditions when the unigrid transpose stuff is on, the code crashes loudly.
       -- BWO, 26 June 2008 */
    ENZO_FAIL("Parameter mismatch: TopGridGravityBoundary = 1 only works with UnigridTranspose = 0");
  }

  /* If the restart dump parameters were set to the previous defaults
     (dtRestartDump = 5 hours), then set back to current default,
     which is no restart dumps. */

  float tol = 1e-6;
  if (ABS(MetaData.dtRestartDump - 3600.0*5) / (3600.0*5) < tol &&
      ABS(MetaData.TimeLastRestartDump) < tol) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, 
	      "==================================================\n"
	      "-> Turning off restart dumps because the previous\n"
	      "default was set to 5 hours but not used.  To avoid this warning,\n"
	      "set dtRestartDump to a negative value.  If you wanted to use\n"
	      "restart dumps, please set it != 18000.\n"
	      "==================================================\n");
    MetaData.dtRestartDump = FLOAT_UNDEFINED;
  }

  /* If refining by must-refine particles, particle mass refinement
     must be turned on. */

  int method;
  bool TurnOnParticleMassRefinement = false;
  for (method = 0; method < MAX_FLAGGING_METHODS; method++) 
    if (CellFlaggingMethod[method] == 8) {
      TurnOnParticleMassRefinement = true;
      break;
    }
  for (method = 0; method < MAX_FLAGGING_METHODS; method++) 
    if (CellFlaggingMethod[method] == 4) {
      TurnOnParticleMassRefinement = false;
      break;
    }
  if (TurnOnParticleMassRefinement) {
    method = 0;
    while (CellFlaggingMethod[method] != INT_UNDEFINED)
      method++;
    CellFlaggingMethod[method] = 4;
  }


  /* If we're refining the region around P3 supernovae,
     MustRefineByParticles must be set.  Check this.  */

  if (PopIIISupernovaMustRefine == TRUE) {
    bool TurnOnParticleMustRefine = true;
    for (method = 0; method < MAX_FLAGGING_METHODS; method++)
      if (CellFlaggingMethod[method] == 8)
	TurnOnParticleMustRefine = false;
    if (TurnOnParticleMustRefine) {
      method = 0;
      while (CellFlaggingMethod[method] != INT_UNDEFINED)
	method++;
      CellFlaggingMethod[method] = 8;
    }

    /* Check if the must refine level is still at the default.  If so,
       break because it's zero!  Won't do anything, and the user will
       be disappointed to find that the simulation didn't refine
       around the SN. */

    if (MustRefineParticlesRefineToLevel == 0)
      ENZO_FAIL("MustRefineParticlesRefineToLevel is still ZERO, and you set"
		"PopIIISupernovaMustRefine.  Set the level or turn off"
		"PopIIISupernovaMustRefine.");
  } // ENDIF PopIIISupernovaMustRefine
//del

  if (TracerParticleOn) {
    ParticleTypeInFile = TRUE;
  }

  if (OutputParticleTypeGrouping && (!ParticleTypeInFile)) {
    OutputParticleTypeGrouping = FALSE;
  }
 
  //  if (WritePotential && ComovingCoordinates && SelfGravity) {
  if (WritePotential && SelfGravity) {
    CopyGravPotential = TRUE;
  }

  // Keep track of number of outputs left until exit.
  MetaData.OutputsLeftBeforeExit = MetaData.NumberOfOutputsBeforeExit;

  if (MyProcessorNumber == ROOT_PROCESSOR) {
 
    if ( MetaData.GlobalDir != NULL ) {
      fprintf(stderr, "Output to Global Dir %s\n", MetaData.GlobalDir);
    }
 
    if ( MetaData.LocalDir != NULL ) {
      fprintf(stderr, "Output to Local Dir %s\n", MetaData.LocalDir);
    }

  }
 
  if ( (MetaData.GlobalDir != NULL) && (MetaData.LocalDir != NULL) ) {
    ENZO_FAIL("Cannot select GlobalDir AND LocalDir!\n");
  }
 
  char *cwd_buffer = new char[MAX_LINE_LENGTH];
  size_t cwd_buffer_len = MAX_LINE_LENGTH;
 
  if ( (MetaData.GlobalDir == NULL) && (MetaData.LocalDir == NULL) ) {
    /*if(getcwd(cwd_buffer, cwd_buffer_len) == NULL) {
      fprintf(stderr, "GETCWD call FAILED\n");
    }
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,"CWD %s\n", cwd_buffer);
    */
    /* No one seems to want GlobalDir to default to abspath(CWD).  I'm leaving
       the code here in case you do. MJT */ 
    strcpy(cwd_buffer, ".");
    MetaData.GlobalDir = cwd_buffer;
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,"Global Dir set to %s\n", cwd_buffer);
  }

  /* Generate unique identifier if one wasn't found. */
  if(MetaData.SimulationUUID == NULL){
    MetaData.SimulationUUID = new char[MAX_LINE_LENGTH];
    get_uuid(MetaData.SimulationUUID);
  }
 
   for (int i=0; i<MetaData.TopGridRank;i++)
    TopGridDx[i]=(DomainRightEdge[i]-DomainLeftEdge[i])/MetaData.TopGridDims[i];

 //  for (int i=0; i<MetaData.TopGridRank; i++)
//      fprintf (stderr, "read  %"ISYM"  %"ISYM" \n", 
// 	      MetaData.LeftFaceBoundaryCondition[i], 
// 	      MetaData.RightFaceBoundaryCondition[i]);

  if (UseCUDA) {
    LoadBalancing = 0; // Should explore how LoadBalancing = 1 gives problems with CUDA
#ifndef ECUDA
    printf("This executable was compiled without CUDA support.\n");
    printf("use \n");
    printf("make cuda-yes\n");
    printf("Exiting.\n");
    my_exit(EXIT_SUCCESS);
#endif
  }


  if (debug) printf("Initialdt in ReadParameterFile = %e\n", *Initialdt);

  //

  CheckShearingBoundaryConsistency(MetaData);

  return SUCCESS;
#endif /* ndef CONFIG_USE_LIBCONFIG */
}
