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
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "StarParticleData.h"
#include "hydro_rk/EOS.h" 
 
/* This variable is declared here and only used in Grid_ReadGrid. */
 


/* function prototypes */
 
int ReadListOfFloats(FILE *fptr, int N, float floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
int CosmologyReadParameters(FILE *fptr, FLOAT *StopTime, FLOAT *InitTime);
int ReadUnits(FILE *fptr);
int InitializeCloudyCooling(FLOAT Time);
int InitializeRateData(FLOAT Time);
int InitializeEquilibriumCoolData(FLOAT Time);
int InitializeGadgetEquilibriumCoolData(FLOAT Time);
int InitializeRadiationFieldData(FLOAT Time);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
 
 
int CheckShearingBoundaryConsistency(TopGridData &MetaData); 

int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt)
{
  /* declarations */

  
  char line[MAX_LINE_LENGTH];
  int dim, ret, int_dummy;
  float TempFloat;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  int comment_count = 0;
 
  /* read until out of lines */
 
  while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL) 
      && (comment_count < 2)) {

    ret = 0;
 
    /* read MetaData parameters */
 
    ret += sscanf(line, "InitialCycleNumber = %"ISYM, &MetaData.CycleNumber);
    ret += sscanf(line, "InitialTime        = %"PSYM, &MetaData.Time);
    ret += sscanf(line, "InitialCPUTime     = %lf", &MetaData.CPUTime);
    ret += sscanf(line, "Initialdt          = %"FSYM, &Initialdt);
 
    ret += sscanf(line, "StopTime    = %"PSYM, &MetaData.StopTime);
    ret += sscanf(line, "StopCycle   = %"ISYM, &MetaData.StopCycle);
    ret += sscanf(line, "StopSteps   = %"ISYM, &MetaData.StopSteps);
    ret += sscanf(line, "StopCPUTime = %"FSYM, &MetaData.StopCPUTime);
    ret += sscanf(line, "ResubmitOn  = %"ISYM, &MetaData.ResubmitOn);
    if (sscanf(line, "ResubmitCommand = %s", dummy) == 1) 
      MetaData.ResubmitCommand = dummy;
 
    ret += sscanf(line, "TimeLastRestartDump = %"PSYM,
		  &MetaData.TimeLastRestartDump);
    ret += sscanf(line, "dtRestartDump       = %"PSYM, &MetaData.dtRestartDump);
    ret += sscanf(line, "TimeLastDataDump    = %"PSYM,
		  &MetaData.TimeLastDataDump);
    ret += sscanf(line, "dtDataDump          = %"PSYM, &MetaData.dtDataDump);
    ret += sscanf(line, "TimeLastHistoryDump = %"PSYM,
		  &MetaData.TimeLastHistoryDump);
    ret += sscanf(line, "dtHistoryDump       = %"PSYM, &MetaData.dtHistoryDump);
 
    ret += sscanf(line, "TracerParticleOn  = %"ISYM, &TracerParticleOn);
    ret += sscanf(line, "ParticleTypeInFile = %"ISYM, &ParticleTypeInFile);
    ret += sscanf(line, "TimeLastTracerParticleDump = %"PSYM,
                  &MetaData.TimeLastTracerParticleDump);
    ret += sscanf(line, "dtTracerParticleDump       = %"PSYM,
                  &MetaData.dtTracerParticleDump);
 
    ret += sscanf(line, "NewMovieLeftEdge  = %"FSYM" %"FSYM" %"FSYM, 
		  MetaData.NewMovieLeftEdge,
		  MetaData.NewMovieLeftEdge+1, 
		  MetaData.NewMovieLeftEdge+2);
    ret += sscanf(line, "NewMovieRightEdge = %"FSYM" %"FSYM" %"FSYM, 
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
 
    /* Subcycle directed output */
    ret += sscanf(line, "SubcycleSkipDataDump = %"ISYM, 
                  &MetaData.SubcycleSkipDataDump);
    ret += sscanf(line, "SubcycleLastDataDump = %"ISYM, 
                  &MetaData.SubcycleLastDataDump);
    ret += sscanf(line, "SubcycleNumber = %"ISYM, 
                  &MetaData.SubcycleNumber);

    ret += sscanf(line,"FileDirectedOutput = %"ISYM,
		  &FileDirectedOutput);

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
        fprintf(stderr, "CubeDump %"ISYM" > maximum allowed.\n", dim);
        ENZO_FAIL("");
      }
    }

    ret += sscanf(line, "LoadBalancing = %"ISYM, &LoadBalancing);
 
    if (sscanf(line, "TimeActionType[%"ISYM"] = %"ISYM, &dim, &int_dummy) == 2) {
      ret++; TimeActionType[dim] = int_dummy;
      if (dim >= MAX_TIME_ACTIONS-1) {
	fprintf(stderr, "Time action %"ISYM" > maximum allowed.\n", dim);
	ENZO_FAIL("");
      }
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
 
    ret += sscanf(line, "ParticleBoundaryType   = %"ISYM,
		  &MetaData.ParticleBoundaryType);
    ret += sscanf(line, "NumberOfParticles      = %"ISYM,
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
    ret += sscanf(line, "HydroMethod            = %"ISYM, &HydroMethod);
    if (HydroMethod==MHD_RK) useMHD = 1;

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
    ret += sscanf(line, "InterpolationMethod    = %"ISYM, &InterpolationMethod);
    ret += sscanf(line, "ConservativeInterpolation = %"ISYM,
		  &ConservativeInterpolation);
    ret += sscanf(line, "MinimumEfficiency      = %"FSYM, &MinimumEfficiency);
    ret += sscanf(line, "MinimumSubgridEdge     = %"ISYM, &MinimumSubgridEdge);
    ret += sscanf(line, "MaximumSubgridSize     = %"ISYM, &MaximumSubgridSize);
    ret += sscanf(line, "NumberOfBufferZones    = %"ISYM, &NumberOfBufferZones);
    ret += sscanf(line, "MustRefineRegionMinRefinementLevel = %"ISYM,
		  &MustRefineRegionMinRefinementLevel);
    ret += sscanf(line, "MetallicityRefinementMinLevel = %"ISYM,
		  &MetallicityRefinementMinLevel);
    ret += sscanf(line, "MetallicityRefinementMinMetallicity      = %"FSYM, 
		  &MetallicityRefinementMinMetallicity);

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
 
    ret += sscanf(line, "SelfGravity           = %"ISYM, &SelfGravity);
    ret += sscanf(line, "GravitationalConstant = %"FSYM, &GravitationalConstant);
    ret += sscanf(line, "S2ParticleSize        = %"FSYM, &S2ParticleSize);
    ret += sscanf(line, "GravityResolution     = %"FSYM, &GravityResolution);
    ret += sscanf(line, "ComputePotential      = %"ISYM, &ComputePotential);
    ret += sscanf(line, "PotentialIterations   = %"ISYM, &PotentialIterations);
    ret += sscanf(line, "WritePotential        = %"ISYM, &WritePotential);
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
    ret += sscanf(line, "RadiativeCooling = %"ISYM, &RadiativeCooling);
    ret += sscanf(line, "GadgetEquilibriumCooling = %"ISYM, &GadgetEquilibriumCooling);
    ret += sscanf(line, "MultiSpecies = %"ISYM, &MultiSpecies);
    if (sscanf(line, "CloudyCoolingGridFile = %s", dummy) == 1) {
      CloudyCoolingData.CloudyCoolingGridFile = dummy;
      ret++;
    }
    ret += sscanf(line, "IncludeCloudyHeating = %"ISYM, &CloudyCoolingData.IncludeCloudyHeating);
    ret += sscanf(line, "IncludeCloudyMMW = %"ISYM, &CloudyCoolingData.IncludeCloudyMMW);
    ret += sscanf(line, "CMBTemperatureFloor = %"ISYM, &CloudyCoolingData.CMBTemperatureFloor);
    ret += sscanf(line, "ConstantTemperatureFloor = %"FSYM, &CloudyCoolingData.ConstantTemperatureFloor);
    ret += sscanf(line, "CloudyMetallicityNormalization = %"FSYM,&CloudyCoolingData.CloudyMetallicityNormalization);
    ret += sscanf(line, "CloudyElectronFractionFactor = %"FSYM,&CloudyCoolingData.CloudyElectronFractionFactor);
    ret += sscanf(line, "MetalCooling = %d", &MetalCooling);
    if (sscanf(line, "MetalCoolingTable = %s", dummy) == 1) 
      MetalCoolingTable = dummy;
    ret += sscanf(line, "RadiationFieldType = %"ISYM, &RadiationFieldType);
    ret += sscanf(line, "AdjustUVBackground = %"ISYM, &AdjustUVBackground);
    ret += sscanf(line, "SetUVBAmplitude = %"FSYM, &SetUVBAmplitude);
    ret += sscanf(line, "SetHeIIHeatingScale = %"FSYM, &SetHeIIHeatingScale);

    ret += sscanf(line, "RadiationFieldLevelRecompute = %"ISYM,
		  &RadiationFieldLevelRecompute);
    ret += sscanf(line, "RadiationSpectrumNormalization = %"FSYM,
		  &CoolData.f3);
    ret += sscanf(line, "RadiationSpectrumSlope = %"FSYM, &CoolData.alpha0);

    if (sscanf(line, "CoolDataParameterFile = %s", dummy) == 1)
      CoolData.ParameterFilename = dummy;

    ret += sscanf(line, "OutputCoolingTime = %"ISYM, &OutputCoolingTime);
    ret += sscanf(line, "OutputTemperature = %"ISYM, &OutputTemperature);

    ret += sscanf(line, "OutputSmoothedDarkMatter = %"ISYM, 
		  &OutputSmoothedDarkMatter);
    ret += sscanf(line, "SmoothedDarkMatterNeighbors = %"ISYM, 
		  &SmoothedDarkMatterNeighbors);

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
    ret += sscanf(line, "RefineByResistiveLengthSafetyFactor = %" FSYM,
		  &RefineByResistiveLengthSafetyFactor);
    ret += sscanf(line, "MustRefineParticlesRefineToLevel = %"ISYM,
                  &MustRefineParticlesRefineToLevel);
    ret += sscanf(line, "ParticleTypeInFile = %"ISYM,
                  &ParticleTypeInFile);

 
    if (sscanf(line, "StaticRefineRegionLevel[%"ISYM"] = %"ISYM,&dim,&int_dummy) == 2){
      if (dim > MAX_STATIC_REGIONS-1) {
        fprintf(stderr, "StaticRegion number %"ISYM" > MAX allowed\n", dim);
        ENZO_FAIL("");
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
 
    ret += sscanf(line, "PartitionNestedGrids = %"ISYM, &PartitionNestedGrids);
 
    ret += sscanf(line, "ExtractFieldsOnly = %"ISYM, &ExtractFieldsOnly);
 
    ret += sscanf(line, "CubeDumpEnabled = %"ISYM, &CubeDumpEnabled);
 
    ret += sscanf(line, "Debug1 = %"ISYM, &debug1);

    ret += sscanf(line, "Debug2 = %"ISYM, &debug2);

    ret += sscanf(line, "MemoryLimit = %"ISYM, &MemoryLimit);

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
    ret += sscanf(line, "MinimumSlopeForRefinement ="
		  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
		  MinimumSlopeForRefinement+0,
		  MinimumSlopeForRefinement+1,
		  MinimumSlopeForRefinement+2,
		  MinimumSlopeForRefinement+3,
		  MinimumSlopeForRefinement+4,
		  MinimumSlopeForRefinement+5,
		  MinimumSlopeForRefinement+6);
    ret += sscanf(line, "MinimumPressureJumpForRefinement = %"FSYM,
		  &MinimumPressureJumpForRefinement);
    ret += sscanf(line, "MinimumShearForRefinement = %"FSYM,
		  &MinimumShearForRefinement);
    ret += sscanf(line, "MinimumEnergyRatioForRefinement = %"FSYM,
		  &MinimumEnergyRatioForRefinement);
    ret += sscanf(line, "ComovingCoordinates = %"ISYM,&ComovingCoordinates);
    ret += sscanf(line, "StarParticleCreation = %"ISYM, &StarParticleCreation);
    ret += sscanf(line, "StarParticleFeedback = %"ISYM, &StarParticleFeedback);
    ret += sscanf(line, "NumberOfParticleAttributes = %"ISYM,
		  &NumberOfParticleAttributes);
    ret += sscanf(line, "AddParticleAttributes = %"ISYM, &AddParticleAttributes);

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
 
    /* Check version number. */
 
    if (sscanf(line, "VersionNumber = %"FSYM, &TempFloat) == 1) {
      ret++;
      if (fabs(TempFloat - VERSION) >= 1.0e-3)
	fprintf(stderr, "Warning: Incorrect version number.\n");
    }
 
    /* Read star particle parameters. */
 
    ret += sscanf(line, "StarMakerOverDensityThreshold = %"FSYM,
		  &StarMakerOverDensityThreshold);
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
 

    ret += sscanf(line, "StarClusterUseMetalField = %"ISYM, 
		  &StarClusterUseMetalField);
    ret += sscanf(line, "StarClusterMinDynamicalTime = %"FSYM, 
		  &StarClusterMinDynamicalTime);
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

    ret += sscanf(line, "MBHUseMetalField = %"ISYM, 
		  &MBHUseMetalField);
    ret += sscanf(line, "MBHMinDynamicalTime = %"FSYM, 
		  &MBHMinDynamicalTime);
    ret += sscanf(line, "MBHFeedbackEnergy = %lf", &MBHFeedbackEnergy);
    ret += sscanf(line, "MBHFeedbackRadius = %"FSYM, &MBHFeedbackRadius);
    ret += sscanf(line, "MBHMinimumMass = %"FSYM, 
		  &MBHMinimumMass);
    ret += sscanf(line, "MBHCombineRadius = %"FSYM,
		  &MBHCombineRadius);
    ret += sscanf(line, "MBHIonizingLuminosity = %lf", 
		  &MBHIonizingLuminosity);

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
    ret += sscanf(line, "MovieTimestepCounter = %"ISYM, &MetaData.TimestepCounter);

    ret += sscanf(line, "MultiMetals = %"ISYM, &MultiMetals);

    ret += sscanf(line, "RadiativeTransfer = %"ISYM, &RadiativeTransfer);
    ret += sscanf(line, "RadiationXRaySecondaryIon = %"ISYM, 
		  &RadiationXRaySecondaryIon);


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
    ret += sscanf(line, "PythonSubcycleSkip = %"ISYM, &PythonSubcycleSkip);

    /* Inline halo finder */

    ret += sscanf(line, "InlineHaloFinder = %"ISYM, &InlineHaloFinder);
    ret += sscanf(line, "HaloFinderSubfind = %"ISYM, &HaloFinderSubfind);
    ret += sscanf(line, "HaloFinderOutputParticleList = %"ISYM, 
		  &HaloFinderOutputParticleList);
    ret += sscanf(line, "HaloFinderLinkingLength = %"FSYM, 
		  &HaloFinderLinkingLength);
    ret += sscanf(line, "HaloFinderMinimumSize = %"ISYM, &HaloFinderMinimumSize);
    ret += sscanf(line, "HaloFinderCycleSkip = %"ISYM, &HaloFinderCycleSkip);
    ret += sscanf(line, "HaloFinderTimestep = %"FSYM, &HaloFinderTimestep);
    ret += sscanf(line, "HaloFinderLastTime = %"PSYM, &HaloFinderLastTime);

    /* This Block for Stanford Hydro */

    ret += sscanf(line, "UseHydro               = %"ISYM, &UseHydro);


    /* Sink particles (for present day star formation) & winds */
    ret += sscanf(line, "SinkMergeDistance = %lf", &SinkMergeDistance);
    ret += sscanf(line, "SinkMergeMass        = %"FSYM, &SinkMergeMass);
    ret += sscanf(line, "StellarWindFeedback  = %"ISYM, &StellarWindFeedback);
    ret += sscanf(line, "StellarWindTurnOnMass = %"FSYM, &StellarWindTurnOnMass);

    //    ret += sscanf(line, "VelAnyl = %"ISYM, &VelAnyl);


    /* Read MHD Paramters */
    ret += sscanf(line, "UseDivergenceCleaning = %d", &UseDivergenceCleaning);
    ret += sscanf(line, "DivergenceCleaningBoundaryBuffer = %d", &DivergenceCleaningBoundaryBuffer);
    ret += sscanf(line, "DivergenceCleaningThreshold = %f", &DivergenceCleaningThreshold);
    ret += sscanf(line, "PoissonApproximationThreshold = %f", &PoissonApproximationThreshold);
    ret += sscanf(line, "AngularVelocity = %f", &AngularVelocity);
    ret += sscanf(line, "VelocityGradient = %f", &VelocityGradient);
    ret += sscanf(line, "UseDrivingField = %d", &UseDrivingField);
    ret += sscanf(line, "DrivingEfficiency = %f", &DrivingEfficiency);

    ret += sscanf(line, "StringKick = %d", &StringKick);
    ret += sscanf(line, "UsePhysicalUnit = %d", &UsePhysicalUnit);
    ret += sscanf(line, "Theta_Limiter = %f", &Theta_Limiter);
    ret += sscanf(line, "RKOrder = %d", &RKOrder);
    ret += sscanf(line, "UseFloor = %d", &UseFloor);
    ret += sscanf(line, "UseViscosity = %d", &UseViscosity);
    ret += sscanf(line, "UseAmbipolarDiffusion = %d", &UseAmbipolarDiffusion);
    ret += sscanf(line, "UseResistivity = %d", &UseResistivity);
    ret += sscanf(line, "SmallRho = %g", &SmallRho);
    ret += sscanf(line, "SmallP = %g", &SmallP);
    ret += sscanf(line, "SmallT = %g", &SmallT);
    ret += sscanf(line, "MaximumAlvenSpeed = %g", &MaximumAlvenSpeed);
    ret += sscanf(line, "Coordinate = %"ISYM, &Coordinate);
    ret += sscanf(line, "RiemannSolver = %"ISYM, &RiemannSolver);
    ret += sscanf(line, "ReconstructionMethod = %"ISYM, &ReconstructionMethod);
    ret += sscanf(line, "EOSType = %"ISYM, &EOSType);
    ret += sscanf(line, "EOSSoundSpeed = %"FSYM, &EOSSoundSpeed);
    ret += sscanf(line, "EOSCriticalDensity = %"FSYM, &EOSCriticalDensity);
    ret += sscanf(line, "EOSGamma = %"FSYM, &EOSGamma);
    ret += sscanf(line, "UseConstantAcceleration = %"ISYM, &UseConstantAcceleration);
    ret += sscanf(line, "ConstantAcceleration = %"GSYM" %"GSYM" %"GSYM, &ConstantAcceleration[0],
		  &ConstantAcceleration[1], &ConstantAcceleration[2]);
    ret += sscanf(line, "Mu = %"GSYM, &Mu);
    ret += sscanf(line, "CoolingCutOffDensity1 = %"GSYM, &CoolingCutOffDensity1);
    ret += sscanf(line, "CoolingCutOffDensity2 = %"GSYM, &CoolingCutOffDensity2);
    ret += sscanf(line, "CoolingCutOffTemperature = %"GSYM, &CoolingCutOffTemperature);
    ret += sscanf(line, "CoolingPowerCutOffDensity1 = %"GSYM, &CoolingPowerCutOffDensity1);
    ret += sscanf(line, "CoolingPowerCutOffDensity2 = %"GSYM, &CoolingPowerCutOffDensity2);
    ret += sscanf(line, "UseH2OnDust           = %"ISYM, &UseH2OnDust);
    ret += sscanf(line, "PhotoelectricHeating  = %lf", &PhotoelectricHeating);

#ifdef ECUDA
    ret += sscanf(line, "UseCUDA = %"ISYM,&UseCUDA);
#endif

    /* If the dummy char space was used, then make another. */
 
    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
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
    if (strstr(line, "TestOrbit")    ) ret++;
    if (strstr(line, "KelvinHelmholtz")     ) ret++;
    if (strstr(line, "KH")                  ) ret++;
    if (strstr(line, "Noh")                 ) ret++;
    if (strstr(line, "ZeldovichPancake")    ) ret++;
    if (strstr(line, "PressurelessCollapse")) ret++;
    if (strstr(line, "AdiabaticExpansion" ) ) ret++;
    if (strstr(line, "CosmologySimulation") ) ret++;
    if (strstr(line, "TestGravity"        ) ) ret++;
    if (strstr(line, "SphericalInfall"    ) ) ret++;
    if (strstr(line, "TestGravitySphere"  ) ) ret++;
    if (strstr(line, "CollapseTest"       ) ) ret++;
    if (strstr(line, "Cosmology"          ) ) ret++;
    if (strstr(line, "SupernovaRestart"   ) ) ret++;
    if (strstr(line, "TracerParticleCreation")) ret++;
    if (strstr(line, "TurbulenceSimulation")) ret++;
    if (strstr(line, "ProtostellarCollapse")) ret++;
    if (strstr(line, "GalaxySimulation")) ret++;
    if (strstr(line, "CoolingTest")) ret++;
    if (strstr(line, "ShearingBox")) ret++;
    if (strstr(line, "PoissonSolverTest")) ret++;
#ifdef TRANSFER
    if (strstr(line, "Radiative")           ) ret++;
    if (strstr(line, "PhotonTest")          ) ret++;
#endif

    if (strstr(line, "\"\"\"")              ) comment_count++;

    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#')
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s", line);
 
  }
 
  /* clean up */
 
  delete [] dummy;
  rewind(fptr);

  OutputTemperature = ((ProblemType == 7) || (ProblemType == 11));
 
  /* If we have turned on Comoving coordinates, read cosmology parameters. */
 
  if (ComovingCoordinates) {

    // Always output temperature in cosmology runs
    OutputTemperature = TRUE;

    if (CosmologyReadParameters(fptr, &MetaData.StopTime, &MetaData.Time)
	== FAIL) {
      fprintf(stderr, "Error in ReadCosmologyParameters.\n");;
      ENZO_FAIL("");
    }
    rewind(fptr);
  }
  else {
    if (ReadUnits(fptr) == FAIL){
      ENZO_FAIL("Error in ReadUnits. ");
    }
    rewind(fptr);
  }

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, 
    TimeUnits = 1.0, VelocityUnits = 1.0, PressureUnits = 1.0;
  double MassUnits = 1.0;
  if (UsePhysicalUnit) {
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, 
	     &MassUnits, MetaData.Time);
    PressureUnits = DensityUnits*pow(VelocityUnits,2);
  }

  /* Change input physical parameters into code units */

  double mh = 1.6726e-24;
  double uheat = pow(VelocityUnits,2)*2.0*mh/TimeUnits;
  PhotoelectricHeating /= uheat;
  StarMakerOverDensityThreshold /= DensityUnits;
  //  StarEnergyFeedbackRate = StarEnergyFeedbackRate/pow(LengthUnits,2)*pow(TimeUnits,3);

  SinkMergeDistance /= LengthUnits;
  SmallRho /= DensityUnits;
  SmallP /= PressureUnits;
  SmallT /= TemperatureUnits;
  MaximumAlvenSpeed /= VelocityUnits;
  float h, cs, dpdrho, dpde;
  EOS(SmallP, SmallRho, SmallEint, h, cs, dpdrho, dpde, EOSType, 1);
  if (debug && (HydroMethod == HD_RK || HydroMethod == MHD_RK))
    printf("smallrho=%g, smallp=%g, smalleint=%g, PressureUnits=%g, MaximumAlvenSpeed=%g\n",
	   SmallRho, SmallP, SmallEint, PressureUnits, MaximumAlvenSpeed);
  for (int i = 0; i < MAX_FLAGGING_METHODS; i++) {
    if (MinimumMassForRefinement[i] != FLOAT_UNDEFINED) {
      MinimumMassForRefinement[i] /= MassUnits;
    }
  }

  if (!ComovingCoordinates && UsePhysicalUnit) {
    for (int i = 0; i < MAX_FLAGGING_METHODS; i++) {
      if (MinimumOverDensityForRefinement[i] != FLOAT_UNDEFINED) {
	MinimumOverDensityForRefinement[i] /= DensityUnits;
      }
    }
  }

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

  if (MultiSpecies > 0)
    if (InitializeRateData(MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in InitializeRateData.");
    }
 
  if (MultiSpecies             == 0 && 
      MetalCooling             == 0 &&
      GadgetEquilibriumCooling == 0 &&
      RadiativeCooling          > 0) {
    if (InitializeEquilibriumCoolData(MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in InitializeEquilibriumCoolData.");
    }
  }

  /* If set, initialize CloudyCooling. */

  if (MetalCooling == CLOUDY_METAL_COOLING) {
    if (InitializeCloudyCooling(MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in InitializeCloudyCooling.");
    }
  }

  /* If using the internal radiation field, initialize it. */
 
  if (RadiationFieldType >= 10 && RadiationFieldType <= 11)
    if (InitializeRadiationFieldData(MetaData.Time) == FAIL) {
	ENZO_FAIL("Error in InitializeRadiationFieldData.");
      }
 
  /* Turn off DualEnergyFormalism for zeus hydro (and a few other things). */
 
  if (HydroMethod == Zeus_Hydro) {
    ConservativeInterpolation = FALSE;
    DualEnergyFormalism       = FALSE;
    //    FluxCorrection            = FALSE;
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
  case 0:  NSpecies = 0; break;
  case 1:  NSpecies = 5; break;
  case 2:  NSpecies = 8; break;
  case 3:  NSpecies = 11; break;
  default: NSpecies = 0; break;
  }

  // Determine color fields (NColor) later inside a grid object.
  // ...

  /* Set the number of particle attributes, if left unset. */
 
  if (NumberOfParticleAttributes == INT_UNDEFINED)
    if (StarParticleCreation || StarParticleFeedback)
      NumberOfParticleAttributes = 3;
    else
      NumberOfParticleAttributes = 0;
 
#ifdef UNUSED
  if (MaximumGravityRefinementLevel == INT_UNDEFINED)
    MaximumGravityRefinementLevel = (RadiativeCooling && SelfGravity
				     && HydroMethod == Zeus_Hydro) ?
       max(MaximumRefinementLevel-2, 5) : MaximumRefinementLevel;
#else
  if (MaximumGravityRefinementLevel == INT_UNDEFINED)
    MaximumGravityRefinementLevel = MaximumRefinementLevel;
#endif
 
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

  if ((MetaData.GravityBoundary != TopGridPeriodic) &&
      (UnigridTranspose)) {
    /* it turns out that Robert Harkness' unigrid transpose stuff is incompatible with the top
       grid isolated gravity boundary conditions.  I'm not 100 percent sure why this is - in the 
       meantime, just double-check to make sure that if one tries to use the isolated boundary
       conditions when the unigrid transpose stuff is on, the code crashes loudly.
       -- BWO, 26 June 2008 */
      if (MyProcessorNumber == ROOT_PROCESSOR){
	fprintf(stderr, "\n\n");
	fprintf(stderr, "  ************************************************************************\n");
	fprintf(stderr, "  ****  D'oh!  At present, you cannot use isolated top grid boundary  ****\n");
	fprintf(stderr, "  ****  conditions with the top grid unigrid bookkeeping scheme.      ****\n");
	fprintf(stderr, "  ****  Consult Brian O'Shea for the details of this wackiness,       ****\n");
	fprintf(stderr, "  ****  and in the meantime enzo DISABLED unigrid tranposition!       ****\n");
	fprintf(stderr, "  ************************************************************************\n");      
	fprintf(stderr, "\n\n");
      }
      UnigridTranspose = FALSE;
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
  bool MustRefineParticles = false;
  for (method = 0; method < MAX_FLAGGING_METHODS; method++)
    if (CellFlaggingMethod[method] == 8) {
      MustRefineParticles = true;
      break;
    }
  if (MustRefineParticles) {
    method = 0;
    while (CellFlaggingMethod[method] != INT_UNDEFINED)
      method++;
    CellFlaggingMethod[method] = 4;
  }

  if (TracerParticleOn) {
    ParticleTypeInFile = TRUE;
  }
 
  if (WritePotential && ComovingCoordinates && SelfGravity) {
    CopyGravPotential = TRUE;
  }

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
    if(getcwd(cwd_buffer, cwd_buffer_len) == NULL) {
      fprintf(stderr, "GETCWD call FAILED\n");
    }
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,"CWD %s\n", cwd_buffer);
    MetaData.GlobalDir = cwd_buffer;
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,"Global Dir set to %s\n", cwd_buffer);
  }
 
   for (int i=0; i<MetaData.TopGridRank;i++)
    TopGridDx[i]=(DomainRightEdge[i]-DomainLeftEdge[i])/MetaData.TopGridDims[i];

 //  for (int i=0; i<MetaData.TopGridRank; i++)
//      fprintf (stderr, "read  %"ISYM"  %"ISYM" \n", 
// 	      MetaData.LeftFaceBoundaryCondition[i], 
// 	      MetaData.RightFaceBoundaryCondition[i]);


   CheckShearingBoundaryConsistency(MetaData);
  return SUCCESS;
}
