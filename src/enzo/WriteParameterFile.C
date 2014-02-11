/***********************************************************************
/
/  WRITES A PARAMETER FILE (i.e. TopGrid data)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       February 29th, 2008
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine writes the parameter file in the argument and sets parameters
//   based on it.
 
#include <stdio.h>
#include <string.h>
#include <time.h>
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

 
/* function prototypes */
 
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int  CosmologyWriteParameters(FILE *fptr, FLOAT StopTime, FLOAT CurrentTime);
int  WriteUnits(FILE *fptr);
int  GetUnits(float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MAssUnits, FLOAT Time);
void get_uuid(char *buffer);
#ifdef TRANSFER
int RadiativeTransferWriteParameters(FILE *fptr);
int WritePhotonSources(FILE *fptr, FLOAT CurrentTime);
#endif /* TRANSFER */
int UpdateLocalDatabase(TopGridData &MetaData, int CurrentTimeID,
                        char *dset_uuid, char *Filename);
 
int WriteParameterFile(FILE *fptr, TopGridData &MetaData, char *name = NULL)
{
 
  MustRefineParticlesMinimumMass *= POW(1/(float(MetaData.TopGridDims[0])
				       *POW(float(RefineBy), float(MustRefineParticlesRefineToLevel))),3);

  int dim;
 
  /* Compute Units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
    VelocityUnits = 1;
  double MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits,  MetaData.Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
 
  float rhou = 1.0, lenu = 1.0, tempu = 1.0, tu = 1.0, velu = 1.0, presu = 1.0;
  double massu = 1.0;
  if (UsePhysicalUnit) {
    GetUnits(&rhou, &lenu, &tempu, &tu, &velu, &massu, MetaData.Time);
    presu = rhou*lenu*lenu/tu/tu;

    /* Change input physical parameters into real units */
    MustRefineParticlesMinimumMass *= massu;
    StarMakerOverDensityThreshold *= rhou;
    //  StarEnergyFeedbackRate = StarEnergyFeedbackRate/pow(LengthUnits,2)*pow(TimeUnits,3);
    
    if (SinkMergeDistance > 1.0)
      SinkMergeDistance *= lenu;
    SmallRho *= rhou;
    SmallP *= presu;
    SmallT *= tempu;
    MaximumAlvenSpeed *= velu;
    EOSSoundSpeed *=  velu;

    /*
    for (int i = 0; i < MAX_FLAGGING_METHODS; i++) {
      if (MinimumMassForRefinement[i] != FLOAT_UNDEFINED) {
	printf("i = %i, MinMass = %g, massu = %g\n",i,MinimumMassForRefinement[i],massu);
	MinimumMassForRefinement[i] *= massu;
	printf("i = %i, MinMass = %g\n",i,MinimumMassForRefinement[i]);
      }
    }
    */

    /* Check ReadParameterFile for the reason why this is commented out. 
       - Ji-hoon Kim in Apr.2010 */
    /*
    if (!ComovingCoordinates && UsePhysicalUnit) {
      for (int i = 0; i < MAX_FLAGGING_METHODS; i++) {
	if (MinimumOverDensityForRefinement[i] != FLOAT_UNDEFINED) 
	  MinimumOverDensityForRefinement[i] *= rhou;
      }
    }
    */

  }


  /* write data to Parameter output file */
 
  /* write MetaData parameters */
 
  fprintf(fptr, "InitialCycleNumber  = %"ISYM"\n", MetaData.CycleNumber);
  fprintf(fptr, "InitialTime         = %"GOUTSYM"\n", MetaData.Time);
  fprintf(fptr, "InitialCPUTime      = %"GSYM"\n\n", MetaData.CPUTime);
 
  fprintf(fptr, "CheckpointRestart   = %"ISYM"\n", CheckpointRestart);
  fprintf(fptr, "StopTime            = %"GOUTSYM"\n", MetaData.StopTime);
  fprintf(fptr, "StopCycle           = %"ISYM"\n", MetaData.StopCycle);
  fprintf(fptr, "StopSteps           = %"ISYM"\n", MetaData.StopSteps);
  fprintf(fptr, "StopCPUTime         = %lg\n", MetaData.StopCPUTime);
  fprintf(fptr, "ResubmitOn          = %"ISYM"\n", MetaData.ResubmitOn);
  fprintf(fptr, "ResubmitCommand     = %s\n\n", MetaData.ResubmitCommand);
 
  fprintf(fptr, "MaximumTopGridTimeStep = %"GSYM"\n", MetaData.MaximumTopGridTimeStep);

  fprintf(fptr, "TimeLastRestartDump = %"GSYM"\n", MetaData.TimeLastRestartDump);
  fprintf(fptr, "dtRestartDump       = %"GSYM"\n", MetaData.dtRestartDump);
  fprintf(fptr, "TimeLastDataDump    = %"GOUTSYM"\n", MetaData.TimeLastDataDump);
  fprintf(fptr, "dtDataDump          = %"GOUTSYM"\n", MetaData.dtDataDump);
  fprintf(fptr, "TimeLastHistoryDump = %"GOUTSYM"\n", MetaData.TimeLastHistoryDump);
  fprintf(fptr, "dtHistoryDump       = %"GOUTSYM"\n\n", MetaData.dtHistoryDump);
 
  fprintf(fptr, "TracerParticleOn           = %"ISYM"\n", TracerParticleOn);
  fprintf(fptr, "TracerParticleOutputVelocity           = %"ISYM"\n", TracerParticleOutputVelocity);

  fprintf(fptr, "TimeLastTracerParticleDump = %"GOUTSYM"\n",
          MetaData.TimeLastTracerParticleDump);
  fprintf(fptr, "dtTracerParticleDump       = %"GOUTSYM"\n",
          MetaData.dtTracerParticleDump);
  fprintf(fptr, "TimeLastInterpolatedDataDump    = %"GOUTSYM"\n", 
	  MetaData.TimeLastInterpolatedDataDump);
  fprintf(fptr, "dtInterpolatedDataDump          = %"GOUTSYM"\n", 
	  MetaData.dtInterpolatedDataDump);
 
  fprintf(fptr, "NewMovieLeftEdge     = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.NewMovieLeftEdge);
  fprintf(fptr, "NewMovieRightEdge    = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MetaData.NewMovieRightEdge);
  fprintf(fptr, "MovieSkipTimestep    = %"ISYM"\n", MovieSkipTimestep);
  fprintf(fptr, "Movie3DVolumes       = %"ISYM"\n", Movie3DVolumes);
  fprintf(fptr, "MovieVertexCentered  = %"ISYM"\n", MovieVertexCentered);
  fprintf(fptr, "NewMovieParticleOn   = %"ISYM"\n", NewMovieParticleOn);
  fprintf(fptr, "MovieDataField       = ");
  WriteListOfInts(fptr, MAX_MOVIE_FIELDS, MovieDataField);
  fprintf(fptr, "NewMovieDumpNumber   = %"ISYM"\n", NewMovieDumpNumber);
  fprintf(fptr, "NewMovieName         = %s\n", NewMovieName);
  fprintf(fptr, "MovieTimestepCounter = %"ISYM"\n", MetaData.MovieTimestepCounter);
  fprintf(fptr, "\n");

  fprintf(fptr, "CycleLastRestartDump = %"ISYM"\n", MetaData.CycleLastRestartDump);
  fprintf(fptr, "CycleSkipRestartDump = %"ISYM"\n", MetaData.CycleSkipRestartDump);
  fprintf(fptr, "CycleLastDataDump    = %"ISYM"\n", MetaData.CycleLastDataDump);
  fprintf(fptr, "CycleSkipDataDump    = %"ISYM"\n", MetaData.CycleSkipDataDump);
  fprintf(fptr, "CycleLastHistoryDump = %"ISYM"\n", MetaData.CycleLastHistoryDump);
  fprintf(fptr, "CycleSkipHistoryDump = %"ISYM"\n\n",
	  MetaData.CycleSkipHistoryDump);


  fprintf(fptr, "PythonTopGridSkip       = %"ISYM"\n", PythonTopGridSkip);
  fprintf(fptr, "PythonSubcycleSkip      = %"ISYM"\n", PythonSubcycleSkip);
  fprintf(fptr, "PythonReloadScript      = %"ISYM"\n", PythonReloadScript);
#ifdef USE_PYTHON
  fprintf(fptr, "NumberOfPythonCalls         = %"ISYM"\n", NumberOfPythonCalls);
  fprintf(fptr, "NumberOfPythonTopGridCalls  = %"ISYM"\n", NumberOfPythonTopGridCalls);
  fprintf(fptr, "NumberOfPythonSubcycleCalls = %"ISYM"\n", NumberOfPythonSubcycleCalls);
#endif

  fprintf(fptr, "TimingCycleSkip             = %"ISYM"\n", TimingCycleSkip);

  fprintf(fptr, "CycleSkipGlobalDataDump = %"ISYM"\n\n", //AK
          MetaData.CycleSkipGlobalDataDump);

  fprintf(fptr, "SubcycleNumber          = %"ISYM"\n", MetaData.SubcycleNumber);
  fprintf(fptr, "SubcycleSkipDataDump    = %"ISYM"\n", MetaData.SubcycleSkipDataDump);
  fprintf(fptr, "SubcycleLastDataDump    = %"ISYM"\n", MetaData.SubcycleLastDataDump);
  fprintf(fptr, "OutputFirstTimeAtLevel = %"ISYM"\n",
	  MetaData.OutputFirstTimeAtLevel);
  fprintf(fptr, "StopFirstTimeAtLevel    = %"ISYM"\n",
	  MetaData.StopFirstTimeAtLevel);
  fprintf(fptr, "NumberOfOutputsBeforeExit = %"ISYM"\n\n",
	  MetaData.NumberOfOutputsBeforeExit);

  fprintf(fptr, "OutputOnDensity = %"ISYM"\n", OutputOnDensity);
  fprintf(fptr, "StartDensityOutputs = %"GSYM"\n", StartDensityOutputs);
  fprintf(fptr, "CurrentDensityOutput = %"GSYM"\n", CurrentDensityOutput);
  fprintf(fptr, "IncrementDensityOutput = %"GSYM"\n\n", IncrementDensityOutput);

  fprintf(fptr, "FileDirectedOutput = %"ISYM"\n", FileDirectedOutput);

  fprintf(fptr, "HierarchyFileInputFormat = %"ISYM"\n", HierarchyFileInputFormat);
  fprintf(fptr, "HierarchyFileOutputFormat = %"ISYM"\n", HierarchyFileOutputFormat);
 
  fprintf(fptr, "RestartDumpNumber   = %"ISYM"\n", MetaData.RestartDumpNumber);
  fprintf(fptr, "DataDumpNumber      = %"ISYM"\n", MetaData.DataDumpNumber);
  fprintf(fptr, "HistoryDumpNumber   = %"ISYM"\n", MetaData.HistoryDumpNumber);
  fprintf(fptr, "TracerParticleDumpNumber = %"ISYM"\n",
          MetaData.TracerParticleDumpNumber);
 
  fprintf(fptr, "RestartDumpName     = %s\n", MetaData.RestartDumpName);
  fprintf(fptr, "DataDumpName        = %s\n", MetaData.DataDumpName);
  fprintf(fptr, "HistoryDumpName     = %s\n", MetaData.HistoryDumpName);
  fprintf(fptr, "TracerParticleDumpName = %s\n",
          MetaData.TracerParticleDumpName);
  fprintf(fptr, "RedshiftDumpName    = %s\n\n", MetaData.RedshiftDumpName);
 
  if (MetaData.RestartDumpDir != NULL)
    fprintf(fptr, "RestartDumpDir        = %s\n", MetaData.RestartDumpDir);
  if (MetaData.DataDumpDir != NULL)
    fprintf(fptr, "DataDumpDir           = %s\n", MetaData.DataDumpDir);
  if (MetaData.HistoryDumpDir != NULL)
    fprintf(fptr, "HistoryDumpDir        = %s\n", MetaData.HistoryDumpDir);
  if (MetaData.TracerParticleDumpDir != NULL)
    fprintf(fptr, "TracerParticleDumpDir = %s\n", MetaData.TracerParticleDumpDir);
  if (MetaData.RedshiftDumpDir != NULL)
    fprintf(fptr, "RedshiftDumpDir       = %s\n\n", MetaData.RedshiftDumpDir);
 
  if (MetaData.LocalDir != NULL)
    fprintf(fptr, "LocalDir            = %s\n", MetaData.LocalDir);
  if (MetaData.GlobalDir != NULL)
    fprintf(fptr, "GlobalDir           = %s\n", MetaData.GlobalDir);
 
  for (dim = 0; dim < MAX_CUBE_DUMPS; dim++)
    if (CubeDumps[dim] != NULL)
      fprintf(fptr, "CubeDump[%"ISYM"]            = %s\n", dim, CubeDumps[dim]);

  fprintf(fptr, "NumberOfGhostZones    = %"ISYM"\n", NumberOfGhostZones);
  fprintf(fptr, "LoadBalancing          = %"ISYM"\n", LoadBalancing);
  fprintf(fptr, "ResetLoadBalancing     = %"ISYM"\n", ResetLoadBalancing);
  fprintf(fptr, "LoadBalancingCycleSkip = %"ISYM"\n", LoadBalancingCycleSkip);
  fprintf(fptr, "LoadBalancingMinLevel  = %"ISYM"\n", LoadBalancingMinLevel);
  fprintf(fptr, "LoadBalancingMaxLevel  = %"ISYM"\n", LoadBalancingMaxLevel);
 
  fprintf(fptr, "ConductionDynamicRebuildHierarchy = %"ISYM"\n", ConductionDynamicRebuildHierarchy);
  fprintf(fptr, "ConductionDynamicRebuildMinLevel  = %"ISYM"\n", ConductionDynamicRebuildMinLevel);
  for (dim = 0;dim < MAX_DEPTH_OF_HIERARCHY;dim++) {
    if (RebuildHierarchyCycleSkip[dim] != 1) {
      fprintf(fptr, "RebuildHierarchyCycleSkip[%"ISYM"] = %"ISYM"\n",
	      dim, RebuildHierarchyCycleSkip[dim]);
    }
  }

  for (dim = 0; dim < MAX_TIME_ACTIONS; dim++)
    if (TimeActionType[dim] > 0) {
      fprintf(fptr, "TimeActionType[%"ISYM"]      = %"ISYM"\n", dim,TimeActionType[dim]);
      if (ComovingCoordinates)
	fprintf(fptr, "TimeActionRedshift[%"ISYM"]  = %"GOUTSYM"\n", dim,
		TimeActionRedshift[dim]);
      else
	fprintf(fptr, "TimeActionTime[%"ISYM"]      = %"GOUTSYM"\n", dim,
		TimeActionRedshift[dim]);
      fprintf(fptr, "TimeActionParameter[%"ISYM"] = %"GSYM"\n", dim,
	      TimeActionParameter[dim]);
    }
 
  fprintf(fptr, "StaticHierarchy     = %"ISYM"\n", MetaData.StaticHierarchy);
 
  fprintf(fptr, "TopGridRank         = %"ISYM"\n", MetaData.TopGridRank);
  fprintf(fptr, "TopGridDimensions   = ");
  WriteListOfInts(fptr, MetaData.TopGridRank, MetaData.TopGridDims);
  fprintf(fptr, "\n");
 
  fprintf(fptr, "TopGridGravityBoundary = %"ISYM"\n", MetaData.GravityBoundary);

#ifdef TRANSFER
  if (MetaData.RadHydroParameterFname != NULL) 
    fprintf(fptr, "RadHydroParamfile = %s\n", MetaData.RadHydroParameterFname);
#endif
  fprintf(fptr, "ImplicitProblem = %"ISYM"\n", ImplicitProblem);
#ifdef EMISSIVITY
  fprintf(fptr, "StarMakerEmissivityField = %"ISYM"\n", StarMakerEmissivityField);
  fprintf(fptr, "uv_param = %"GSYM"\n", uv_param);
#endif
  fprintf(fptr, "RadiativeTransferFLD   = %"ISYM"\n", RadiativeTransferFLD);

  fprintf(fptr, "ParticleBoundaryType   = %"ISYM"\n",MetaData.ParticleBoundaryType);
  fprintf(fptr, "NumberOfParticles      = %"PISYM" (do not modify)\n",
	  MetaData.NumberOfParticles);
 
  fprintf(fptr, "CourantSafetyNumber    = %"FSYM"\n",
	  MetaData.CourantSafetyNumber);
  fprintf(fptr, "PPMFlatteningParameter = %"ISYM"\n",
	  MetaData.PPMFlatteningParameter);
  fprintf(fptr, "PPMDiffusionParameter  = %"ISYM"\n",
	  MetaData.PPMDiffusionParameter);
  fprintf(fptr, "PPMSteepeningParameter = %"ISYM"\n\n",
	  MetaData.PPMSteepeningParameter);
 
  /* write global Parameters */
 
  fprintf(fptr, "ProblemType                    = %"ISYM"\n", ProblemType);
#ifdef NEW_PROBLEM_TYPES
  fprintf(fptr, "ProblemTypeName                = %s\n", ProblemTypeName);
#endif
  fprintf(fptr, "HydroMethod                    = %"ISYM"\n", HydroMethod);
  fprintf(fptr, "huge_number                    = %e\n", huge_number);
  fprintf(fptr, "tiny_number                    = %e\n", tiny_number);
  fprintf(fptr, "Gamma                          = %"GOUTSYM"\n", Gamma);
  fprintf(fptr, "PressureFree                   = %"ISYM"\n", PressureFree);
  fprintf(fptr, "RefineBy                       = %"ISYM"\n", RefineBy);
  fprintf(fptr, "MaximumRefinementLevel         = %"ISYM"\n", MaximumRefinementLevel);
  fprintf(fptr, "MaximumGravityRefinementLevel  = %"ISYM"\n",
	  MaximumGravityRefinementLevel);
  fprintf(fptr, "MaximumParticleRefinementLevel = %"ISYM"\n",
	  MaximumParticleRefinementLevel);
  fprintf(fptr, "CellFlaggingMethod             = ");
  WriteListOfInts(fptr, MAX_FLAGGING_METHODS, CellFlaggingMethod);
  fprintf(fptr, "FluxCorrection                 = %"ISYM"\n", FluxCorrection);
  fprintf(fptr, "InterpolationMethod            = %"ISYM"\n", InterpolationMethod);
  fprintf(fptr, "ConservativeInterpolation      = %"ISYM"\n", ConservativeInterpolation);
  fprintf(fptr, "MinimumEfficiency              = %"GSYM"\n", MinimumEfficiency);
  fprintf(fptr, "SubgridSizeAutoAdjust          = %"ISYM"\n", SubgridSizeAutoAdjust);
  fprintf(fptr, "OptimalSubgridsPerProcessor    = %"ISYM"\n", 
	  OptimalSubgridsPerProcessor);
  fprintf(fptr, "MinimumSubgridEdge             = %"ISYM"\n", MinimumSubgridEdge);
  fprintf(fptr, "MaximumSubgridSize             = %"ISYM"\n", MaximumSubgridSize);
  fprintf(fptr, "NumberOfBufferZones            = %"ISYM"\n\n", NumberOfBufferZones);

  fprintf(fptr, "FastSiblingLocatorEntireDomain      = %"ISYM"\n", 
	  FastSiblingLocatorEntireDomain);
  fprintf(fptr, "MustRefineRegionMinRefinementLevel  = %"ISYM"\n", 
	  MustRefineRegionMinRefinementLevel);
  fprintf(fptr, "MetallicityRefinementMinLevel       = %"ISYM"\n", 
	  MetallicityRefinementMinLevel);
  fprintf(fptr, "MetallicityRefinementMinMetallicity = %"GSYM"\n", 
	  MetallicityRefinementMinMetallicity);
  fprintf(fptr, "MetallicityRefinementMinDensity     = %"GSYM"\n", 
	  MetallicityRefinementMinDensity);
 
  fprintf(fptr, "DomainLeftEdge         = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, DomainLeftEdge);
  fprintf(fptr, "DomainRightEdge        = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, DomainRightEdge);
  fprintf(fptr, "GridVelocity           = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, GridVelocity);
  fprintf(fptr, "RefineRegionAutoAdjust = %"ISYM"\n", RefineRegionAutoAdjust);
  fprintf(fptr, "RefineRegionLeftEdge   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, RefineRegionLeftEdge);
  fprintf(fptr, "RefineRegionRightEdge  = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, RefineRegionRightEdge);
  fprintf(fptr, "MustRefineRegionLeftEdge   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MustRefineRegionLeftEdge);
  fprintf(fptr, "MustRefineRegionRightEdge  = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, MustRefineRegionRightEdge);
  fprintf(fptr, "RefineRegionTimeType   = %d\n", RefineRegionTimeType);
  if (RefineRegionFile != NULL)
    fprintf(fptr, "RefineRegionFile       = %s\n", RefineRegionFile);
  fprintf(fptr, "\n");
  fprintf(fptr, "\n");

  if (DatabaseLocation != NULL)
    fprintf(fptr, "DatabaseLocation       = %s\n", DatabaseLocation);
  fprintf(fptr, "\n");
  fprintf(fptr, "\n");
 
  for (dim = 0; dim < MAX_NUMBER_OF_BARYON_FIELDS; dim++) {
    if (DataLabel[dim])
      fprintf(fptr, "DataLabel[%"ISYM"]              = %s\n", dim, DataLabel[dim]);
    if (DataUnits[dim])
      fprintf(fptr, "DataUnits[%"ISYM"]              = %s\n", dim, DataUnits[dim]);
    if (DataLabel[dim]) {
      if ((strstr(DataLabel[dim], "Density") != NULL) ||
	  (strstr(DataLabel[dim], "Colour") != NULL))
	fprintf(fptr, "#DataCGSConversionFactor[%"ISYM"] = %"GSYM"\n", dim, DensityUnits);
      if (strstr(DataLabel[dim], "velocity") != NULL)
	fprintf(fptr, "#DataCGSConversionFactor[%"ISYM"] = %"GSYM"\n", dim, VelocityUnits);
    }
  }
  fprintf(fptr, "#TimeUnits                 = %"GSYM"\n", TimeUnits);
  fprintf(fptr, "#TemperatureUnits          = %"GSYM"\n", TemperatureUnits);
  fprintf(fptr, "\n");
 
  fprintf(fptr, "UniformGravity             = %"ISYM"\n", UniformGravity);
  fprintf(fptr, "UniformGravityDirection    = %"ISYM"\n", UniformGravityDirection);
  fprintf(fptr, "UniformGravityConstant     = %"GSYM"\n", UniformGravityConstant);
 
  fprintf(fptr, "PointSourceGravity           = %"ISYM"\n",PointSourceGravity);
  fprintf(fptr, "PointSourceGravityPosition   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, PointSourceGravityPosition);
  fprintf(fptr, "PointSourceGravityConstant   = %"GSYM"\n",
	  PointSourceGravityConstant);
  fprintf(fptr, "PointSourceGravityCoreRadius = %"GSYM"\n\n",
	  PointSourceGravityCoreRadius);

  fprintf(fptr, "ExternalGravity           = %"ISYM"\n",ExternalGravity); 
  fprintf(fptr, "ExternalGravityConstant     = %"FSYM"\n",ExternalGravityConstant);
  fprintf(fptr, "ExternalGravityRadius     = %"FSYM"\n",ExternalGravityRadius); 
  fprintf(fptr, "ExternalGravityDensity     = %"FSYM"\n",ExternalGravityDensity);
  fprintf(fptr, "ExternalGravityPosition   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, ExternalGravityPosition);
  fprintf(fptr, "ExternalGravityOrientation   = ");
  WriteListOfFloats(fptr, MetaData.TopGridRank, ExternalGravityOrientation);

  fprintf(fptr, "SelfGravity                    = %"ISYM"\n", SelfGravity);
  fprintf(fptr, "SelfGravityGasOff              = %"ISYM"\n", SelfGravityGasOff);
  fprintf(fptr, "AccretionKernal                = %"ISYM"\n", AccretionKernal);
  fprintf(fptr, "GravitationalConstant          = %e\n",
	  GravitationalConstant);
  fprintf(fptr, "S2ParticleSize                 = %"GSYM"\n", S2ParticleSize);
  fprintf(fptr, "GravityResolution              = %"GSYM"\n", GravityResolution);
  fprintf(fptr, "ComputePotential               = %"ISYM"\n", ComputePotential);
  fprintf(fptr, "PotentialIterations            = %"ISYM"\n", PotentialIterations);
  fprintf(fptr, "WritePotential                 = %"ISYM"\n", WritePotential);
  fprintf(fptr, "ParticleSubgridDepositMode     = %"ISYM"\n", ParticleSubgridDepositMode);
  fprintf(fptr, "BaryonSelfGravityApproximation = %"ISYM"\n",
	  BaryonSelfGravityApproximation);

  fprintf(fptr, "InlineHaloFinder               = %"ISYM"\n", InlineHaloFinder);
  fprintf(fptr, "HaloFinderSubfind              = %"ISYM"\n", HaloFinderSubfind);
  fprintf(fptr, "HaloFinderCycleSkip            = %"ISYM"\n", 
	  HaloFinderCycleSkip);
  fprintf(fptr, "HaloFinderRunAfterOutput       = %"ISYM"\n", 
	  HaloFinderRunAfterOutput);
  fprintf(fptr, "HaloFinderOutputParticleList   = %"ISYM"\n", 
	  HaloFinderOutputParticleList);
  fprintf(fptr, "HaloFinderMinimumSize          = %"ISYM"\n", 
	  HaloFinderMinimumSize);
  fprintf(fptr, "HaloFinderLinkingLength        = %"FSYM"\n",
	  HaloFinderLinkingLength);
  fprintf(fptr, "HaloFinderTimestep             = %"FSYM"\n",
	  HaloFinderTimestep);
  fprintf(fptr, "HaloFinderLastTime             = %"PSYM"\n\n", 
	  HaloFinderLastTime);
 
  fprintf(fptr, "GreensFunctionMaxNumber     = %"ISYM"\n", GreensFunctionMaxNumber);
  fprintf(fptr, "GreensFunctionMaxSize       = %"ISYM"\n", GreensFunctionMaxSize);
 
  fprintf(fptr, "DualEnergyFormalism         = %"ISYM"\n", DualEnergyFormalism);
  fprintf(fptr, "DualEnergyFormalismEta1     = %e\n", DualEnergyFormalismEta1);
  fprintf(fptr, "DualEnergyFormalismEta2     = %e\n", DualEnergyFormalismEta2);
  fprintf(fptr, "ParticleCourantSafetyNumber = %"FSYM"\n\n", ParticleCourantSafetyNumber);
  fprintf(fptr, "RootGridCourantSafetyNumber = %"FSYM"\n\n", RootGridCourantSafetyNumber);
  fprintf(fptr, "RandomForcing                  = %"ISYM"\n", RandomForcing);
  fprintf(fptr, "RandomForcingEdot              = %"GSYM"\n", RandomForcingEdot);
#ifdef USE_GRACKLE
  /* Grackle chemistry parameters */
  fprintf(fptr, "use_grackle                 = %"ISYM"\n", grackle_chemistry.use_grackle);
  fprintf(fptr, "with_radiative_cooling      = %"ISYM"\n", grackle_chemistry.with_radiative_cooling);
  fprintf(fptr, "grackle_data_file           = %s\n", grackle_chemistry.grackle_data_file);
  fprintf(fptr, "UVbackground                = %"ISYM"\n", grackle_chemistry.UVbackground);
  fprintf(fptr, "Compton_xray_heating        = %"ISYM"\n", grackle_chemistry.Compton_xray_heating);
  fprintf(fptr, "LWbackground_intensity      = %"FSYM"\n", grackle_chemistry.LWbackground_intensity);
  fprintf(fptr, "LWbackground_sawtooth_suppression = %"ISYM"\n", grackle_chemistry.LWbackground_sawtooth_suppression);
  /********************************/
#endif
  fprintf(fptr, "RadiativeCooling               = %"ISYM"\n", RadiativeCooling);
  fprintf(fptr, "RadiativeCoolingModel          = %"ISYM"\n", RadiativeCoolingModel);
  fprintf(fptr, "GadgetEquilibriumCooling       = %"ISYM"\n", GadgetEquilibriumCooling);
  fprintf(fptr, "MultiSpecies                   = %"ISYM"\n", MultiSpecies);
  fprintf(fptr, "CIECooling                     = %"ISYM"\n", CIECooling);
  fprintf(fptr, "H2OpticalDepthApproximation    = %"ISYM"\n", H2OpticalDepthApproximation);
  fprintf(fptr, "ThreeBodyRate                  = %"ISYM"\n", ThreeBodyRate);
  fprintf(fptr, "H2FormationOnDust              = %"ISYM"\n", H2FormationOnDust);
  fprintf(fptr, "CloudyCoolingGridFile          = %s\n", CloudyCoolingData.CloudyCoolingGridFile);
  fprintf(fptr, "IncludeCloudyHeating           = %"ISYM"\n", CloudyCoolingData.IncludeCloudyHeating);
  fprintf(fptr, "CMBTemperatureFloor            = %"ISYM"\n", CloudyCoolingData.CMBTemperatureFloor);
  fprintf(fptr, "CloudyElectronFractionFactor   = %"FSYM"\n", CloudyCoolingData.CloudyElectronFractionFactor);
  fprintf(fptr, "MetalCooling                   = %"ISYM"\n", MetalCooling);
  fprintf(fptr, "MetalCoolingTable              = %s\n", MetalCoolingTable);
  fprintf(fptr, "RadiativeTransfer              = %"ISYM"\n", RadiativeTransfer);
  fprintf(fptr, "RadiationXRaySecondaryIon      = %"ISYM"\n", RadiationXRaySecondaryIon);
  fprintf(fptr, "RadiationXRayComptonHeating    = %"ISYM"\n", RadiationXRayComptonHeating);
  fprintf(fptr, "ShockMethod                    = %"ISYM"\n", ShockMethod);
  fprintf(fptr, "ShockTemperatureFloor          = %"FSYM"\n", ShockTemperatureFloor);
  fprintf(fptr, "StorePreShockFields            = %"ISYM"\n", StorePreShockFields);
  fprintf(fptr, "FindShocksOnlyOnOutput         = %"ISYM"\n", FindShocksOnlyOnOutput);
  fprintf(fptr, "RadiationFieldType             = %"ISYM"\n", RadiationFieldType);
  fprintf(fptr, "TabulatedLWBackground          = %"ISYM"\n", TabulatedLWBackground);
  fprintf(fptr, "AdjustUVBackground             = %"ISYM"\n", AdjustUVBackground);
  fprintf(fptr, "AdjustUVBackgroundHighRedshift = %"ISYM"\n", AdjustUVBackgroundHighRedshift);
  fprintf(fptr, "SetUVBAmplitude                = %"GSYM"\n", SetUVBAmplitude);
  fprintf(fptr, "SetHeIIHeatingScale            = %"GSYM"\n", SetHeIIHeatingScale);
  fprintf(fptr, "RadiationFieldLevelRecompute   = %"ISYM"\n", RadiationFieldLevelRecompute);
  fprintf(fptr, "RadiationFieldRedshift         = %"FSYM"\n", RadiationFieldRedshift);
  fprintf(fptr, "RadiationShield                = %"ISYM"\n", RadiationData.RadiationShield);
  fprintf(fptr, "RadiationSpectrumNormalization = %"GSYM"\n", CoolData.f3);
  fprintf(fptr, "RadiationSpectrumSlope         = %"GSYM"\n", CoolData.alpha0);
  fprintf(fptr, "CoolDataf0to3                  = %"FSYM"\n", CoolData.f0to3);
  fprintf(fptr, "RadiationRedshiftOn            = %"FSYM"\n", CoolData.RadiationRedshiftOn);
  fprintf(fptr, "RadiationRedshiftOff           = %"FSYM"\n", CoolData.RadiationRedshiftOff);
  fprintf(fptr, "RadiationRedshiftFullOn        = %"FSYM"\n", CoolData.RadiationRedshiftFullOn);
  fprintf(fptr, "RadiationRedshiftDropOff       = %"FSYM"\n", CoolData.RadiationRedshiftDropOff);
  fprintf(fptr, "HydrogenFractionByMass         = %"FSYM"\n", CoolData.HydrogenFractionByMass);
  fprintf(fptr, "DeuteriumToHydrogenRatio       = %"FSYM"\n", CoolData.DeuteriumToHydrogenRatio);
  fprintf(fptr, "SolarMetalFractionByMass       = %"FSYM"\n", CoolData.SolarMetalFractionByMass);
  fprintf(fptr, "NumberOfTemperatureBins        = %"ISYM"\n", CoolData.NumberOfTemperatureBins);
  fprintf(fptr, "CoolDataIh2co                  = %"ISYM"\n", CoolData.ih2co);
  fprintf(fptr, "CoolDataIpiht                  = %"ISYM"\n", CoolData.ipiht);
  fprintf(fptr, "TemperatureStart               = %"FSYM"\n", CoolData.TemperatureStart);
  fprintf(fptr, "TemperatureEnd                 = %"FSYM"\n", CoolData.TemperatureEnd);
  fprintf(fptr, "CoolDataCompXray               = %"FSYM"\n", CoolData.comp_xray);
  fprintf(fptr, "CoolDataTempXray               = %"FSYM"\n", CoolData.temp_xray);
  fprintf(fptr, "NumberOfDustTemperatureBins    = %"ISYM"\n", RateData.NumberOfDustTemperatureBins);
  fprintf(fptr, "DustTemperatureStart           = %"FSYM"\n", RateData.DustTemperatureStart);
  fprintf(fptr, "DustTemperatureEnd             = %"FSYM"\n", RateData.DustTemperatureEnd);
  fprintf(fptr, "PhotoelectricHeating           = %"ISYM"\n", PhotoelectricHeating);
  fprintf(fptr, "PhotoelectricHeatingRate       = %"GSYM"\n", PhotoelectricHeatingRate);
  
  fprintf(fptr, "VelAnyl                        = %"ISYM"\n", VelAnyl);
  fprintf(fptr, "BAnyl                          = %"ISYM"\n", BAnyl);

  // Negative number means that it was flagged from the command line.  Don't propagate.
  if (OutputCoolingTime < 0)
    fprintf(fptr, "OutputCoolingTime              = %"ISYM"\n", 0);
  else
    fprintf(fptr, "OutputCoolingTime              = %"ISYM"\n", OutputCoolingTime);
  fprintf(fptr, "OutputTemperature              = %"ISYM"\n", OutputTemperature);
  fprintf(fptr, "OutputDustTemperature          = %"ISYM"\n", OutputDustTemperature);

  // Negative number means that it was flagged from the command line.  Don't propagate.
  if (OutputSmoothedDarkMatter < 0)
    fprintf(fptr, "OutputSmoothedDarkMatter       = %"ISYM"\n", 0);
  else
    fprintf(fptr, "OutputSmoothedDarkMatter       = %"ISYM"\n", 
	    OutputSmoothedDarkMatter);
  fprintf(fptr, "SmoothedDarkMatterNeighbors    = %"ISYM"\n", 
	  SmoothedDarkMatterNeighbors);
  fprintf(fptr, "OutputGriddedStarParticle      = %"ISYM"\n", 
	  OutputGriddedStarParticle);
 
  fprintf(fptr, "ZEUSLinearArtificialViscosity    = %"GSYM"\n",
	  ZEUSLinearArtificialViscosity);
  fprintf(fptr, "ZEUSQuadraticArtificialViscosity = %"GSYM"\n",
	  ZEUSQuadraticArtificialViscosity);
  fprintf(fptr, "UseMinimumPressureSupport        = %"ISYM"\n",
	  UseMinimumPressureSupport);
  fprintf(fptr, "MinimumPressureSupportParameter  = %"FSYM"\n",
	  MinimumPressureSupportParameter);
  fprintf(fptr, "RefineByJeansLengthSafetyFactor  = %"FSYM"\n",
	  RefineByJeansLengthSafetyFactor);
  fprintf(fptr, "JeansRefinementColdTemperature  = %"FSYM"\n",
	  JeansRefinementColdTemperature);
  fprintf(fptr, "RefineByResistiveLengthSafetyFactor  = %"FSYM"\n", 
	  RefineByResistiveLengthSafetyFactor);
  fprintf(fptr, "MustRefineParticlesRefineToLevel = %"ISYM"\n",
          MustRefineParticlesRefineToLevel);
  fprintf(fptr, "MustRefineParticlesRefineToLevelAutoAdjust = %"ISYM"\n",
          MustRefineParticlesRefineToLevelAutoAdjust);
  fprintf(fptr, "MustRefineParticlesMinimumMass = %"FSYM"\n",
          MustRefineParticlesMinimumMass);
  fprintf(fptr, "ParticleTypeInFile               = %"ISYM"\n",
          ParticleTypeInFile);
  fprintf(fptr, "WriteGhostZones                  = %"ISYM"\n",
          WriteGhostZones);
  fprintf(fptr, "ReadGhostZones                   = %"ISYM"\n",
          ReadGhostZones);
  fprintf(fptr, "OutputParticleTypeGrouping       = %"ISYM"\n",
          OutputParticleTypeGrouping);
  fprintf(fptr, "MoveParticlesBetweenSiblings     = %"ISYM"\n",
	  MoveParticlesBetweenSiblings);
  fprintf(fptr, "ParticleSplitterIterations       = %"ISYM"\n",
	  ParticleSplitterIterations);
  fprintf(fptr, "ParticleSplitterChildrenParticleSeparation     = %"FSYM"\n",
	  ParticleSplitterChildrenParticleSeparation);
  fprintf(fptr, "ResetMagneticField               = %"ISYM"\n",
	  ResetMagneticField);
  fprintf(fptr, "ResetMagneticFieldAmplitude      = %"GSYM" %"GSYM" %"GSYM"\n", 
	  ResetMagneticFieldAmplitude[0],
	  ResetMagneticFieldAmplitude[1],
	  ResetMagneticFieldAmplitude[2]);

  for (int ireg = 0; ireg < MAX_STATIC_REGIONS; ireg++){
    if (AvoidRefineRegionLevel[ireg] != INT_UNDEFINED) {
      fprintf(fptr, "AvoidRefineRegionLevel[%"ISYM"] = %"ISYM"\n", ireg,
	      AvoidRefineRegionLevel[ireg]);
      fprintf(fptr, "AvoidRefineRegionLeftEdge[%"ISYM"] = ", ireg);
      WriteListOfFloats(fptr, MAX_DIMENSION, AvoidRefineRegionLeftEdge[ireg]);
      fprintf(fptr, "AvoidRefineRegionRightEdge[%"ISYM"] = ", ireg);
      WriteListOfFloats(fptr, MAX_DIMENSION, AvoidRefineRegionRightEdge[ireg]);
    }
  }


  fprintf(fptr, "MultiRefineRegionMaximumOuterLevel  = %"ISYM"\n",
          MultiRefineRegionMaximumOuterLevel);
  fprintf(fptr, "MultiRefineRegionMinimumOuterLevel  = %"ISYM"\n",
          MultiRefineRegionMinimumOuterLevel);
  
  for (int ireg = 0; ireg < MAX_STATIC_REGIONS; ireg++){
    if (MultiRefineRegionGeometry[ireg] >= 0) {
      fprintf(fptr, "MultiRefineRegionMaximumLevel[%"ISYM"] = %"ISYM"\n", ireg,
              MultiRefineRegionMaximumLevel[ireg]);

      fprintf(fptr, "MultiRefineRegionMinimumLevel[%"ISYM"] = %"ISYM"\n", ireg,
              MultiRefineRegionMinimumLevel[ireg]);

      fprintf(fptr, "MultiRefineRegionGeometry[%"ISYM"] = %"ISYM"\n", ireg,
              MultiRefineRegionGeometry[ireg]);

      fprintf(fptr, "MultiRefineRegionRadius[%"ISYM"] = %"GSYM"\n", ireg,
              MultiRefineRegionRadius[ireg]);

      fprintf(fptr, "MultiRefineRegionWidth[%"ISYM"] = %"GSYM"\n",
              MultiRefineRegionWidth[ireg]);

      fprintf(fptr, "MultiRefineRegionStaggeredRefinement[%"ISYM"] =%"GSYM"\n",
              MultiRefineRegionStaggeredRefinement[ireg]);

      fprintf(fptr, "MultiRefineRegionLeftEdge[%"ISYM"] = ", ireg);
      WriteListOfFloats(fptr, MAX_DIMENSION, MultiRefineRegionLeftEdge[ireg]);

      fprintf(fptr, "MultiRefineRegionRightEdge[%"ISYM"] = ", ireg);
      WriteListOfFloats(fptr, MAX_DIMENSION, MultiRefineRegionRightEdge[ireg]);

      fprintf(fptr, "MultiRefineRegionCenter[%"ISYM"] = ", ireg);
      WriteListOfFloats(fptr, MAX_DIMENSION, MultiRefineRegionCenter[ireg]);

      fprintf(fptr, "MultiRefineRegionOrientation[%"ISYM"] = ", ireg);
      WriteListOfFloats(fptr, MAX_DIMENSION, MultiRefineRegionOrientation[ireg]);

      fprintf(fptr, "\n");
    }
  }

  for (int ireg = 0; ireg < MAX_STATIC_REGIONS; ireg++){
    if (StaticRefineRegionLevel[ireg] != INT_UNDEFINED) {
      fprintf(fptr, "StaticRefineRegionLevel[%"ISYM"] = %"ISYM"\n", ireg,
	      StaticRefineRegionLevel[ireg]);
      fprintf(fptr, "StaticRefineRegionLeftEdge[%"ISYM"] = ", ireg);
      WriteListOfFloats(fptr, MAX_DIMENSION, StaticRefineRegionLeftEdge[ireg]);
      fprintf(fptr, "StaticRefineRegionRightEdge[%"ISYM"] = ", ireg);
      WriteListOfFloats(fptr, MAX_DIMENSION, StaticRefineRegionRightEdge[ireg]);
      fprintf(fptr, "\n");
    }
  }
 
  fprintf(fptr, "ParallelRootGridIO              = %"ISYM"\n", ParallelRootGridIO);
  fprintf(fptr, "ParallelParticleIO              = %"ISYM"\n", ParallelParticleIO);
  fprintf(fptr, "Unigrid                         = %"ISYM"\n", Unigrid);
  fprintf(fptr, "UnigridTranspose                = %"ISYM"\n", UnigridTranspose);
  fprintf(fptr, "NumberOfRootGridTilesPerDimensionPerProcessor = %"ISYM"\n", 
	  NumberOfRootGridTilesPerDimensionPerProcessor);
  fprintf(fptr, "PartitionNestedGrids            = %"ISYM"\n", PartitionNestedGrids);
  fprintf(fptr, "ExtractFieldsOnly               = %"ISYM"\n", ExtractFieldsOnly);
  fprintf(fptr, "CubeDumpEnabled                 = %"ISYM"\n", CubeDumpEnabled);
 
  fprintf(fptr, "Debug1                          = %"ISYM"\n", debug1);
  fprintf(fptr, "Debug2                          = %"ISYM"\n", debug2);

  fprintf(fptr, "MemoryLimit                     = %lld\n", MemoryLimit);

#ifdef STAGE_INPUT
  fprintf(fptr, "StageInput                      = %"ISYM"\n", StageInput);
  fprintf(fptr, "LocalPath                       = %s\n", LocalPath);
  fprintf(fptr, "GlobalPath                      = %s\n", GlobalPath);
#endif

#ifdef OOC_BOUNDARY

  fprintf(fptr, "ExternalBoundaryIO              = %"ISYM"\n", ExternalBoundaryIO);
  fprintf(fptr, "ExternalBoundaryTypeIO          = %"ISYM"\n", ExternalBoundaryTypeIO);
  fprintf(fptr, "ExternalBoundaryValueIO         = %"ISYM"\n", ExternalBoundaryValueIO);
  fprintf(fptr, "SimpleConstantBoundary          = %"ISYM"\n", SimpleConstantBoundary);

#endif
 

  fprintf(fptr, "SlopeFlaggingFields ="
	  " %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",
	  SlopeFlaggingFields[0], 
	  SlopeFlaggingFields[1],
	  SlopeFlaggingFields[2], 
	  SlopeFlaggingFields[3],
	  SlopeFlaggingFields[4],
	  SlopeFlaggingFields[5],
	  SlopeFlaggingFields[6]);

  fprintf(fptr, "MinimumSlopeForRefinement ="
	  " %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
	  MinimumSlopeForRefinement[0],
	  MinimumSlopeForRefinement[1],
	  MinimumSlopeForRefinement[2],
	  MinimumSlopeForRefinement[3],
	  MinimumSlopeForRefinement[4],
	  MinimumSlopeForRefinement[5],
	  MinimumSlopeForRefinement[6]);

  fprintf(fptr, "SecondDerivativeFlaggingFields ="
	  " %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",
	  SecondDerivativeFlaggingFields[0], 
	  SecondDerivativeFlaggingFields[1],
	  SecondDerivativeFlaggingFields[2], 
	  SecondDerivativeFlaggingFields[3],
	  SecondDerivativeFlaggingFields[4],
	  SecondDerivativeFlaggingFields[5],
	  SecondDerivativeFlaggingFields[6]);

  fprintf(fptr, "MinimumSecondDerivativeForRefinement ="
	  " %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
	  MinimumSecondDerivativeForRefinement[0],
	  MinimumSecondDerivativeForRefinement[1],
	  MinimumSecondDerivativeForRefinement[2],
	  MinimumSecondDerivativeForRefinement[3],
	  MinimumSecondDerivativeForRefinement[4],
	  MinimumSecondDerivativeForRefinement[5],
	  MinimumSecondDerivativeForRefinement[6]);

  fprintf(fptr, "SecondDerivativeEpsilon = %"GSYM"\n",
	  SecondDerivativeEpsilon);

  fprintf(fptr, "MinimumOverDensityForRefinement ="
	  " %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n",
	  MinimumOverDensityForRefinement[0],
	  MinimumOverDensityForRefinement[1],
	  MinimumOverDensityForRefinement[2],
	  MinimumOverDensityForRefinement[3],
	  MinimumOverDensityForRefinement[4],
	  MinimumOverDensityForRefinement[5],
	  MinimumOverDensityForRefinement[6]);

  fprintf(fptr, "MinimumMassForRefinement ="
	  " %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM" %.9"GSYM"\n",
	  MinimumMassForRefinement[0]*massu,
	  MinimumMassForRefinement[1]*massu,
	  MinimumMassForRefinement[2]*massu,
	  MinimumMassForRefinement[3]*massu,
	  MinimumMassForRefinement[4]*massu,
	  MinimumMassForRefinement[5]*massu,
	  MinimumMassForRefinement[6]*massu);

  fprintf(fptr, "MinimumMassForRefinementLevelExponent ="
	  " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",
	  MinimumMassForRefinementLevelExponent[0],
	  MinimumMassForRefinementLevelExponent[1],
	  MinimumMassForRefinementLevelExponent[2],
	  MinimumMassForRefinementLevelExponent[3],
	  MinimumMassForRefinementLevelExponent[4],
	  MinimumMassForRefinementLevelExponent[5],
	  MinimumMassForRefinementLevelExponent[6]);

  fprintf(fptr, "MinimumShearForRefinement             = %e\n",
	  MinimumShearForRefinement);
  fprintf(fptr, "OldShearMethod                        = %"ISYM"\n",
	  OldShearMethod);
  fprintf(fptr, "MinimumPressureJumpForRefinement      = %e\n",
	  MinimumPressureJumpForRefinement);
  fprintf(fptr, "MinimumEnergyRatioForRefinement       = %e\n",
	  MinimumEnergyRatioForRefinement);
  fprintf(fptr, "ShockwaveRefinementMinMach            = %"FSYM"\n",
         ShockwaveRefinementMinMach);
  fprintf(fptr, "ShockwaveRefinementMinVelocity        = %"FSYM"\n",
         ShockwaveRefinementMinVelocity);
  fprintf(fptr, "ShockwaveRefinementMaxLevel           = %"ISYM"\n",
         ShockwaveRefinementMaxLevel);
  fprintf(fptr, "ComovingCoordinates                   = %"ISYM"\n",
	  ComovingCoordinates);
  fprintf(fptr, "StarParticleCreation                  = %"ISYM"\n",
	  StarParticleCreation);
  fprintf(fptr, "BigStarFormation                      = %"ISYM"\n",
	  BigStarFormation);
  fprintf(fptr, "BigStarFormationDone                  = %"ISYM"\n",
	  BigStarFormationDone);
  fprintf(fptr, "BigStarSeparation                     = %"FSYM"\n",
	  BigStarSeparation);
  fprintf(fptr, "SimpleQ                               = %lg\n",
	  SimpleQ);
  fprintf(fptr, "SimpleRampTime                        = %"FSYM"\n",
	  SimpleRampTime);
  fprintf(fptr, "StarFormationOncePerRootGridTimeStep  = %"ISYM"\n",
	  StarFormationOncePerRootGridTimeStep);
  fprintf(fptr, "StarParticleFeedback                  = %"ISYM"\n",
	  StarParticleFeedback);
  fprintf(fptr, "NumberOfParticleAttributes            = %"ISYM"\n",
	  NumberOfParticleAttributes);
  fprintf(fptr, "AddParticleAttributes                 = %"ISYM"\n", 
	  AddParticleAttributes);

    /* Sink particles (for present day star formation) & winds */
  fprintf(fptr, "SinkMergeDistance                     = %"FSYM"\n", 
	  SinkMergeDistance);
  fprintf(fptr, "SinkMergeMass                         = %"FSYM"\n", 
	  SinkMergeMass);
  fprintf(fptr, "StellarWindFeedback                   = %"ISYM"\n", 
	  StellarWindFeedback);
  fprintf(fptr, "StellarWindTurnOnMass                 = %"FSYM"\n", 
	  StellarWindTurnOnMass);
  fprintf(fptr, "MSStellarWindTurnOnMass                 = %"FSYM"\n", 
	  MSStellarWindTurnOnMass);


  fprintf(fptr, "StarMakerOverDensityThreshold         = %"GSYM"\n",
	  StarMakerOverDensityThreshold);
  fprintf(fptr, "StarMakerSHDensityThreshold           = %"GSYM"\n",
      StarMakerSHDensityThreshold);
  fprintf(fptr, "StarMakerTimeIndependentFormation     = %"ISYM"\n",
	  StarMakerTimeIndependentFormation);
  fprintf(fptr, "StarMakerMassEfficiency               = %"GSYM"\n",
	  StarMakerMassEfficiency);
  fprintf(fptr, "StarMakerMinimumMass                  = %"GSYM"\n",
          StarMakerMinimumMass);
  fprintf(fptr, "StarMakerMinimumDynamicalTime         = %"GSYM"\n",
          StarMakerMinimumDynamicalTime);
  fprintf(fptr, "StarMassEjectionFraction              = %"GSYM"\n",
          StarMassEjectionFraction);
  fprintf(fptr, "StarMetalYield                        = %"GSYM"\n",
          StarMetalYield);
  fprintf(fptr, "StarEnergyToThermalFeedback           = %"GSYM"\n",
          StarEnergyToThermalFeedback);
  fprintf(fptr, "StarEnergyToStellarUV                 = %"GSYM"\n",
          StarEnergyToStellarUV);
  fprintf(fptr, "StarEnergyToQuasarUV                  = %"GSYM"\n",
          StarEnergyToQuasarUV);
  fprintf(fptr, "StarFeedbackDistRadius                = %"ISYM"\n",
          StarFeedbackDistRadius);
  fprintf(fptr, "StarFeedbackDistCellStep              = %"ISYM"\n",
          StarFeedbackDistCellStep);
  fprintf(fptr, "StarMakerTypeIaSNe                    = %"ISYM"\n",
	  StarMakerTypeIaSNe);
  fprintf(fptr, "StarMakerPlanetaryNebulae             = %"ISYM"\n",
	  StarMakerPlanetaryNebulae);
  fprintf(fptr, "MultiMetals                           = %"ISYM"\n\n",
          MultiMetals);
  fprintf(fptr, "IsotropicConduction                   = %"ISYM"\n", IsotropicConduction);
  fprintf(fptr, "AnisotropicConduction                 = %"ISYM"\n", AnisotropicConduction);
  fprintf(fptr, "IsotropicConductionSpitzerFraction    = %"FSYM"\n", IsotropicConductionSpitzerFraction);
  fprintf(fptr, "AnisotropicConductionSpitzerFraction  = %"FSYM"\n", AnisotropicConductionSpitzerFraction);
  fprintf(fptr, "ConductionCourantSafetyNumber   = %"FSYM"\n", ConductionCourantSafetyNumber);
  fprintf(fptr, "SpeedOfLightTimeStepLimit             = %"ISYM"\n", SpeedOfLightTimeStepLimit);

  fprintf(fptr, "RefineByJeansLengthUnits              = %"ISYM"\n",RefineByJeansLengthUnits);
  fprintf(fptr, "IsothermalSoundSpeed                  = %"GSYM"\n",IsothermalSoundSpeed);
          
  fprintf(fptr, "StarClusterUseMetalField              = %"ISYM"\n",
	  StarClusterUseMetalField);
  fprintf(fptr, "StarClusterUnresolvedModel            = %"ISYM"\n",
	  StarClusterUnresolvedModel);
  fprintf(fptr, "StarClusterMinDynamicalTime           = %"GSYM"\n",
          StarClusterMinDynamicalTime);
  fprintf(fptr, "StarClusterIonizingLuminosity         = %lg\n",
          StarClusterIonizingLuminosity);
  fprintf(fptr, "StarClusterHeliumIonization           = %"ISYM"\n",
	  StarClusterHeliumIonization);
  fprintf(fptr, "StarClusterSNEnergy                   = %lg\n",
          StarClusterSNEnergy);
  fprintf(fptr, "StarClusterSNRadius                   = %"GSYM"\n",
          StarClusterSNRadius);
  fprintf(fptr, "StarClusterFormEfficiency             = %"GSYM"\n",
          StarClusterFormEfficiency);
  fprintf(fptr, "StarClusterMinimumMass                = %"GSYM"\n",
          StarClusterMinimumMass);
  fprintf(fptr, "StarClusterCombineRadius              = %"GSYM"\n",
          StarClusterCombineRadius);
  fprintf(fptr, "StarClusterRegionLeftEdge             = %"FSYM" %"FSYM" %"FSYM"\n",
          StarClusterRegionLeftEdge[0], StarClusterRegionLeftEdge[1],
	  StarClusterRegionLeftEdge[2]);
  fprintf(fptr, "StarClusterRegionRightEdge            = %"FSYM" %"FSYM" %"FSYM"\n",
          StarClusterRegionRightEdge[0], StarClusterRegionRightEdge[1],
	  StarClusterRegionRightEdge[2]);
  fprintf(fptr, "PopIIIStarMass                        = %"GSYM"\n",
          PopIIIStarMass);
  fprintf(fptr, "PopIIIInitialMassFunction             = %"ISYM"\n",
          PopIIIInitialMassFunction);
  fprintf(fptr, "PopIIIInitialMassFunctionSeed         = %"ISYM"\n",
          PopIIIInitialMassFunctionSeed);
  fprintf(fptr, "PopIIIInitialMassFunctionCalls        = %"ISYM"\n",
          PopIIIInitialMassFunctionCalls);
  fprintf(fptr, "PopIIIMassRange                       = %"FSYM" %"FSYM"\n",
          PopIIILowerMassCutoff, PopIIIUpperMassCutoff);
  fprintf(fptr, "PopIIIInitialMassFunctionSlope        = %"FSYM"\n",
          PopIIIInitialMassFunctionSlope);
  fprintf(fptr, "PopIIIHeliumIonization                = %"ISYM"\n",
	  PopIIIHeliumIonization);
  fprintf(fptr, "PopIIIBlackHoles                      = %"ISYM"\n",
          PopIIIBlackHoles);
  fprintf(fptr, "PopIIIBHLuminosityEfficiency          = %"FSYM"\n",
          PopIIIBHLuminosityEfficiency);
  fprintf(fptr, "PopIIIOverDensityThreshold            = %"GSYM"\n",
          PopIIIOverDensityThreshold);
  fprintf(fptr, "PopIIIH2CriticalFraction              = %"GSYM"\n",
          PopIIIH2CriticalFraction);
  fprintf(fptr, "PopIIIMetalCriticalFraction           = %"GSYM"\n",
          PopIIIMetalCriticalFraction);
  fprintf(fptr, "PopIIISupernovaRadius                 = %"GSYM"\n",
          PopIIISupernovaRadius);
  fprintf(fptr, "PopIIISupernovaUseColour              = %"ISYM"\n",
          PopIIISupernovaUseColour);
  fprintf(fptr, "PopIIISupernovaMustRefine             = %"ISYM"\n",
          PopIIISupernovaMustRefine);
  fprintf(fptr, "PopIIISupernovaMustRefineResolution   = %"ISYM"\n\n",
          PopIIISupernovaMustRefineResolution);

  fprintf(fptr, "PopIIIColorDensityThreshold           = %"GSYM"\n",
          PopIIIColorDensityThreshold);
  fprintf(fptr, "PopIIIColorMass                       = %"GSYM"\n\n",
          PopIIIColorMass);

  fprintf(fptr, "MBHAccretion                          = %"ISYM"\n", MBHAccretion);
  fprintf(fptr, "MBHAccretionRadius                    = %"GSYM"\n", MBHAccretionRadius);
  fprintf(fptr, "MBHAccretingMassRatio                 = %"GSYM"\n", MBHAccretingMassRatio);
  fprintf(fptr, "MBHAccretionFixedTemperature          = %"GSYM"\n", MBHAccretionFixedTemperature);
  fprintf(fptr, "MBHAccretionFixedRate                 = %"GSYM"\n", MBHAccretionFixedRate);
  fprintf(fptr, "MBHTurnOffStarFormation               = %"ISYM"\n", MBHTurnOffStarFormation);
  fprintf(fptr, "MBHCombineRadius                      = %"GSYM"\n", MBHCombineRadius);
  fprintf(fptr, "MBHMinDynamicalTime                   = %"GSYM"\n", MBHMinDynamicalTime);
  fprintf(fptr, "MBHMinimumMass                        = %"GSYM"\n\n", MBHMinimumMass);

  fprintf(fptr, "MBHFeedback                           = %"ISYM"\n", MBHFeedback);
  fprintf(fptr, "MBHFeedbackRadiativeEfficiency        = %"GSYM"\n", MBHFeedbackRadiativeEfficiency);
  fprintf(fptr, "MBHFeedbackEnergyCoupling             = %"GSYM"\n", MBHFeedbackEnergyCoupling);
  fprintf(fptr, "MBHFeedbackMassEjectionFraction       = %"GSYM"\n", MBHFeedbackMassEjectionFraction);
  fprintf(fptr, "MBHFeedbackMetalYield                 = %"GSYM"\n", MBHFeedbackMetalYield);
  fprintf(fptr, "MBHFeedbackThermalRadius              = %"GSYM"\n", MBHFeedbackThermalRadius);
  fprintf(fptr, "MBHFeedbackJetsThresholdMass          = %"GSYM"\n\n", MBHFeedbackJetsThresholdMass);

  fprintf(fptr, "MBHParticleIO                         = %"ISYM"\n",
	  MBHParticleIO);
  if (MBHParticleIOFilename != NULL)
    fprintf(fptr, "MBHParticleIOFilename               = %s\n", MBHParticleIOFilename);
  if (MBHInsertLocationFilename != NULL)
    fprintf(fptr, "MBHInsertLocationFilename           = %s\n\n", MBHInsertLocationFilename);

  fprintf(fptr, "ClusterSMBHFeedback           = %"ISYM"\n", ClusterSMBHFeedback);
  fprintf(fptr, "ClusterSMBHJetMdot            = %"FSYM"\n", ClusterSMBHJetMdot);
  fprintf(fptr, "ClusterSMBHJetVelocity        = %"FSYM"\n", ClusterSMBHJetVelocity);
  fprintf(fptr, "ClusterSMBHJetRadius          = %"FSYM"\n", ClusterSMBHJetRadius);
  fprintf(fptr, "ClusterSMBHJetLaunchOffset    = %"FSYM"\n", ClusterSMBHJetLaunchOffset);
  fprintf(fptr, "ClusterSMBHStartTime          = %"FSYM"\n", ClusterSMBHStartTime);
  fprintf(fptr, "ClusterSMBHTramp              = %"FSYM"\n", ClusterSMBHTramp);
  fprintf(fptr, "ClusterSMBHJetOpenAngleRadius = %"FSYM"\n", ClusterSMBHJetOpenAngleRadius);
  fprintf(fptr, "ClusterSMBHFastJetRadius      = %"FSYM"\n", ClusterSMBHFastJetRadius);
  fprintf(fptr, "ClusterSMBHFastJetVelocity    = %"FSYM"\n", ClusterSMBHFastJetVelocity);
  fprintf(fptr, "ClusterSMBHJetEdot            = %"FSYM"\n", ClusterSMBHJetEdot);
  fprintf(fptr, "ClusterSMBHKineticFraction    = %"FSYM"\n", ClusterSMBHKineticFraction);
  fprintf(fptr, "ClusterSMBHJetAngleTheta      = %"FSYM"\n", ClusterSMBHJetAngleTheta);
  fprintf(fptr, "ClusterSMBHJetAnglePhi        = %"FSYM"\n", ClusterSMBHJetAnglePhi);
  fprintf(fptr, "ClusterSMBHJetPrecessionPeriod= %"FSYM"\n", ClusterSMBHJetPrecessionPeriod);
  fprintf(fptr, "ClusterSMBHCalculateGasMass   = %"ISYM"\n", ClusterSMBHCalculateGasMass);
  fprintf(fptr, "ClusterSMBHFeedbackSwitch     = %"ISYM"\n", ClusterSMBHFeedbackSwitch);
  fprintf(fptr, "ClusterSMBHEnoughColdGas      = %"FSYM"\n", ClusterSMBHEnoughColdGas);
  fprintf(fptr, "ClusterSMBHAccretionTime      = %"FSYM"\n", ClusterSMBHAccretionTime);
  fprintf(fptr, "ClusterSMBHJetDim             = %"ISYM"\n", ClusterSMBHJetDim);
  fprintf(fptr, "ClusterSMBHAccretionEpsilon   = %"FSYM"\n", ClusterSMBHAccretionEpsilon);
  fprintf(fptr, "H2StarMakerEfficiency              = %"GSYM"\n", H2StarMakerEfficiency);
  fprintf(fptr, "H2StarMakerNumberDensityThreshold  = %"GSYM"\n", H2StarMakerNumberDensityThreshold);
  fprintf(fptr, "H2StarMakerMinimumMass             = %"GSYM"\n", H2StarMakerMinimumMass);
  fprintf(fptr, "H2StarMakerMinimumH2FractionForStarFormation = %"GSYM"\n", H2StarMakerMinimumH2FractionForStarFormation);
  fprintf(fptr, "H2StarMakerStochastic              = %"ISYM"\n", H2StarMakerStochastic);
  fprintf(fptr, "H2StarMakerUseSobolevColumn        = %"ISYM"\n", H2StarMakerUseSobolevColumn);
  fprintf(fptr, "H2StarMakerSigmaOverR              = %"GSYM"\n", H2StarMakerSigmaOverR);
  fprintf(fptr, "H2StarMakerAssumeColdWarmPressureBalance = %"ISYM"\n", H2StarMakerAssumeColdWarmPressureBalance);
  fprintf(fptr, "H2StarMakerH2DissociationFlux_MW   = %"GSYM"\n", H2StarMakerH2DissociationFlux_MW);
  fprintf(fptr, "H2StarMakerH2FloorInColdGas        = %"GSYM"\n\n", H2StarMakerH2FloorInColdGas);
  fprintf(fptr, "H2StarMakerColdGasTemperature      = %"GSYM"\n\n", H2StarMakerColdGasTemperature);

  /* Most Stanford additions: */

  fprintf(fptr, "UseHydro                   = %"ISYM"\n", UseHydro);
  fprintf(fptr, "Theta_Limiter              = %f\n", Theta_Limiter);
  fprintf(fptr, "RiemannSolver              = %d\n", RiemannSolver);
  fprintf(fptr, "RiemannSolverFallback      = %d\n", RiemannSolverFallback);
  fprintf(fptr, "ConservativeReconstruction = %d\n", ConservativeReconstruction);
  fprintf(fptr, "PositiveReconstruction     = %d\n", PositiveReconstruction);
  fprintf(fptr, "ReconstructionMethod       = %d\n", ReconstructionMethod);
  fprintf(fptr, "RKOrder                    = %d\n", RKOrder);
  fprintf(fptr, "UsePhysicalUnit            = %d\n", UsePhysicalUnit);
  fprintf(fptr, "UseFloor                   = %d\n", UseFloor);
  fprintf(fptr, "UseViscosity               = %d\n", UseViscosity);
  fprintf(fptr, "ViscosityCoefficient       = %g\n", ViscosityCoefficient);  
  fprintf(fptr, "UseAmbipolarDiffusion      = %d\n", UseAmbipolarDiffusion);
  fprintf(fptr, "UseResistivity             = %d\n", UseResistivity);
  fprintf(fptr, "SmallRho                   = %g\n", SmallRho);
  fprintf(fptr, "SmallP                     = %g\n", SmallP);
  fprintf(fptr, "SmallT                     = %g\n", SmallT);
  fprintf(fptr, "MaximumAlvenSpeed          = %g\n", MaximumAlvenSpeed);
  fprintf(fptr, "Coordinate                 = %d\n", Coordinate);
  fprintf(fptr, "EOSType                    = %d\n", EOSType);
  fprintf(fptr, "EOSSoundSpeed              = %g\n", EOSSoundSpeed);
  fprintf(fptr, "EOSCriticalDensity         = %g\n", EOSCriticalDensity);
  fprintf(fptr, "EOSGamma                   = %g\n", EOSGamma); 
  fprintf(fptr, "Mu                         = %g\n", Mu);
  fprintf(fptr, "DivBDampingLength          = %g\n", DivBDampingLength);
  fprintf(fptr, "CoolingCutOffDensity1      = %g\n", CoolingCutOffDensity1);
  fprintf(fptr, "CoolingCutOffDensity2      = %g\n", CoolingCutOffDensity2);
  fprintf(fptr, "CoolingCutOffTemperature   = %g\n", CoolingCutOffTemperature);
  fprintf(fptr, "CoolingPowerCutOffDensity1 = %g\n", CoolingPowerCutOffDensity1);
  fprintf(fptr, "CoolingPowerCutOffDensity2 = %g\n", CoolingPowerCutOffDensity2);
  fprintf(fptr, "UseConstantAcceleration    = %d\n", UseConstantAcceleration);
  fprintf(fptr, "ConstantAcceleration       = %g %g %g\n", ConstantAcceleration[0],
	  ConstantAcceleration[1], ConstantAcceleration[2]);

  fprintf(fptr, "UseDrivingField            = %d\n", UseDrivingField);
  fprintf(fptr, "DrivingEfficiency          = %f\n", DrivingEfficiency);
#ifdef ECUDA
  fprintf(fptr, "UseCUDA = %"ISYM"\n", UseCUDA);
#endif

  /* Poisson Solver */

  fprintf(fptr, "DivergenceCleaningBoundaryBuffer = %"ISYM"\n",
	  DivergenceCleaningBoundaryBuffer);
  fprintf(fptr, "UseDivergenceCleaning            = %d\n", UseDivergenceCleaning);
  fprintf(fptr, "DivergenceCleaningThreshold      = %g\n", 
	  DivergenceCleaningThreshold);
  fprintf(fptr, "PoissonApproximationThreshold    = %g\n", 
	  PoissonApproximationThreshold);
  fprintf(fptr, "PoissonBoundaryType    = %d\n", 
	  PoissonBoundaryType);

  /* Gas Drag */ 
  fprintf(fptr, "UseGasDrag                       = %"ISYM"\n",UseGasDrag);
  fprintf(fptr, "GasDragCoefficient               = %"FSYM"\n",GasDragCoefficient);

  /* Shearing Box Boundary parameters */
  fprintf(fptr, "AngularVelocity              = %"FSYM"\n",AngularVelocity);
  fprintf(fptr, "VelocityGradient             = %"FSYM"\n",VelocityGradient);
  fprintf(fptr, "ShearingVelocityDirection    = %"ISYM"\n",ShearingVelocityDirection);
  fprintf(fptr, "ShearingBoxProblemType    = %"ISYM"\n\n", ShearingBoxProblemType);

  /* write data which defines the boundary conditions */
 
  fprintf(fptr, "LeftFaceBoundaryCondition  = ");
  WriteListOfInts(fptr, MetaData.TopGridRank,
		  (int*) MetaData.LeftFaceBoundaryCondition);
  fprintf(fptr, "RightFaceBoundaryCondition = ");
  WriteListOfInts(fptr, MetaData.TopGridRank,
		  (int*) MetaData.RightFaceBoundaryCondition);
  if (MetaData.BoundaryConditionName)
    fprintf(fptr, "BoundaryConditionName      = %s\n\n",
	    MetaData.BoundaryConditionName);
 

  /* If appropriate, write Cosmology data. */
 
  if (ComovingCoordinates) {
    if (CosmologyWriteParameters(fptr, MetaData.StopTime, MetaData.Time) ==
	FAIL) {
      ENZO_FAIL("Error in CosmologyWriteParameters.\n");
    }
  }
  else {
    if (WriteUnits(fptr) == FAIL) {
      ENZO_FAIL("Error in WriteUnits.\n");
    }
  }

  /* If radiative transfer, write parameters.  If photon test, write
     source data. */

#ifdef TRANSFER
  if (RadiativeTransferWriteParameters(fptr) == FAIL) {
    ENZO_FAIL("Error in RadiativeTransferWriteParameters.\n");
  }

  if (ProblemType == 50)
    if (WritePhotonSources(fptr, MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in WritePhotonSources.\n");
    }
#endif

  if (UsePhysicalUnit) {
    /* Change input physical parameters into code units */

    StarMakerOverDensityThreshold /= rhou;
 
    if (SinkMergeDistance > 1.0)
      SinkMergeDistance /= lenu;
    SmallRho /= rhou;
    SmallP /= presu;
    SmallT /= tempu;
    MaximumAlvenSpeed /= velu;
    EOSSoundSpeed /=  velu;
    MustRefineParticlesMinimumMass /= massu;
    /*
    for (int i = 0; i < MAX_FLAGGING_METHODS; i++) {
      if (MinimumMassForRefinement[i] != FLOAT_UNDEFINED) {
	MinimumMassForRefinement[i] /= massu;
      }
    }
    */

    /*
    if (!ComovingCoordinates && UsePhysicalUnit) {
      for (int i = 0; i < MAX_FLAGGING_METHODS; i++) {
	if (MinimumOverDensityForRefinement[i] != FLOAT_UNDEFINED) 
	  MinimumOverDensityForRefinement[i] /= rhou;
      }
    }
    */

  }

  MustRefineParticlesMinimumMass /= POW(1/(float(MetaData.TopGridDims[0])
				       *POW(float(RefineBy), float(MustRefineParticlesRefineToLevel))),3);
  //MHDCT variables
  fprintf(fptr, "MHDCTSlopeLimiter          = %"ISYM"\n", MHDCTSlopeLimiter);
  fprintf(fptr, "MHDCTDualEnergyMethod          = %"ISYM"\n", MHDCTDualEnergyMethod);
  fprintf(fptr, "MHDPowellSource          = %"ISYM"\n", MHDCTPowellSource);
  fprintf(fptr, "MHDCTUseSpecificEnergy          = %"ISYM"\n", MHDCTUseSpecificEnergy);
  fprintf(fptr, "WriteBoundary          = %"ISYM"\n", WriteBoundary);
  fprintf(fptr,"CT_AthenaDissipation          =%"GSYM"\n",CT_AthenaDissipation);
  fprintf(fptr,"MHD_WriteElectric             =%"ISYM"\n",MHD_WriteElectric);
  fprintf(fptr,"tiny_pressure                 =%"GSYM"\n",tiny_pressure);
  fprintf(fptr,"MHD_CT_Method                 =%"ISYM"\n",MHD_CT_Method);
  fprintf(fptr,"NumberOfGhostZones           =%"ISYM"\n",NumberOfGhostZones);
  fprintf(fptr,"IsothermalSoundSpeed          =%"GSYM"\n",IsothermalSoundSpeed);
  fprintf(fptr,"FixedTimestep          =%"GSYM"\n",FixedTimestep);
  fprintf(fptr,"MHD_ProjectB                  =%"ISYM"\n",MHD_ProjectB);
  fprintf(fptr,"MHD_ProjectE                  =%"ISYM"\n",MHD_ProjectE);
  fprintf(fptr,"ProcessorTopology             =%"ISYM" %"ISYM" %"ISYM"\n",
          ProcessorTopology[0],ProcessorTopology[1],ProcessorTopology[2]);
  fprintf(fptr,"EquationOfState               =%"ISYM"\n",EquationOfState);

  fprintf(fptr, "CorrectParentBoundaryFlux          = %d\n", CorrectParentBoundaryFlux);

  /* Output current time */
  time_t ID;
  ID = time(NULL);
  fprintf(fptr, "CurrentTimeIdentifier = %"ISYM"\n", int(ID));

  /* If the simulation was given a name, write that. */
  if(MetaData.MetaDataIdentifier != NULL){
    fprintf(fptr, "MetaDataIdentifier              = %s\n",
	    MetaData.MetaDataIdentifier);
  }
  /* Write unique simulation identifier. */
  fprintf(fptr, "MetaDataSimulationUUID          = %s\n", MetaData.SimulationUUID);
  /* Give this dataset a unique identifier. */
  char dset_uuid[MAX_LINE_LENGTH];
  get_uuid(dset_uuid);
  fprintf(fptr, "MetaDataDatasetUUID             = %s\n", dset_uuid);
  /* If the restart data had a UUID, write that. */
  if(MetaData.RestartDatasetUUID != NULL){
    fprintf(fptr, "MetaDataRestartDatasetUUID      = %s\n",
	    MetaData.RestartDatasetUUID);
  }
  if(MetaData.InitialConditionsUUID != NULL){

    fprintf(fptr, "MetaDataInitialConditionsUUID   = %s\n",
	    MetaData.InitialConditionsUUID);
  }

  /* write version info */
 
  fprintf(fptr, "VersionNumber              = %"FSYM"\n\n", VERSION);

  if (name != NULL)
    UpdateLocalDatabase(MetaData, ID, dset_uuid, name);
 
  return SUCCESS;
}
