/***********************************************************************
/
/  SETS THE DEFAULT GLOBAL VALUES
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       February 29th, 2008
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//

#include "preincludes.h" 
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "TopGridData.h"
#include "StarParticleData.h"
#include "phys_constants.h"
 
/* character strings */
 
char DefaultDimUnits[] = "cm";
char *DefaultDimLabel[] = {"x", "y", "z"};
 
char DefaultRestartName[] = "restart";
char DefaultDataName[] = "data";
char DefaultHistoryName[] = "history";
char DefaultRedshiftName[] = "RedshiftOutput";
char DefaultNewMovieName[] = "MoviePack";
char DefaultTracerParticleName[] = "TracerOutput";
 
char DefaultRestartDir[] = "RS";
char DefaultDataDir[] = "DD";
char DefaultHistoryDir[] = "HD";
char DefaultRedshiftDir[] = "RD";
char DefaultTracerParticleDir[] = "TD";
char DefaultExtraName[] = "ExtraDumpXX";
char DefaultExtraDir[]="ED00";
 
 
int SetDefaultGlobalValues(TopGridData &MetaData)
{
 
  /* declarations */
 
  int dim, i, j;
 
  huge_number               = 1.0e+20;
  tiny_number               = 1.0e-20;

  /* set the default MetaData values. */
 
  MetaData.CycleNumber     = 0;
  MetaData.SubcycleNumber     = 0;
  MetaData.Time            = 0.0;
  MetaData.CPUTime         = 0.0;
 
  MetaData.StopTime        = FLOAT_UNDEFINED;  // This must be set be the user
  MetaData.StopCycle       = 100000;            // 10000 timesteps
  MetaData.StopSteps       = 10000;            // 10000 timesteps
  MetaData.StopCPUTime     = 720.0*3600.0;     // 30 days
  MetaData.ResubmitOn      = FALSE;
  MetaData.ResubmitCommand = NULL;
 
  MetaData.MaximumTopGridTimeStep = huge_number;

  MetaData.TimeLastRestartDump = 0.0;
  MetaData.dtRestartDump       = FLOAT_UNDEFINED;
  MetaData.TimeLastDataDump    = FLOAT_UNDEFINED;
  MetaData.dtDataDump          = 0.0;
  MetaData.TimeLastHistoryDump = FLOAT_UNDEFINED;
  MetaData.dtHistoryDump       = 0.0;
  MetaData.TimeLastTracerParticleDump = FLOAT_UNDEFINED;
  MetaData.dtTracerParticleDump       = 0.0;
  MetaData.TimeLastInterpolatedDataDump    = FLOAT_UNDEFINED;
  MetaData.dtInterpolatedDataDump          = 0.0;
  MetaData.WroteData           = FALSE;
 
  MetaData.CycleLastRestartDump = 0;
  MetaData.CycleSkipRestartDump = 0;
  MetaData.CycleLastDataDump    = INT_UNDEFINED;
  MetaData.CycleSkipDataDump    = 0;
  MetaData.SubcycleLastDataDump    = 0;
  MetaData.SubcycleSkipDataDump    = 0;
  MetaData.CycleLastHistoryDump = INT_UNDEFINED;
  MetaData.CycleSkipHistoryDump = 0;
  MetaData.CycleSkipGlobalDataDump = 0; //AK
 
  MetaData.OutputFirstTimeAtLevel = 0; // zero is off
  MetaData.StopFirstTimeAtLevel   = 0; // zero is off
 
  MetaData.NumberOfOutputsBeforeExit = 0;
  MetaData.OutputsLeftBeforeExit     = 0;

  MetaData.RestartDumpNumber   = 0;            // starting restart id number
  MetaData.RestartDumpName     = DefaultRestartName;
  MetaData.RestartDumpDir      = DefaultRestartDir;
  MetaData.DataDumpNumber      = 0;
  MetaData.DataDumpName        = DefaultDataName;
  MetaData.DataDumpDir         = DefaultDataDir;
  MetaData.HistoryDumpNumber   = 0;
  MetaData.HistoryDumpName     = DefaultHistoryName;
  MetaData.HistoryDumpDir      = DefaultHistoryDir;
  MetaData.TracerParticleDumpNumber = 0;
  MetaData.TracerParticleDumpName   = DefaultTracerParticleName;
  MetaData.TracerParticleDumpDir    = DefaultTracerParticleDir;
  //MetaData.RedshiftDumpNumber  = 0;
  MetaData.RedshiftDumpName    = DefaultRedshiftName;
  MetaData.RedshiftDumpDir     = DefaultRedshiftDir;
  MetaData.ExtraDumpDir        = DefaultExtraDir;
  MetaData.ExtraDumpName        = DefaultExtraName;

  MetaData.MetaDataIdentifier    = NULL;
  MetaData.SimulationUUID        = NULL;
  MetaData.RestartDatasetUUID    = NULL;
  MetaData.InitialConditionsUUID = NULL;

  MetaData.LocalDir            = NULL;
  MetaData.GlobalDir           = NULL;

  NumberOfGhostZones = 3;
  LoadBalancing = 1;     //On, memory equalization method
  LoadBalancingCycleSkip = 10;  // Load balance root grids every 10 cycles
  ResetLoadBalancing = FALSE;
  CoresPerNode = 1;
  PreviousMaxTask = 0;
  LoadBalancingMinLevel = 0;     //All Levels
  LoadBalancingMaxLevel = MAX_DEPTH_OF_HIERARCHY;  //All Levels

  FileDirectedOutput = 1;

  // Default Hierarchy File IO settings (1 = ASCII; 2 = HDF5+ASCII)
  HierarchyFileInputFormat = 1;
  HierarchyFileOutputFormat = 1;

  ConductionDynamicRebuildHierarchy = FALSE;
  ConductionDynamicRebuildMinLevel = 0;
  for (i = 0;i < MAX_DEPTH_OF_HIERARCHY;i++) {
    RebuildHierarchyCycleSkip[i] = 1;
  }

  for (i = 0; i < MAX_TIME_ACTIONS; i++) {
    TimeActionType[i]      = 0;
    TimeActionRedshift[i]  = -1;
    TimeActionTime[i]      = 0;
    TimeActionParameter[i] = FLOAT_UNDEFINED;
  }
 
  for (i = 0; i < MAX_CUBE_DUMPS; i++) {
    CubeDumps[i] = NULL;
  }
 
  MetaData.StaticHierarchy     = TRUE;
  FastSiblingLocatorEntireDomain = TRUE;
 
  MetaData.TopGridRank = INT_UNDEFINED;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    MetaData.TopGridDims[dim]                = INT_UNDEFINED;
    MetaData.LeftFaceBoundaryCondition[dim]  = reflecting;
    MetaData.RightFaceBoundaryCondition[dim] = reflecting;
  }
  MetaData.BoundaryConditionName = NULL;
 
  MetaData.GravityBoundary        = TopGridPeriodic;

#ifdef TRANSFER
  MetaData.RadHydroParameterFname = NULL;
#endif
 
  MetaData.ParticleBoundaryType   = periodic;  // only one implemented!
  MetaData.NumberOfParticles      = 0;         // no particles
 
  MetaData.CourantSafetyNumber    = 0.6;
  MetaData.PPMFlatteningParameter = 0;    // off
  MetaData.PPMDiffusionParameter  = 0;    // off
  MetaData.PPMSteepeningParameter = 0;    // off

  MetaData.FirstTimestepAfterRestart = TRUE;
 
  /* set the default global data. */
  CheckpointRestart         = 0;

  // Debug flag set in main
  ProblemType               = 0;                 // None
  HydroMethod               = PPM_DirectEuler;   //
  Gamma                     = 5.0/3.0;           // 5/3
  PressureFree              = FALSE;             // use pressure (duh)
  QuantumPressure           = FALSE;             // FDM: no fuzzy dark matter
  FDMMass           = 1.0;             // FDM: FDM mass 1e-22 eV

  RefineBy                  = 2;                 // Refinement factor
  MaximumRefinementLevel    = 2;                 // three levels (w/ topgrid)
  MaximumGravityRefinementLevel = INT_UNDEFINED;
  MaximumParticleRefinementLevel = -1;            // unused if negative
  MaximumRefinementLevelPhysicalScale = -1.0;     // unused if negative - move max refinement level to match this physical scale (in pc)
  MustRefineRegionMinRefinementLevel = -1;        // unused if negative
  MetallicityRefinementMinLevel = -1;
  MetallicityRefinementMinMetallicity = 1.0e-5;
  MetallicityRefinementMinDensity = FLOAT_UNDEFINED;
  FluxCorrection            = TRUE;

  UseCoolingTimestep = FALSE;
  CoolingTimestepSafetyFactor = 0.1;

  InterpolationMethod       = SecondOrderA;      // ?
  ConservativeInterpolation = TRUE;              // true for ppm
  MinimumEfficiency         = 0.2;               // between 0-1, usually ~0.1
  MinimumSubgridEdge        = 6;                 // min for acceptable subgrid
  MaximumSubgridSize        = 32768;             // max for acceptable subgrid
  CriticalGridRatio         = 3.0;              // max grid ratio

  SubgridSizeAutoAdjust     = TRUE; // true for adjusting maxsize and minedge
  OptimalSubgridsPerProcessor = 16;    // Subgrids per processor
  NumberOfBufferZones       = 1;
 
  for (i = 0; i < MAX_FLAGGING_METHODS; i++) {
    MinimumSlopeForRefinement[i]= 0.3;
    SlopeFlaggingFields[i] = INT_UNDEFINED;
    CellFlaggingMethod[i]       = INT_UNDEFINED;
    MinimumMassForRefinement[i] = FLOAT_UNDEFINED;   // usually set by:
    MinimumOverDensityForRefinement[i]       = 1.5;
    MinimumMassForRefinementLevelExponent[i] = 0;
    MinimumSecondDerivativeForRefinement[i]= 0.3;
    SecondDerivativeFlaggingFields[i] = INT_UNDEFINED;
  }
  SecondDerivativeEpsilon = 1.0e-2;
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    DomainLeftEdge[dim]             = 0.0;
    DomainRightEdge[dim]            = 1.0;
    GridVelocity[dim]               = 0.0;
    DimLabels[dim]                  = DefaultDimLabel[dim];
    DimUnits[dim]                   = DefaultDimUnits;
    RefineRegionLeftEdge[dim]       = FLOAT_UNDEFINED;
    RefineRegionRightEdge[dim]      = FLOAT_UNDEFINED;
    RefineRegionAutoAdjust          = FALSE;
    MetaData.NewMovieLeftEdge[dim]  = 0.0;
    MetaData.NewMovieRightEdge[dim] = 1.0;
    PointSourceGravityPosition[dim] = 0.0;
    MustRefineRegionLeftEdge[dim]   = 0.0;
    MustRefineRegionRightEdge[dim]  = 1.0;
    MustRefineParticlesLeftEdge[dim] = 0.0;
    MustRefineParticlesRightEdge[dim] = 0.0;
    CoolingRefineRegionLeftEdge[dim]   = 0.0;
    CoolingRefineRegionRightEdge[dim]  = 1.0;
    DiskGravityPosition[dim]        = 0.0;
    DiskGravityAngularMomentum[dim] = 0.0;
    DiskGravityDarkMatterCOM[dim]   = 0.0;
    GalaxySimulationRPSWindVelocity[dim] = 0.0;
    GalaxySimulationPreWindVelocity[dim] = 0.0;
    StellarWindCenterPosition[dim] = 0.5;
  }
  if( MAX_DIMENSION > 0 ) DiskGravityAngularMomentum[MAX_DIMENSION-1] = 1.0; 

  MultiRefineRegionMaximumOuterLevel = INT_UNDEFINED;
  MultiRefineRegionMinimumOuterLevel = INT_UNDEFINED;
  for (i = 0; i < MAX_STATIC_REGIONS; i++) {
    MultiRefineRegionMaximumLevel[i] = INT_UNDEFINED;
    MultiRefineRegionMinimumLevel[i] = 0;
    MultiRefineRegionGeometry[i] = -1; 
    MultiRefineRegionRadius[i] = INT_UNDEFINED;
    MultiRefineRegionWidth[i] = 3.0;
    MultiRefineRegionStaggeredRefinement[i] = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      MultiRefineRegionLeftEdge[i][dim] = FLOAT_UNDEFINED;
      MultiRefineRegionRightEdge[i][dim] = FLOAT_UNDEFINED;
      MultiRefineRegionCenter[i][dim]         = FLOAT_UNDEFINED;
      MultiRefineRegionOrientation[i][dim]    = FLOAT_UNDEFINED;
    }
  }

  for (i = 0; i < MAX_STATIC_REGIONS; i++) {
    StaticRefineRegionLevel[i] = INT_UNDEFINED;
    AvoidRefineRegionLevel[i]  = INT_UNDEFINED;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      AvoidRefineRegionLeftEdge[i][dim] = FLOAT_UNDEFINED;
      AvoidRefineRegionRightEdge[i][dim] = FLOAT_UNDEFINED;
    }
  }

  /* For evolving refinement regions. */
  RefineRegionFile = NULL;
  RefineRegionTimeType = -1; /* 0=time bins 1=redshift bins*/
  for (i = 0; i < MAX_REFINE_REGIONS; i++) {
    EvolveRefineRegionTime[i] = FLOAT_UNDEFINED;
    for (j = 0; j < MAX_DIMENSION; j++) {
      EvolveRefineRegionLeftEdge[i][j]  = FLOAT_UNDEFINED;
      EvolveRefineRegionRightEdge[i][j] = FLOAT_UNDEFINED;
    }
  }

  /* for cooling refinement regions. */

  UseCoolingRefineRegion = FALSE;
  EvolveCoolingRefineRegion = FALSE;
  CoolingRefineRegionFile = NULL;
  CoolingRefineRegionTimeType = -1; /* 0=time bins 1=redshift bins*/

  for (i = 0; i < MAX_REFINE_REGIONS; i++) {
    EvolveCoolingRefineRegionTime[i] = FLOAT_UNDEFINED;
    for (j = 0; j < MAX_DIMENSION; j++) {
      EvolveCoolingRefineRegionLeftEdge[i][j]  = FLOAT_UNDEFINED;
      EvolveCoolingRefineRegionRightEdge[i][j] = FLOAT_UNDEFINED;
    }
  }


  /* For evolving MustRefine regions. */
  MustRefineRegionFile = NULL;
  MustRefineRegionTimeType = -1; /* 0=time bins 1=redshift bins*/
  for (i = 0; i < MAX_REFINE_REGIONS; i++) {
    EvolveMustRefineRegionTime[i] = FLOAT_UNDEFINED;
    EvolveMustRefineRegionMinLevel[i] = INT_UNDEFINED;
    for (j = 0; j < MAX_DIMENSION; j++) {
      EvolveMustRefineRegionLeftEdge[i][j]  = FLOAT_UNDEFINED;
      EvolveMustRefineRegionRightEdge[i][j] = FLOAT_UNDEFINED;
    }
  }

  /* For storing mass of density fields living domain */
  StoreDomainBoundaryMassFlux = 0;
  BoundaryMassFluxFilename    = "boundary_mass_flux.dat";
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++){
    BoundaryMassFluxFieldNumbers[i] = -1;
    BoundaryMassFluxContainer[i] = 0.0;
  }

  DatabaseLocation = NULL;

  ParallelRootGridIO          = FALSE;
  ParallelParticleIO          = FALSE;
  Unigrid                     = FALSE;
  UnigridTranspose            = 2;
  NumberOfRootGridTilesPerDimensionPerProcessor = 1;
  PartitionNestedGrids        = FALSE;
  ExtractFieldsOnly           = TRUE;
  for (i = 0; i < MAX_DIMENSION; i++) {
    UserDefinedRootGridLayout[i] = INT_UNDEFINED;
  }

  ExternalBoundaryIO          = FALSE;
  ExternalBoundaryTypeIO      = FALSE;
  ExternalBoundaryValueIO     = FALSE;
  SimpleConstantBoundary      = FALSE;

  debug1                      = 0;
  debug2                      = 0;

  TracerParticleOn            = FALSE;
  TracerParticleOutputVelocity = FALSE;

  OutputOnDensity                  = 0;
  StartDensityOutputs              = 999;
  CurrentDensityOutput             = 999;
  IncrementDensityOutput           = 999;
  CurrentMaximumDensity            = -999;
  StopFirstTimeAtDensity              = 0.0;
  StopFirstTimeAtMetalEnrichedDensity = 0.0;
  CurrentMaximumMetalEnrichedDensity  = -999;
  EnrichedMetalFraction               = 1.e-8;

  CubeDumpEnabled             = 0;

#ifdef STAGE_INPUT
  StageInput                  = 0;
#endif

  First_Pass                  = 0;

  MemoryLimit                 = 4000000000L;
 
  ExternalGravity             = FALSE;             // off
  ExternalGravityDensity      = 0.0;
  ExternalGravityRadius       = 0.0;

  // for externalgravity type 2 - 4
  ExternalGravityTimeOn       = 0.0;               // time in code units to turn on time varying external potential
  ExternalGravityTimeOff      = huge_number;       // time in code units to turn off time varying external potential
  ExternalGravityTime         = NULL;
  for (dim = 0 ; dim < MAX_DIMENSION; dim ++){
    ExternalGravityTimePositions[dim] = NULL;
  }
  ExternalGravityMass         = 0.0; // mass in solar masses of time varying point source potential (EG = 4)


  UniformGravity              = FALSE;             // off
  UniformGravityDirection     = 0;                 // x-direction
  UniformGravityConstant      = 1.0;
 
  PointSourceGravity           = FALSE;             // off
  PointSourceGravityConstant   = 1.0;
  PointSourceGravityCoreRadius = 0.0;

  DiskGravity                        = FALSE;
  DiskGravityStellarDiskMass         = 1.0E11;      // Solar Masses
  DiskGravityStellarDiskScaleHeightR = 4.0E-3;      // Mpc
  DiskGravityStellarDiskScaleHeightz = 2.5E-4;      // Mpc
  DiskGravityStellarBulgeMass        = 1.0E10;      // Solar Masses
  DiskGravityStellarBulgeR           = 4.0E-4;      // Mpc
  DiskGravityDarkMatterR             = 2.3E-2;      // Mpc
  DiskGravityDarkMatterDensity       = 3.81323E-25; // CGS
  DiskGravityDarkMatterMassInterior  = 0.0;         // solar masses - Only used when above is -1
  DiskGravityDarkMatterMassInteriorR = 0.0;         // Mpc - Only used in conjuction with above
  DiskGravityDarkMatterUpdateCOM     = 0;           // when live DM is used, compute COM every root time step
  DiskGravityDarkMatterRefineCore    = 0;           // if > 0, update must refine region to be around halo center with this being 1/2 side length in code units

  DiskGravityDoublePower             = FALSE;       // use general doubnle power law DM dis - use with particles
  DiskGravityDoublePowerR            = NULL;
  DiskGravityDoublePowerMass         = NULL;
  DiskGravityDoublePowerPot          = NULL;
  DiskGravityDarkMatterCutoffR       = huge_number;
  DiskGravityDarkMatterAlpha         = 1.0;
  DiskGravityDarkMatterBeta          = 3.0;
  DiskGravityDarkMatterGamma         = 1.0;
  DiskGravityDarkMatterDelta         = -1.0;
  DiskGravityDarkMatterRDecay        = huge_number;

  SelfGravity                 = FALSE;             // off
  SelfGravityGasOff           = FALSE;             // off
  AccretionKernal             = FALSE;             // off
  CopyGravPotential           = FALSE;             // off
  PotentialIterations         = 4;                 // ~4 is reasonable
  GravitationalConstant       = 4*pi;              // G = 1
  ComputePotential            = FALSE;
  WritePotential              = FALSE;
  ParticleSubgridDepositMode  = CIC_DEPOSIT_SMALL;

  GalaxySimulationRPSWind = 0;
  GalaxySimulationRPSWindShockSpeed = 0.0;
  GalaxySimulationRPSWindDelay = 0.0;

  GalaxySimulationRPSWindDensity = 1.0;
  GalaxySimulationRPSWindTotalEnergy = 1.0;
  GalaxySimulationRPSWindPressure = 1.0;

  GalaxySimulationPreWindDensity = 1.0;
  GalaxySimulationPreWindTotalEnergy = 1.0;

  GalaxySimulationInitialStellarDist = 0;          // only works in individual star maker - initialize a stellar distribution if ON

  DualEnergyFormalism         = FALSE;             // off
  DualEnergyFormalismEta1     = 0.001;             // typical 0.001
  DualEnergyFormalismEta2     = 0.1;               // 0.08-0.1
  ParticleCourantSafetyNumber = 0.5;
  RootGridCourantSafetyNumber = 1.0;
  RandomForcing               = FALSE;             // off //AK
  RandomForcingEdot           = -1.0;              //AK
  RandomForcingMachNumber     = 0.0;               //AK

  DrivenFlowProfile = 0; // off
  DrivenFlowSeed = 0;
  DrivenFlowWeight = 0.0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    DrivenFlowAlpha[dim] = 0;
    DrivenFlowBandWidth[dim] = 0.0;
    DrivenFlowAutoCorrl[dim] = 0.0;
    DrivenFlowVelocity[dim] = 0.0;
    DrivenFlowDomainLength[dim] = 0.0;
  }

  UseSGSModel = 0; // off
  SGSFilterStencil = 0; // the one-dimensional stencil of the complete filter 
  SGSNeedJacobians = 0; // set automatically in ReadParameter file 
  SGSNeedMixedFilteredQuantities = 0; // set automatically in ReadParameter file
  SGSFilterWidth = 0.; // off, i.e. use grid-scale quantities
  for (i = 0; i < 4; i++)
    // discrete filter weights of explicit filter
    SGSFilterWeights[i] = 0.;
  SGScoeffERS2M2Star = 0.0; // off
  SGScoeffEVStarEnS2Star = 0.0; // off
  SGScoeffEnS2StarTrace = 0.0; // off
  SGScoeffNLemfCompr = 0.0; // off
  SGScoeffNLu = 0.0; // off
  SGScoeffNLuNormedEnS2Star = 0.0; // off
  SGScoeffNLb = 0.0; // off

  RadiativeCooling            = FALSE;             // off
  RadiativeCoolingModel       = 1;                 //1=cool_rates.in table lookup
                                                   //3=Koyama&Inutsuka 2002
  GadgetEquilibriumCooling    = FALSE;             // off
  RadiativeTransfer           = 0;                 // off
  RadiativeTransferFLD        = 0;                 // off
  ImplicitProblem             = 0;                 // off
  StarMakerEmissivityField    = 0;                 // off
  uv_param                    = 1.1e-5;            // consistent with Razoumov Norman 2002

  MultiSpecies                = FALSE;             // off
  NoMultiSpeciesButColors     = FALSE;             // off
  ThreeBodyRate               = 0;                 // ABN02
  CIECooling                  = 1;
  H2OpticalDepthApproximation = 1;
  H2FormationOnDust           = FALSE;
  GloverChemistryModel        = 0;                 // 0ff
  CRModel                     = 0;                 // off
  CRDiffusion                 = 0;                 // off
  CRkappa                     = 0.0;
  CRCourantSafetyNumber       = 0.5;
  CRFeedback                  = 0.0;               // no stellar feedback into CRs
  CRdensFloor                 = 0.0;               // off
  CRmaxSoundSpeed             = 0.0;               // off 
  CRgamma                     = 4.0/3.0;           // relativistic, adiabatic gas
  CosmologySimulationUniformCR= 1e-20;             // FIXME
  ShockMethod                 = 0;                 // off
  ShockTemperatureFloor       = 1.0;               // Set to 1K
  StorePreShockFields         = 0;
  FindShocksOnlyOnOutput      = 0;                 // Find at every cycle and 
                                                   // during output by default.
  RadiationFieldType          = 0;
  RadiationFieldRedshift      = FLOAT_UNDEFINED;
  TabulatedLWBackground       = 0;
  RadiationFieldLevelRecompute = 0;
  RadiationData.RadiationShield = 0;
  AdjustUVBackground          = 1;
  AdjustUVBackgroundHighRedshift = 0;
  SetUVBAmplitude             = 1.0;
  SetHeIIHeatingScale         = 1.8;
  PhotoelectricHeating	      = 0;
  PhotoelectricHeatingRate    = 8.5e-26;           // ergs cm-3 s-1
  RadiationXRaySecondaryIon   = 0;
  RadiationXRayComptonHeating = 0;

  CoolData.alpha0             = 1.5;               // radiation spectral slope
  CoolData.f3                 = 1.0e-21;           // radiation normalization
  CoolData.f0to3                    = 0.1;
  CoolData.RadiationRedshiftOn      = 7.0;
  CoolData.RadiationRedshiftOff     = 0.0;
  CoolData.RadiationRedshiftFullOn  = 6.0;
  CoolData.RadiationRedshiftDropOff = 0.0;
  CoolData.HydrogenFractionByMass   = 0.76;
  /* The DToHRatio is by mass in the code, so multiply by 2. */
  CoolData.DeuteriumToHydrogenRatio = 2.0*3.4e-5; // Burles & Tytler 1998

  /*
     Previously, the solar metal mass fraction was 0.02041.  
     This is close to 0.0194 of Anders & Grevesse (1989), but significantly 
     higher than the more recent value of 0.0122 from Asplund et al. (2005).
     Now, the solar metal mass fraction has been set to 0.01295, 
     which is consistent with the abundances used in Cloudy when generating the 
     Grackle cooling tables.
  */
  CoolData.SolarMetalFractionByMass = 0.01295; // Cloudy v13 abundances

  CoolData.NumberOfTemperatureBins = 600;
  CoolData.ih2co                   = 1;
  CoolData.ipiht                   = 1;
  CoolData.TemperatureStart        = 1.0;
  CoolData.TemperatureEnd          = 1.0e8;
  CoolData.comp_xray               = 0;
  CoolData.temp_xray               = 0;
  RateData.CaseBRecombination      = 0;   // default to case A rates
  RateData.NumberOfDustTemperatureBins = 250;
  RateData.DustTemperatureStart    = 1.0;
  RateData.DustTemperatureEnd      = 1500.0;

  CloudyCoolingData.CloudyCoolingGridRank          = 0;
  CloudyCoolingData.CloudyCoolingGridFile          = "";
  CloudyCoolingData.IncludeCloudyHeating           = 0;
  CloudyCoolingData.CMBTemperatureFloor            = 1;         // use CMB floor.
  CloudyCoolingData.CloudyElectronFractionFactor = 9.153959e-3; // calculated using Cloudy 07.02 abundances

  use_grackle = FALSE;

  OutputCoolingTime = FALSE;
  OutputTemperature = FALSE;
  OutputDustTemperature = FALSE;

  OutputSmoothedDarkMatter = FALSE;
  SmoothedDarkMatterNeighbors = 32;

  OutputGriddedStarParticle = FALSE;

  ZEUSLinearArtificialViscosity    = 0.0;
  ZEUSQuadraticArtificialViscosity = 2.0;
  UseMinimumPressureSupport        = FALSE;
  MinimumPressureSupportParameter  = 100.0;
 
  //MinimumSlopeForRefinement        = 0.3;          // 30% change in value
  MinimumShearForRefinement        = 1.0;          //AK
  OldShearMethod                   = 0;            
  MinimumPressureJumpForRefinement = 0.33;         // As in PPM method paper
  MinimumEnergyRatioForRefinement  = 0.1;          // conservative!
  RefineByJeansLengthSafetyFactor  = 4.0;
  JeansRefinementColdTemperature  = -1.0;
  RefineByResistiveLengthSafetyFactor  = 2.0;
  ShockwaveRefinementMinMach = 1.3; // Only above M=1.3
  ShockwaveRefinementMinVelocity = 1.0e7; //1000 km/s
  ShockwaveRefinementMaxLevel = 0; 
  MustRefineParticlesRefineToLevel = 0;
  MustRefineParticlesCreateParticles = 0;
  MustRefineParticlesRefineToLevelAutoAdjust = FALSE;
  MustRefineParticlesMinimumMass   = 0.0;
  MustRefineParticlesBufferSize    = 8;
  ComovingCoordinates              = FALSE;        // No comoving coordinates
  StarParticleCreation             = FALSE;
  StarParticleFeedback             = FALSE;
  StarParticleRadiativeFeedback    = FALSE;
  BigStarFormation                 = FALSE;
  BigStarFormationDone             = FALSE;
  BigStarSeparation                = 0.25;
  SimpleQ                          = 1e50;
  SimpleRampTime                   = 0.1;
  StarFormationOncePerRootGridTimeStep = FALSE;
  StarMakerTypeIaSNe               = FALSE;
  StarMakerTypeIISNeMetalField     = FALSE;
  StarMakerPlanetaryNebulae        = FALSE;
  StarMakerUseOverDensityThreshold = TRUE;
  StarMakerMaximumFractionCell     = 0.5;
  StarMakerOverDensityThreshold    = 100;          // times mean total density
  StarMakerSHDensityThreshold      = 7e-26;        // cgs density for rho_crit in Springel & Hernquist star_maker5
  StarMakerTimeIndependentFormation = FALSE;
  StarMakerMassEfficiency          = 1;
  StarMakerMinimumMass             = 1.0e9;        // in solar masses
  StarMakerMinimumDynamicalTime    = 1.0e6;        // in years
  StarMassEjectionFraction         = 0.25;
  StarMetalYield                   = 0.02;
  StarEnergyToThermalFeedback      = 1.0e-5;
  StarEnergyToStellarUV            = 3.0e-6;
  StarEnergyToQuasarUV             = 5.0e-6;
  StarFeedbackKineticFraction      = 0.0;
  StarMakerExplosionDelayTime      = -1.0; //in Myrs if >= 0.0
  StarFeedbackDistRadius           = 0;
  StarFeedbackDistCellStep         = 0;
  StarFeedbackDistTotalCells       = 1;
  StarFeedbackPreSN                = TRUE;
  StarMakerUseJeansMass            = TRUE;
  MultiMetals                      = FALSE;
  NumberOfParticleAttributes       = INT_UNDEFINED;
  NumberOfParticleTableIDs         = MAX_NUMBER_OF_PARTICLE_TABLE_POSITIONS;
  ParticleAttributeTableStartIndex = INT_UNDEFINED;
  ParticleTypeInFile               = TRUE;
  ReadGhostZones                   = FALSE;
  WriteGhostZones                  = FALSE;
  OutputParticleTypeGrouping       = FALSE;

  IsotropicConduction = FALSE;
  AnisotropicConduction = FALSE;
  IsotropicConductionSpitzerFraction = 0.0;
  AnisotropicConductionSpitzerFraction = 0.0;
  ConductionCourantSafetyNumber = 0.5;
  SpeedOfLightTimeStepLimit = FALSE;

  ClusterSMBHFeedback              = FALSE;
  ClusterSMBHJetMdot               = 3.0;
  ClusterSMBHJetVelocity           = 10000.0;
  ClusterSMBHJetRadius             = 6.0;
  ClusterSMBHJetLaunchOffset       = 10.0;
  ClusterSMBHStartTime             = 1.0;
  ClusterSMBHTramp                 = 0.1;
  ClusterSMBHJetOpenAngleRadius    = 0.0;
  ClusterSMBHFastJetRadius         = 0.1;
  ClusterSMBHFastJetVelocity       = 10000.0;
  ClusterSMBHJetEdot               = 1.0;
  ClusterSMBHKineticFraction       = 1.0;
  ClusterSMBHJetAngleTheta         = 0.0;
  ClusterSMBHJetAnglePhi           = 0.0;
  ClusterSMBHJetPrecessionPeriod   = 0.0;
  ClusterSMBHCalculateGasMass      = 1;
  ClusterSMBHFeedbackSwitch        = FALSE;
  ClusterSMBHEnoughColdGas         = 1.0e7;
  ClusterSMBHAccretionTime         = 5.0;
  ClusterSMBHJetDim                = 2;
  ClusterSMBHAccretionEpsilon      = 0.001;
  ClusterSMBHDiskRadius            = 0.5;
  ClusterSMBHBCG                   = 1.0;
  ClusterSMBHMass                  = 0.0;

  EllipticalGalaxyRe               = 0.0;
  OldStarFeedbackAlpha             = 0.0;
  SNIaFeedbackEnergy               = 1.0;

  StellarWindSpeed                 = 0.0;   // in cgs
  StellarWindDensity               = 1.0; // in code unit
  StellarWindRadius                = 0.01;  // in code unit
  StellarWindTemperature           = 100.0;  // in K

  PythonTopGridSkip                = 0;
  PythonSubcycleSkip               = 1;
  PythonReloadScript               = FALSE;
  
  // EnzoTiming Dump Frequency
  TimingCycleSkip                  = 1;

  InlineHaloFinder                 = FALSE;
  HaloFinderSubfind                = FALSE;
  HaloFinderOutputParticleList     = FALSE;
  HaloFinderRunAfterOutput         = TRUE;
  HaloFinderMinimumSize            = 50;
  HaloFinderLinkingLength          = 0.2;
  HaloFinderCycleSkip              = 3;
  HaloFinderTimestep               = FLOAT_UNDEFINED;
  HaloFinderLastTime               = 0.0;

  StarClusterUseMetalField         = FALSE;
  StarClusterUnresolvedModel       = FALSE;
  StarClusterHeliumIonization      = FALSE;
  StarClusterMinDynamicalTime      = 10e6;         // in years
  StarClusterIonizingLuminosity    = 1e47;         // ph/s / Msun
  StarClusterSNEnergy              = 6.8e48;       // erg / Msun (Woosley&Weaver86)
  StarClusterSNRadius              = 10;           // pc
  StarClusterFormEfficiency        = 0.1;
  StarClusterMinimumMass           = 1000;         // Msun
  StarClusterCombineRadius         = 10;           // pc
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    StarClusterRegionLeftEdge[dim] = 0.0;
    StarClusterRegionRightEdge[dim] = 1.0;
  }
 
  MixSpeciesAndColors           = 1;            //Enable SNColour field to be advected as species in MHD
 
  PopIIIStarMass                   = 100;
  PopIIIInitialMassFunction        = FALSE;
  PopIIIInitialMassFunctionSeed    = INT_UNDEFINED;
  PopIIIInitialMassFunctionCalls   = 0;
  PopIIIHeliumIonization           = FALSE;
  PopIIILowerMassCutoff            = 1.0;
  PopIIIUpperMassCutoff            = 300.0;
  PopIIIInitialMassFunctionSlope   = -1.3;         // high mass slope
  PopIIIBlackHoles                 = FALSE;
  PopIIIBHLuminosityEfficiency     = 0.1;
  PopIIIOverDensityThreshold       = 1e6;          // times mean total density
  PopIIIH2CriticalFraction         = 5e-4;
  PopIIIMetalCriticalFraction      = 1e-4;
  PopIIISupernovaRadius            = 1;            // pc
  PopIIISupernovaUseColour         = FALSE;
  PopIIISupernovaMustRefine        = FALSE;
  PopIIISupernovaMustRefineResolution = 32;
  PopIIIColorDensityThreshold      = 1e6;          // times mean total density
  PopIIIColorMass                  = 1e6;          // total mass to color
  PopIIIUseHypernova               = TRUE;         // TRUE for HN yields, FALSE for CCSN
  PopIIISupernovaExplosions        = TRUE;         // TRUE for supernova energy injection
  PopIIIOutputOnFeedback           = FALSE;        // TRUE to output at creation and supernova
  IMFData                          = NULL;

  MBHAccretion                     = FALSE;        // 1: Bondi rate, 2: fix temperature, 3: fix rate, 4: Bondi with v_rel=0, 5: Bondi with v_rel=0 and vorticity
  MBHAccretionRadius               = 50;           // pc
  MBHAccretingMassRatio            = 1.0;          // 100%, check Star_CalculateMassAccretion.C
  MBHAccretionFixedTemperature     = 3e5;          // K,       for MBHAccretion = 2
  MBHAccretionFixedRate            = 1e-3;         // Msun/yr, for MBHAccretion = 3
  MBHTurnOffStarFormation          = FALSE;        // check Grid_StarParticleHandler.C
  MBHCombineRadius                 = 50;           // pc
  MBHMinDynamicalTime              = 10e6;         // in years
  MBHMinimumMass                   = 1e3;          // Msun

  MBHFeedback                      = FALSE;        // 1: isotropic thermal, 2: jet along z, 3: jet along L, 4: jet along L with 10deg noise, 5: jet along random direction
  MBHFeedbackRadiativeEfficiency   = 0.1;          // Shakura & Sunyaev (1973)
  MBHFeedbackEnergyCoupling        = 0.05;         // Springel (2005), Di Matteo (2005)
  MBHFeedbackMassEjectionFraction  = 0.1;          // 10%, check Star_CalculateFeedbackParameters.C
  MBHFeedbackMetalYield            = 0.02;         // 2%, check Star_CalculateFeedbackParameters.C
  MBHFeedbackThermalRadius         = 50;           // pc
  MBHFeedbackJetsThresholdMass     = 10;           // Msun

  /* Star Class MBH Paricle IO (PARTICLE_TYPE_MBH) */
  MBHParticleIO                    = FALSE;
  MBHParticleIOFilename            = (char*) "mbh_particle_io.dat";
  MBHInsertLocationFilename        = (char*) "mbh_insert_location.in";
  OutputWhenJetsHaveNotEjected     = FALSE;

  H2StarMakerEfficiency = 0.01;
  H2StarMakerNumberDensityThreshold = 0.0;
  H2StarMakerMinimumMass = 0.0;
  H2StarMakerMinimumH2FractionForStarFormation = 1e-5;
  H2StarMakerStochastic = 1;
  H2StarMakerUseSobolevColumn = 0;
  H2StarMakerSigmaOverR = 1.0/30.0;
  H2StarMakerAssumeColdWarmPressureBalance = 0;
  H2StarMakerH2DissociationFlux_MW = 1.0;
  H2StarMakerH2FloorInColdGas = 0.0;
  H2StarMakerColdGasTemperature = 1e4;

  // Paramters for Individual Star Star formation and Feedback

  /* PopIII Stars */
  TypeIILowerMass = 11.0;
  TypeIIUpperMass = 40.0;
  PISNLowerMass   = 140.0;
  PISNUpperMass   = 260.0;
  PopIIIPISNEnergy = -1.0;

  UseFUVBackground                     =    1;     // apply FUV background from UVB for Pe heating - only HM2012
  MetalMixingExperiment                =    0;     // turn on to add by-hand injection events to examine enrichment (see )
  IndividualStarRefineToLevel          =   -1;     // Level to force-refine to aroudn star particles (region buffer size in radius)
  IndividualStarRefineForRadiation     =    0;     // By default, above refine to level ONLY for mass / energy feedback. Also do when stars are ionizing if ON
  IndividualStarRefineBufferSize       =    4;
  /* IndividualStar: Star Formation */
  IndividualStarTemperatureLimit       =   -1;
  IndividualStarTemperatureLimitFactor =    2;
  IndividualStarICSupernovaRate        =   0.0;     // rate (# / yr) of supernova for initial driving - off if zero
  IndividualStarICSupernovaTime        =  10.0;     // (Myr) length of time to do SN driving if true - rate dec linearly starting at 1/2 this time
  IndividualStarICSupernovaR           =   -1.0;     // radius limit for SN - set to galaxy scale radius
  IndividualStarICSupernovaZ           =   -1.0;     // vertical height limit for supernova - set to Galaxy scale height
  IndividualStarICSupernovaMethod      =     1;     // 1: uniform disk with R and Z - 2: uniform sphere using R only
  IndividualStarICSupernovaFromFile    =     0;     // if ON, loads SNR as a function of time from a file.
  IndividualStarICSupernovaInjectionMethod = 1;   // 1 = thermal only, spherical ball - 2: thermal and kinetic mix - Simpson et. al. 2015
  for (i = 0; i < 3; i ++)
      IndividualStarICSupernovaPos[i] = 0.5;

  IndividualStarWDFixedLifetime      =  -1.0;         // debugging parameter - fixed lifetime in Myr before SNIa if > 0
  IndividualStarCreationStencilSize  =     3;         // n x n cell region (on each side) to sample for star formation
  IndividualStarCheckVelocityDiv     =     1;         // use velocity divergence in SF check
  IndividualStarICLifetimeMode       =     0;         // 0 - use interpolated lifetime, 1 - set to now, 2 - from file
  IndividualStarRefineTime           =    0.1;        // for particles with any SN, treat as must refine for this long after death (in Myr)
  // StarParticleOverdensityThreshold is used as primary density threshold parameter
  IndividualStarSecondaryOverDensityThreshold = -1;  // in cc - if < 0, set to over density thresh in ReadParamFile
  IndividualStarTemperatureThreshold = 1.0E4;       // threshold for star formation (T < T_thresh)
  IndividualStarMassFraction         =   0.5;         // Maximum fraction of region that can be converted into stars in one timestep
  IndividualStarSFGasMassThreshold   = 200.0;         // for SF algorithm 1, size of mass chunk that will be 100% converted to stars
  IndividualStarPopIIIFormation      =     0;         // flag to allow for Pop III star formation in gas below a Z threshold
                                                      // also tracks PopIII metal enrichment w/ additional field (but ONLY when tags are written to file)

  IndividualStarRProcessModel        =     0;         // Include an R-process model (TBD) and independent tracer field
  IndividualStarTrackAGBMetalDensity =     0;         // Track separate AGB metal mass field (only for outut chemical tags)
  IndividualStarTrackSNMetalDensity  =     0;         // Track spearate SNII AND SNIa mass fields (only for output chemical tags)

      // for popIII SF for individual stars, use the same PopIII IMF flags that already exist
      // 
  IndividualStarIMF                  =     0;         // 0: salpeter, 1: kroupa, 2: chabrier
  IndividualStarIgnoreNegativeMass   = 0; // Experimental!!!
  IndividualStarIMFCalls             =     0;         // Do not touch - number of calls to IMF so far in simulation
  IndividualStarSalpeterSlope        = -1.35;         // slope
  IndividualStarKroupaAlpha1         =  -0.3;         // kroupa slope over mass range
  IndividualStarKroupaAlpha2         =  -1.3;         // '' over 
  IndividualStarKroupaAlpha3         =  -2.3;         // '' over
  IndividualStarIMFLowerMassCutoff   =   1.0;         // Solar masses
  IndividualStarIMFUpperMassCutoff   = 100.0;         // Solar masses
  IndividualStarIMFMassFloor         = IndividualStarIMFLowerMassCutoff;         // If this is above the lower mass cutoff, stars below this mass will get grouped together into a single particle
  IndividualStarVelocityDispersion   =   1.0;         // initial velocity disperion of stars in SF region (km/s)
  IndividualStarIMFSeed              = INT_UNDEFINED; // random number seed for IMF sampling

  /* IndividualStar: Stellar Feedback (non-radiation) */
  IndividualStarFeedbackOverlapSample = 32;         // number of points per cell to compute fractional overlap in feedback routine
  IndividualStarFeedbackStencilSize   = 3;          // Size of feedback injection region (radius in number of cells)
  IndividualStarFeedbackRadius        = -1.0;       // Used (over above) if > 0. Feedback radius in pc (not # of grid zones)

  IndividualStarStellarWinds          = 1;          // on or off
  IndividualStarWindTemperature       = 1.0E6;      // temperature cap on stellar wind source region (K)
  IndividualStarUseWindMixingModel    = 0;          // account for unresolved mixing at wind/ISM shell interface and allow mass loading
  IndividualStarStellarWindVelocity   = -1;         // when < 0, use Leithener et. al. model for stellar wind velocities
                                                    // when > 0, uniform wind velocity for all stars in km / s
                                                    // when = 0, use this to do mass deposition without energy injection
  IndividualStarMaximumStellarWindVelocity = 3000.0; // km/s - maximum wind velocity in all cases
  IndividualStarAGBWindVelocity       = 20.0;       // km/s - wind velocity for all AGB stars

  IndividualStarAGBThreshold            = 8.0;      // solar masses - stars below this mass have winds at end of life only
  IndividualStarSNIIMassCutoff          = 8.0;      // Solar masses - stars above this go core collapse supernova
  IndividualStarDirectCollapseThreshold = 25.0;     // solar masses - no SNII and no ejecta above this mass
  IndividualStarSupernovaEnergy         = 1;        // when < 0, use factor x mc^2 for supernova energy injection
                                                    // when > 0, constant supernova energy in units of 10^51 erg
  IndividualStarPrintSNStats            = 0;        // print out grid and density info for each SN explosion

  IndividualStarUseSNIa                 = 0;        // on or off
  IndividualStarDTDSlope                = 1.20;     // beta (positive) - Default from Maoz et. al. 2012
  IndividualStarWDMinimumMass           = 1.7;      // Solar masses - min MS projenitor mass that forms WD
  IndividualStarWDMaximumMass           = 8.0;      // solar masses - max MS projenitor mass that forms WD
  IndividualStarSNIaMinimumMass         = 3.0;      // Solar masses - min MS projenitor mass that goes SNIa after WD formation
  IndividualStarSNIaMaximumMass         = 8.0;      // solar masses - max MS projenitor mass that goes SNIa after WD formation
  IndividualStarSNIaFraction            = 0.043;    // Fraction of MS stars that can be SNIa progenitors that will go SNIa in hubble time

  /* IndividualStar: Yields Tracking */
  IndividualStarFollowStellarYields  = 0;           // on or off
  IndividualStarSurfaceAbundances    = 0;           // return surface abundances in addition to yields
  IndividualStarExtrapolateYields    = 0;           // on or off - extrapolate yields from NuGrid set by scaling. If off, use PARSEC wind yields
  IndividualStarOutputChemicalTags   = 0;           // on or off - if ON does not tag particles, but instead outputs them to a file
  IndividualStarChemicalTagFilename  = NULL;        // filename for above
  IndividualStarSaveTablePositions   = 1;           // save table positions as particle attributes to save time

  /* IndividualStar: Stellar Feedback - Radiation */
  IndividualStarRadiationMinimumMass = 8.0;         // Solar masses - Stars above this are posible rad sources
  IndividualStarOTRadiationMass      = 8.0;         // Solar masses - Stars above this are allowed to have optically thin radiation
  IndividualStarIonizingRadiationMinimumMass = 8.0; // Solar masses - stars above this are allowed ionizing radiation
  IndividualStarFUVHeating           = 0;           // on or off - include Bakes & Tielens Pe Heating using optically thin FUV
  IndividualStarLWRadiation          = 0;           // on or off - include optically thin LW photons from stars
  IndividualStarFUVTemperatureCutoff = 2.0E4;       // K - if FUV heating is on, heat up to this temperature
  IndividualStarBlackBodyOnly        = 0;           // on or off - On = BB spectrum only - Off = OSTAR2002 when applicable

  IndividualStarBlackBodyq0Factors[0]  = 0.1 ;      // if OSTAR is ON, adjust black body to be continious
  IndividualStarBlackBodyq0Factors[1]  = 3.2 ;      // factors are for q0 and q1 with first
  IndividualStarBlackBodyq1Factors[0]  = 0.001 ;    // applying to low mass stars off of the grid
  IndividualStarBlackBodyq1Factors[1]  = 4.000 ;    // and second high mass off grid.

  IndividualStarBlackBodyFUVFactors[0] = 1.0;   // same as above, but for FUV luminosities
  IndividualStarBlackBodyFUVFactors[1] = 1.0;

  IndividualStarBlackBodyLWFactors[0]  = 1.0;     // same, but for LW
  IndividualStarBlackBodyLWFactors[1]  = 1.0;

  PhotoelectricHeatingDustModel = 1; // 0 - no shielding, linear in metallicity ; 1 - approx local shield, dust to gas ratio model
  PhotoelectricHeatingDustModelEfficiency = 0.0;    // Pe heating efficiency - <= 0 uses fit to Wolfire et. al. 2003 at solar radius

  /* Stellar Yields Parameters */
  StellarYieldsNumberOfSpecies       = INT_UNDEFINED; // number of species to follow - optional, calculated automatically if left undefined
  StellarYieldsScaledSolarInitialAbundances = 0;    // use solar abundances to set initial mass fractions, linearly scaled by metalliticy

  StellarAbundancesFilename = NULL;

  for (i = 0; i < MAX_STELLAR_YIELDS; i++){
    StellarYieldsAtomicNumbers[i] = -1;
  }
  for (i = 0; i < MAX_STELLAR_YIELDS; i++){
    StellarYieldsResetAtomicNumbers[i] = -1;
  }
  ResetStellarAbundances = 0;

  // chemical evolution test star
  ChemicalEvolutionTestNumberOfStars   = 1;
  ChemicalEvolutionTestStarMass        = 20.0   ; // solar masses
  ChemicalEvolutionTestStarMetallicity = 0.0001 ;
  ChemicalEvolutionTestStarFormed      = FALSE  ;
  ChemicalEvolutionTestStarLifetime    = 0      ; // 0 is off > 0 lifetime in Myr
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    ChemicalEvolutionTestStarPosition[dim] = 0.5; // code units - center of box
    ChemicalEvolutionTestStarVelocity[dim] = 0.0;
  }

  StarMakerMinimumMassRamp = 0;
  StarMakerMinimumMassRampStartTime = FLOAT_UNDEFINED;
  StarMakerMinimumMassRampStartMass = FLOAT_UNDEFINED;
  StarMakerMinimumMassRampEndTime   = FLOAT_UNDEFINED;
  StarMakerMinimumMassRampEndMass   = FLOAT_UNDEFINED;

  NumberOfParticleAttributes       = INT_UNDEFINED;
  AddParticleAttributes            = FALSE;
  LastSupernovaTime                = FLOAT_UNDEFINED;

  for (i = 0; i<MAX_MOVIE_FIELDS; i++)
    MovieDataField[i] = INT_UNDEFINED;
  MovieSkipTimestep = INT_UNDEFINED;
  NewMovieName = DefaultNewMovieName;
  NewMovieDumpNumber = 0;
  NewMovieParticleOn = FALSE;
  Movie3DVolumes  = FALSE;
  MovieVertexCentered = FALSE;
  MetaData.MovieTimestepCounter      = 0;

  ran1_init = 0;
  rand_init = 0;

  SinkMergeDistance     = 1e16;
  SinkMergeMass         = 0.1;
  TotalSinkMass         = 0.0;
  StellarWindFeedback   = 0;
  StellarWindTurnOnMass = 0.1;
  MSStellarWindTurnOnMass = 10.0;

  UseHydro		     = 1;
  Coordinate		     = Cartesian;
  NSpecies		     = INT_UNDEFINED;
  NColor		     = INT_UNDEFINED;
  Theta_Limiter		     = 1.5;
  UsePhysicalUnit	     = 0;
  NEQ_HYDRO		     = 5;
  NEQ_MHD		     = 9;
  SmallRho		     = 1e-30;
  SmallP		     = 1e-35;
  SmallEint		     = 1e-30;
  SmallT		     = 1e-10;
  MaximumAlvenSpeed	     = 1e30;
  RiemannSolver		     = INT_UNDEFINED;
  RiemannSolverFallback      = 1;
  ReconstructionMethod	     = INT_UNDEFINED;
  PositiveReconstruction     = FALSE;
  ConservativeReconstruction = 0;
  EOSType		     = 0;
  EOSSoundSpeed		     = 2.65e4;
  EOSCriticalDensity	     = 1e-13;
  EOSGamma		     = 1.667;
  Mu			     = 0.6;
  DivBDampingLength          = 1.;
  UseCUDA		     = 0;
  UseFloor		     = 0;
  UseViscosity		     = 0;
  ViscosityCoefficient       = 0.;
  UseAmbipolarDiffusion	     = 0;
  UseResistivity	     = 0;

  StringKick = 0;
  StringKickDimension = 0;

  iden	= 0;
  ivx	= 1;
  ivy	= 2;
  ivz	= 3;
  ietot = 4;
  ieint = 0;
  iBx	= 5;
  iBy	= 6;
  iBz	= 7;
  iPhi	= 8;
  iD	= 0;
  iS1	= 1;
  iS2	= 2;
  iS3	= 3;
  iEtot = 4;
  iEint = 0;

  UsePoissonDivergenceCleaning		   = 0;
  PoissonDivergenceCleaningBoundaryBuffer = 0;
  PoissonDivergenceCleaningThreshold	   = 0.001;
  PoissonApproximationThreshold	   = 0.001;
  PoissonBoundaryType	   = 0;

  UseDrivingField   = 0;
  DrivingEfficiency = 1.0;

  /* End of Stanford Hydro additions */

  /* test problem values */
  TestProblemData.HydrogenFractionByMass = 0.76;

  /* The DToHRatio is by mass in the code, so multiply by 2. */
  TestProblemData.DeuteriumToHydrogenRatio = 2.0*3.4e-5; // Burles & Tytler 1998 

  // multispecies default values assume completely neutral gas with primordial D/H ratio
  TestProblemData.MultiSpecies = 0;
  TestProblemData.HI_Fraction = 1.0;
  TestProblemData.HII_Fraction = tiny_number;
  TestProblemData.HeI_Fraction = 1.0;
  TestProblemData.HeII_Fraction = tiny_number;
  TestProblemData.HeIII_Fraction = tiny_number;
  TestProblemData.HM_Fraction = tiny_number;
  TestProblemData.H2I_Fraction = tiny_number;
  TestProblemData.H2II_Fraction = tiny_number;
  TestProblemData.DI_Fraction = 2.0*3.4e-5;
  TestProblemData.DII_Fraction = tiny_number;
  TestProblemData.HDI_Fraction = tiny_number;

  // This is for ionized gas (within a supernova blast, for example)
  TestProblemData.HI_Fraction_Inner = 1.0;
  TestProblemData.HII_Fraction_Inner = tiny_number;
  TestProblemData.HeI_Fraction_Inner = 1.0;
  TestProblemData.HeII_Fraction_Inner = tiny_number;
  TestProblemData.HeIII_Fraction_Inner = tiny_number;
  TestProblemData.HM_Fraction_Inner = tiny_number;
  TestProblemData.H2I_Fraction_Inner = tiny_number;
  TestProblemData.H2II_Fraction_Inner = tiny_number;
  TestProblemData.DI_Fraction_Inner = 2.0*3.4e-5;
  TestProblemData.DII_Fraction_Inner = tiny_number;
  TestProblemData.HDI_Fraction_Inner = tiny_number;

  TestProblemData.UseMetallicityField = 0;
  TestProblemData.MetallicityField_Fraction = tiny_number;
  TestProblemData.MetallicitySNIaField_Fraction = tiny_number;
  TestProblemData.MetallicitySNIIField_Fraction = tiny_number;

  TestProblemData.UseMassInjection = 0;
  TestProblemData.InitialHydrogenMass = tiny_number;
  TestProblemData.InitialDeuteriumMass = tiny_number;
  TestProblemData.InitialHeliumMass = tiny_number;
  TestProblemData.InitialMetalMass = tiny_number;

  TestProblemData.MultiMetals = 0;
  TestProblemData.MultiMetalsField1_Fraction = tiny_number;
  TestProblemData.MultiMetalsField2_Fraction = tiny_number;

  // AJE alpha process, iron, s and r process elements
  TestProblemData.NI_Fraction  = tiny_number;
  TestProblemData.MgI_Fraction = tiny_number;
  TestProblemData.FeI_Fraction = tiny_number;
  TestProblemData.YI_Fraction  = tiny_number;
  TestProblemData.BaI_Fraction = tiny_number;
  TestProblemData.LaI_Fraction = tiny_number;
  TestProblemData.EuI_Fraction = tiny_number;
  // same for second set of chemical tracers:
  TestProblemData.CI_Fraction_2 = tiny_number;
  TestProblemData.NI_Fraction_2 = tiny_number;
  TestProblemData.OI_Fraction_2 = tiny_number;
  TestProblemData.MgI_Fraction_2 = tiny_number;
  TestProblemData.SiI_Fraction_2 = tiny_number;
  TestProblemData.FeI_Fraction_2 = tiny_number;
  TestProblemData.YI_Fraction_2  = tiny_number;
  TestProblemData.LaI_Fraction_2 = tiny_number;
  TestProblemData.BaI_Fraction_2 = tiny_number;
  TestProblemData.EuI_Fraction_2 = tiny_number;

  for (int i = 0; i < MAX_STELLAR_YIELDS; i++){
    TestProblemData.ChemicalTracerSpecies_Fractions[i]   = tiny_number;
    TestProblemData.ChemicalTracerSpecies_Fractions_2[i] = tiny_number;
  }

  TestProblemData.GloverChemistryModel = 0;
  // This is for the gas in the surrounding medium, for the blast wave problem.
  TestProblemData.CI_Fraction = tiny_number;
  TestProblemData.CII_Fraction = tiny_number;
  TestProblemData.OI_Fraction = tiny_number;
  TestProblemData.OII_Fraction = tiny_number;
  TestProblemData.SiI_Fraction = tiny_number;
  TestProblemData.SiII_Fraction = tiny_number;
  TestProblemData.SiIII_Fraction = tiny_number;
  TestProblemData.CHI_Fraction = tiny_number;
  TestProblemData.CH2I_Fraction = tiny_number;
  TestProblemData.CH3II_Fraction = tiny_number;
  TestProblemData.C2I_Fraction = tiny_number;
  TestProblemData.COI_Fraction = tiny_number;
  TestProblemData.HCOII_Fraction = tiny_number;
  TestProblemData.OHI_Fraction = tiny_number;
  TestProblemData.H2OI_Fraction = tiny_number;
  TestProblemData.O2I_Fraction = tiny_number;

  // This is for the gas in the region where the blast wave energy is injected
  TestProblemData.CI_Fraction_Inner = tiny_number;
  TestProblemData.CII_Fraction_Inner = tiny_number;
  TestProblemData.OI_Fraction_Inner = tiny_number;
  TestProblemData.OII_Fraction_Inner = tiny_number;
  TestProblemData.SiI_Fraction_Inner = tiny_number;
  TestProblemData.SiII_Fraction_Inner = tiny_number;
  TestProblemData.SiIII_Fraction_Inner = tiny_number;
  TestProblemData.CHI_Fraction_Inner = tiny_number;
  TestProblemData.CH2I_Fraction_Inner = tiny_number;
  TestProblemData.CH3II_Fraction_Inner = tiny_number;
  TestProblemData.C2I_Fraction_Inner = tiny_number;
  TestProblemData.COI_Fraction_Inner = tiny_number;
  TestProblemData.HCOII_Fraction_Inner = tiny_number;
  TestProblemData.OHI_Fraction_Inner = tiny_number;
  TestProblemData.H2OI_Fraction_Inner = tiny_number;
  TestProblemData.O2I_Fraction_Inner = tiny_number;

  TestProblemData.MinimumHNumberDensity = 1;
  TestProblemData.MaximumHNumberDensity = 1e6;
  TestProblemData.MinimumMetallicity    = 1e-6;
  TestProblemData.MaximumMetallicity    = 1;
  TestProblemData.MinimumTemperature    = 10;
  TestProblemData.MaximumTemperature    = 1e7;
  TestProblemData.ResetEnergies         = 1;

  // This should only be false for analysis.
  // It could also be used (cautiously) for other purposes.
  LoadGridDataAtStart = TRUE;

  IsothermalSoundSpeed = 1.0;

  MetalCooling = FALSE;
  MetalCoolingTable = (char*) "metal_cool.dat";

#ifdef USE_PYTHON
  NumberOfPythonCalls = 0;
  NumberOfPythonTopGridCalls = 0;
  NumberOfPythonSubcycleCalls = 0;
  grid_dictionary = PyDict_New();
  old_grid_dictionary = PyDict_New();
  hierarchy_information = PyDict_New();
  yt_parameter_file = PyDict_New();
  conversion_factors = PyDict_New();
  my_processor = PyLong_FromLong((Eint) MyProcessorNumber);
#endif

  /* Some stateful variables for EvolveLevel */
  for(i = 0; i < MAX_DEPTH_OF_HIERARCHY; i++) {
    LevelCycleCount[i] = LevelSubCycleCount[i] = 0;
    dtRebuildHierarchy[i] = -1.0;
    TimeSinceRebuildHierarchy[i] = 0.0;
    dtThisLevelSoFar[i] = dtThisLevel[i] = 0.0;
  }

  /* Shearing Boundary Conditions variables */

  AngularVelocity=0.001;
  VelocityGradient=1.0;
  ShearingBoundaryDirection=-1;
  ShearingVelocityDirection=-1;
  ShearingBoxProblemType = 0; 
  UseMHD=0;
  MaxVelocityIndex = 3;

  //MHDCT variables
  MHDCTSlopeLimiter = 1;
  MHDCTDualEnergyMethod = INT_UNDEFINED;
  MHDCTPowellSource = 0;
  MHDCTUseSpecificEnergy = TRUE;
  WriteBoundary             = FALSE;
  CT_AthenaDissipation = 0.1;
  MHD_WriteElectric = TRUE;
  tiny_pressure = tiny_number;
  MHD_CT_Method = 2;
  NumberOfGhostZones = 3;
  IsothermalSoundSpeed = 1.0;
  MHD_ProjectB = FALSE;
  MHD_ProjectE = TRUE;
  UseMHDCT = FALSE;
  EquationOfState = 0;
  for(int dccdbg=0; dccdbg<MAX_EXTRA_OUTPUTS;dccdbg++) ExtraOutputs[dccdbg]=INT_UNDEFINED;
  WriteAcceleration = FALSE;

  CorrectParentBoundaryFlux = FALSE; //Corrects (or doesn't) parent flux when subgrid shares an exposed face with parent.

  for(int dccdbg=0; dccdbg<MAX_EXTRA_OUTPUTS;dccdbg++) ExtraOutputs[dccdbg]=INT_UNDEFINED;
  MoveParticlesBetweenSiblings = TRUE;

  /* Particle Splitter */

  ParticleSplitterIterations = FALSE;
  ParticleSplitterChildrenParticleSeparation = 1.0;
  ParticleSplitterRandomSeed = 131180;
  ParticleSplitterMustRefine = FALSE;
  ParticleSplitterMustRefineIDFile = NULL;
  for(int i = 0; i < MAX_SPLIT_ITERATIONS; i++)
    ParticleSplitterFraction[i] = 1.0;
  for(int i = 0; i < MAX_DIMENSION; i++)
    ParticleSplitterCenter[i] = -1.0;
  for(int i = 0; i < MAX_SPLIT_ITERATIONS; i++)
    ParticleSplitterCenterRegion[i] = -1.0;

  /* Magnetic Field Resetter */

  ResetMagneticField = FALSE;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    ResetMagneticFieldAmplitude[dim] = 0.0;   // in Gauss
  }  

  VelAnyl                     = 0;
  BAnyl                       = 0;
  WriteExternalAccel          = 0;

  for (i = 0; i < MAX_ACTIVE_PARTICLE_TYPES; i++)
    EnabledActiveParticles[i] = NULL;
  EnabledActiveParticlesCount = 0;
  UnfulfilledStarFormationMass = 0;
  NextActiveParticleID = INT_UNDEFINED;
  NumberOfActiveParticles = 0;
  ActiveParticleDensityThreshold = 1e8; //in cm^-3
  //SmartStar Feedback modes
  SmartStarFeedback = FALSE;
  SmartStarEddingtonCap = FALSE;
  SmartStarBHFeedback = FALSE;
  SmartStarBHRadiativeFeedback = FALSE;
  SmartStarBHJetFeedback = FALSE;
  SmartStarBHThermalFeedback = FALSE;
  SmartStarStellarRadiativeFeedback = FALSE;
  
  //SmartStar Feedback parameters - should be as minimal as possible
  SmartStarFeedbackEnergyCoupling = 0.016666;
  SmartStarFeedbackJetsThresholdMass = 1.0;
  SmartStarJetVelocity = 1e-1; //as a fraction of clight
  SmartStarSuperEddingtonAdjustment = TRUE;
  SmartStarSpin = 0.7;
  SmartStarSMSLifetime = 1e6; //1 Myr
  /* Gas drag parameters */
  UseGasDrag = 0;
  GasDragCoefficient = 0.;

  /* Supernova magnetic feedback */
  /* Default == 0 -> no magnetic field contribution */
  UseMagneticSupernovaFeedback = 0;   // 0 = off; 1+ = on;  2 = Radius and duration calculated during runtime
  MagneticSupernovaRadius = 300.0;      // Total injection radius of magnetic field in parsecs
  MagneticSupernovaDuration = 5e4  ;    // Total duration of magnetic feedback in years
  MagneticSupernovaEnergy = 1.0e51;    // Total energy (ergs) injected per star particle (supernova)

  return SUCCESS;
}
