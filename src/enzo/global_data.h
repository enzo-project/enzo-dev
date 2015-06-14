/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       February, 2004
/  modified:   Robert Harkness, August 12th 2006
/
/  PURPOSE:
/    This is the global data, which should be held to a minimum.  Any changes
/    in this file require changes in: WriteGlobalData,
/    ReadGlobalData and InitializeNew.  
/    This file is dual-purposed:
/        1) read with    DEFINE_STORAGE defined for the (single) definition
/        2) read without DEFINE_STORAGE defined for external linkage
/
************************************************************************/
#ifndef GLOBAL_DATA_DEFINED__
#define GLOBAL_DATA_DEFINED__

#include <stdio.h>
#ifdef MEMORY_POOL
#include "MemoryPool.h"
#endif
#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

#ifdef NEW_PROBLEM_TYPES
class EnzoProblemType;
#endif

/* Load Balancing.  Currently only memory count method implemented
                          0 = off
                          1 = Equalize processor memory count
                         2 = Load balance only on a node
*/
EXTERN int NumberOfGhostZones;
EXTERN int LoadBalancing;
EXTERN int LoadBalancingCycleSkip;
EXTERN int ResetLoadBalancing;
EXTERN int CoresPerNode;
EXTERN int PreviousMaxTask;
EXTERN int LoadBalancingMinLevel;
EXTERN int LoadBalancingMaxLevel;

/* FileDirectedOutput checks for file existence: 
   stopNow (writes, stops),   outputNow, subgridcycleCount */
EXTERN int FileDirectedOutput;

/* These two flags determine the format of the hierarchy file for
   input and output:

   HierarchyFileInputFormat  = 0 -- HDF5 (default)
   HierarchyFileInputFormat  = 1 -- ASCII

   HierarchyFileOutputFormat = 0 -- HDF5 (default)
   HierarchyFileOutputFormat = 1 -- ASCII
   HierarchyFileOutputFormat = 2 -- both HDF5 and ASCII
*/
EXTERN int HierarchyFileInputFormat;
EXTERN int HierarchyFileOutputFormat;

/* LevelLookupTable is read in from the HDF5 hierarchy file. Its
   purpose is to allow one to quickly determine while level an input
   grid is on from its grid ID. This is needed to access the
   corresponding hierarchy dataset. The array is only valid during the
   reading of the HDF5 hierarchy file (i.e. in ReadAllData and
   below). */
EXTERN int *LevelLookupTable;
/* This is useful for loops over all grids. */
EXTERN int TotalNumberOfGrids;

/* debugging, extraction flags */

EXTERN int debug;
EXTERN int debug1;
EXTERN int debug2;
EXTERN int extract;

/* Problem: 00 = None                    01 = ShockTube
            02 = WavePool                03 = ShockPool  
	    04 = Double-Mach reflection  05 = ShockInABox
	    06 = Implosion               07 = SedovBlast
	    08 = KelvinHelmholtz instability
	    20 = 1D Zeldovich Pancake    21 = 1D pressureless collapse
	    22 = Adiabatic expansion     23 = TestGravity
            24 = Spherical infall        25 = TestGravitySphere
	    26 = GravityEquilibriumTest  27 = CollapseTest
	    28 = TestGravityMotion
	    30 = Cosmology simulation
	    50 = ThermalInstability simulation
	    51 = ThermalPancake test
	    60 = TurbulenceSimulation
	                                                                  */
EXTERN int CheckpointRestart;
EXTERN int WriteGhostZones;
EXTERN int ReadGhostZones;
EXTERN int ProblemType;
#ifdef NEW_PROBLEM_TYPES
EXTERN char *ProblemTypeName;
EXTERN EnzoProblemType *CurrentProblemType;
#endif

/* Hydrodynamics method:
       0 - PPM_DE      1 - PPM_LR (not working)    2 - ZEUS    3 - RK hydro   4 - RK MHD    */

EXTERN hydro_method HydroMethod;

/* Large and small numbers (i.e. compared to any real quantity).  This may
   be machine and problem dependent. */

EXTERN float huge_number, tiny_number;

/* Gamma: Ideal gas law constant. */

EXTERN float Gamma;

/* Flag indicating if the gas is pressureless. */

EXTERN int PressureFree;

/* Factor to refine by */

EXTERN int RefineBy;

/* Maximum refinement level (0 = topgrid). */

EXTERN int MaximumRefinementLevel;
EXTERN int MaximumGravityRefinementLevel;
EXTERN int MaximumParticleRefinementLevel;
EXTERN int FastSiblingLocatorEntireDomain;

/* Cell Flagging method:  0 = None
                          1 = FlagCellsToBeRefinedBySlope
			  2 = FlagCellsToBeRefinedByMass (baryon only)
			  3 = FlagCellsToBeRefinedByShocks
			  4 = FlagCellsToBeRefinedByMass (particles only)
	     (disabled)	  5 = FlagCellsToBeRefinedByOverdensity (baryon only)
			  6 = FlagCellsToBeRefinedByJeansLength
                          7 = FlagCellsToBeRefinedByCoolingTime
                          8 = FlagCellsToBeRefinedByMustRefineParticles
                          9 = FlagCellsToBeRefinedByShear
			 11 = FlagCellsToBeRefinedByResistiveLength
                         12 = FlagCellsToBeRefinedByMustRefineRegion
			 13 = FlagCellsToBeRefinedByMetallicity
       15 = FlagCellsToBeRefinedBySecondDerivative
 */

EXTERN int CellFlaggingMethod[MAX_FLAGGING_METHODS];

/* left and right boundaries of the 'must refine region'
   for CellFlaggingMethod = 10 */

EXTERN FLOAT MustRefineRegionLeftEdge[MAX_DIMENSION];  // left edge
EXTERN FLOAT MustRefineRegionRightEdge[MAX_DIMENSION];  // right edge

/* left and right boundaries of the 'avoid refine region'
   for CellFlaggingMethod = 101 */

EXTERN int   AvoidRefineRegionLevel[MAX_STATIC_REGIONS];
EXTERN FLOAT AvoidRefineRegionLeftEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];
EXTERN FLOAT AvoidRefineRegionRightEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];

/* specifies the level to which FlagCellsToBeRefinedByMustRefineRegion
   will refine up to (does not prevent refinement to higher levels) */

EXTERN int MustRefineRegionMinRefinementLevel;

/* specifies the level to which FlagGridCellsToBeRefinedByMetallicity
   will refine up to (does not prevent refinement to higher levels) */
EXTERN int MetallicityRefinementMinLevel;

/* threshold metallicity and density for FlagGridCellsToBeRefinedByMetallicity */
EXTERN float MetallicityRefinementMinMetallicity;
EXTERN float MetallicityRefinementMinDensity;

/* Velocity to limit timesteps */

EXTERN float TimestepSafetyVelocity;

/* Flag indicating if the flux correction should be applied. */

EXTERN int FluxCorrection;

/* This specifies the interpolation method (see typedefs.h). */

EXTERN interpolation_type InterpolationMethod;
EXTERN int ConservativeInterpolation;

/* This is the minimum efficiency of combined grid needs to achieve in
   order to be considered better than the two grids from which it formed. */

EXTERN float MinimumEfficiency;

/* This flag will automatically adjust MinimumSubgridEdge and
   MaximumSubgridSize.  It will select MaximumSubgridSize from
   OptimalSubgridPerProcessor. */

EXTERN int SubgridSizeAutoAdjust;
EXTERN int OptimalSubgridsPerProcessor;

/* This is the minimum allowable edge size for a new subgrid (>=4) */

EXTERN int MinimumSubgridEdge;

/* This is the maximum allowable size for a new subgrid (>=2000) */

EXTERN int MaximumSubgridSize;

/* The number of zones that will be refined around each flagged zone. */

EXTERN int NumberOfBufferZones;

/* The left and right boundaries of the entire computational domain. */

EXTERN FLOAT DomainLeftEdge[MAX_DIMENSION], DomainRightEdge[MAX_DIMENSION];

/* Velocity of entire computational domain. */

EXTERN float GridVelocity[MAX_DIMENSION];

/* HDF names for labels and scales. */

EXTERN char *DimUnits[MAX_DIMENSION], *DimLabels[MAX_DIMENSION];
EXTERN char *DataLabel[MAX_NUMBER_OF_BARYON_FIELDS];
EXTERN char *DataUnits[MAX_NUMBER_OF_BARYON_FIELDS];

/* Region in which refinement is allowed (in problem space). */

EXTERN FLOAT RefineRegionLeftEdge[MAX_DIMENSION], 
             RefineRegionRightEdge[MAX_DIMENSION];
EXTERN int RefineRegionAutoAdjust;

EXTERN int MultiRefineRegion;
EXTERN FLOAT MultiRefineRegionLeftEdge[MAX_STATIC_REGIONS][MAX_DIMENSION], 
             MultiRefineRegionRightEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];
EXTERN int MultiRefineRegionGeometry[MAX_STATIC_REGIONS];
EXTERN FLOAT MultiRefineRegionCenter[MAX_STATIC_REGIONS][MAX_DIMENSION];
EXTERN FLOAT MultiRefineRegionOrientation[MAX_STATIC_REGIONS][MAX_DIMENSION];
EXTERN FLOAT MultiRefineRegionRadius[MAX_STATIC_REGIONS];
EXTERN FLOAT MultiRefineRegionWidth[MAX_STATIC_REGIONS];
EXTERN int MultiRefineRegionMaximumLevel[MAX_STATIC_REGIONS];
EXTERN int MultiRefineRegionMinimumLevel[MAX_STATIC_REGIONS];
EXTERN int MultiRefineRegionMaximumOuterLevel;
EXTERN int MultiRefineRegionMinimumOuterLevel;
EXTERN FLOAT MultiRefineRegionStaggeredRefinement[MAX_STATIC_REGIONS];

/* Uniform gravity: on/off flag, direction, and strength. */

EXTERN int UniformGravity, UniformGravityDirection;
EXTERN float UniformGravityConstant;

/* point source gravity: on/off flag position, and strength. */

EXTERN int PointSourceGravity;
EXTERN FLOAT PointSourceGravityPosition[MAX_DIMENSION];
EXTERN float PointSourceGravityConstant;
EXTERN float PointSourceGravityCoreRadius;

/* disk gravity */
EXTERN int DiskGravity;
EXTERN FLOAT DiskGravityPosition[MAX_DIMENSION],
             DiskGravityAngularMomentum[MAX_DIMENSION];
EXTERN float DiskGravityStellarDiskMass;
EXTERN float DiskGravityStellarDiskScaleHeightR;
EXTERN float DiskGravityStellarDiskScaleHeightz;
EXTERN float DiskGravityStellarBulgeMass;
EXTERN float DiskGravityStellarBulgeR;
EXTERN float DiskGravityDarkMatterR;
EXTERN float DiskGravityDarkMatterDensity;

/* SelfGravity (TRUE or FALSE) */

EXTERN int SelfGravity;
EXTERN int SelfGravityGasOff;
EXTERN int AccretionKernal;

/* CopyGravPotential (TRUE or FALSE) */

EXTERN int CopyGravPotential;

/* Number of iterations to solve the potential */

EXTERN int PotentialIterations;

/* Flag indicating whether or not to use the baryon self-gravity approximation
   (subgrid cells influence are approximated by their projection to the
   current grid). */

EXTERN int BaryonSelfGravityApproximation;

/* Coefficient in front of source term in Poisson's equations.
   (i.e. Del^phi = GravitationConstant * density, usually 4*Pi*G). */

EXTERN float GravitationalConstant;

/* S2 Particle size in top grid cell units (usually around 3).  The S2
   particle is S(r) = A*(a/2-r) (if r < a/2, 0 otherwise).  The constant
   A depends on the dimension: 1D) 4/a^2,  2D) 24/(Pi*a^3)  3D) 48/(Pi*a^3). */

EXTERN float S2ParticleSize;

/* Gravity resolution factor is a float indicating the comparative resolution
   of the gravitational computation compared to the grid (1-2 or so). */

EXTERN float GravityResolution;

/* Flag to indicate if gravitational potential field should be computed
   and stored. */

EXTERN int ComputePotential;

/* Flag to indicate output for gravitational potential field. */

EXTERN int WritePotential;

/* Parameter to control how particles in a subgrid are deposited in
   the target grid.  Options are: 
     CIC_DEPOSIT - cloud in cell using cloud size equal to target grid size
     CIC_DEPOSIT_SMALL - CIC using cloud size equal to source grid size
     NGP_DEPOSIT - nearest grid point */

EXTERN int ParticleSubgridDepositMode;

/* Maximum number of GreensFunctions that will be stored in any time.
   This number must be less than MAX_NUMBER_OF_GREENS_FUNCTIONS. */

EXTERN int GreensFunctionMaxNumber;

/* Maximum number of words associated with GreensFunction storage
   (Not currently implemented). */

EXTERN int GreensFunctionMaxSize;

/* Dual energy formalism (TRUE or FALSE). */

EXTERN int DualEnergyFormalism;

/* Two parameters for the dual energy formalism. */

EXTERN float DualEnergyFormalismEta1;
EXTERN float DualEnergyFormalismEta2;

/* This is the particle equivalent of the Courant factor.  It is the maximum
   number of cells a particle is allowed to travel in a single timestep. */

EXTERN float ParticleCourantSafetyNumber;

/* This is a parameter to control root grid time steps, and is basically
   a hack to ensure that star particles don't get ejected out of grids. */

EXTERN float RootGridCourantSafetyNumber;

/* Radiative cooling on/off flag and associated data. */

EXTERN int RadiativeCooling;
EXTERN CoolDataType CoolData;
EXTERN int RadiativeCoolingModel;

/* Cloudy cooling parameters and data. */

EXTERN CloudyCoolingDataType CloudyCoolingData;

/* Gadget Equilibrium cooling on/off flag */

EXTERN int GadgetEquilibriumCooling;

/* Random Forcing on/off flag and associated data. */ //AK

EXTERN int     RandomForcing;
EXTERN FLOAT   RandomForcingEdot;
EXTERN FLOAT   RandomForcingMachNumber;  //#####
EXTERN fpos_t  BaryonFileNamePosition;

/* Multi-species rate equation flag and associated data. */

EXTERN int MultiSpecies;
EXTERN int NoMultiSpeciesButColors;
EXTERN int ThreeBodyRate;
EXTERN RateDataType RateData;
EXTERN int H2FormationOnDust;

/* Glover chemistry/cooling network flags */
EXTERN int GloverChemistryModel;  // 0 is off, on is 1-7, excluding 6
EXTERN int GloverRadiationBackground; // 1: low Z, 2: ISM
EXTERN int GloverOpticalDepth; // 0: opticaly thin, 1: single-cell

/* Multi-element metallicity field flag and count. */

EXTERN int MultiMetals;

/* Cosmic Ray Model
 * 0: Off - default
 * 1: On, (two fluid model)
 */
EXTERN int CRModel;
/* Cosmic Ray Diffusion
 * 0: Off - default
 * 1: On, CRkappa is constant across grid
 */
EXTERN int CRDiffusion;
/* Cosmic Ray Feedback
 *    0.0 -- No CR feedback
 *    1.0 -- All feedback into CR field
 */
EXTERN float CRFeedback;
EXTERN float CRkappa;
EXTERN float CRCourantSafetyNumber;
EXTERN float CRdensFloor;
EXTERN float CRmaxSoundSpeed;
EXTERN float CRgamma;
EXTERN float CosmologySimulationUniformCR; // FIXME

/* Shock Finding Method
 * 0: Off - default
 * 1: temperature unsplit 
 * 2: temperature split 
 * 3: velocity unsplit
 * 4: velocity split
 */
EXTERN int ShockMethod; 
EXTERN float ShockTemperatureFloor;
EXTERN int StorePreShockFields;
EXTERN int FindShocksOnlyOnOutput;


/* Type of radiation field. 
   0 - none,                    1 - Haardt & Madau alpha=-1.5
   2 - H&M alpha = -1.8       
   10 - homogenous internal radiation field (a la Renyue's work) */

EXTERN int RadiationFieldType;
EXTERN int AdjustUVBackground; 
EXTERN int AdjustUVBackgroundHighRedshift; 
EXTERN float SetUVBAmplitude;
EXTERN float SetHeIIHeatingScale;
EXTERN RadiationFieldDataType RadiationData;
EXTERN int RadiationFieldLevelRecompute;
EXTERN int RadiationXRaySecondaryIon;
EXTERN int RadiationXRayComptonHeating;
EXTERN int TabulatedLWBackground;
EXTERN float RadiationFieldRedshift;

/* Photoelectric cooling turn on/off */

EXTERN int PhotoelectricHeating;
EXTERN float PhotoelectricHeatingRate;

/* Output cooling time with grid data. */

EXTERN int OutputCoolingTime;

/* Output temperature with grid data. */

EXTERN int OutputTemperature;

/* Output dust temperature with grid data. */

EXTERN int OutputDustTemperature;

/* Output smoothed dark matter fields. */

EXTERN int OutputSmoothedDarkMatter;
EXTERN int SmoothedDarkMatterNeighbors;

/* Output gridded star particle fields. */

EXTERN int OutputGriddedStarParticle;

/* ZEUS Hydro artificial viscosity parameters (C1, C2 of Stone & Norman). */

EXTERN float ZEUSLinearArtificialViscosity;
EXTERN float ZEUSQuadraticArtificialViscosity;

/* Parameters for MinimumPressureSupport. */

EXTERN int UseMinimumPressureSupport;
EXTERN float MinimumPressureSupportParameter;

/* Parameters for statically refined regions. */
EXTERN FLOAT StaticRefineRegionLeftEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];
EXTERN FLOAT StaticRefineRegionRightEdge[MAX_STATIC_REGIONS][MAX_DIMENSION];
EXTERN int   StaticRefineRegionLevel[MAX_STATIC_REGIONS];

/* Evolving refinement region. */
EXTERN char *RefineRegionFile;
EXTERN int RefineRegionTimeType; // 0=time 1=redshift
EXTERN int EvolveRefineRegionNtimes;
EXTERN FLOAT EvolveRefineRegionTime[MAX_REFINE_REGIONS]; // time bins
EXTERN FLOAT EvolveRefineRegionLeftEdge[MAX_REFINE_REGIONS][3]; // left corners
EXTERN FLOAT EvolveRefineRegionRightEdge[MAX_REFINE_REGIONS][3]; // right corners

/* Processor identifier for this thread/processor */

EXTERN int MyProcessorNumber;
EXTERN int NumberOfProcessors;
EXTERN float CommunicationTime;

/* Parameter to indicate if top grid should do parallel IO
   (currently only works for ProblemType == 30). */

EXTERN int ParallelRootGridIO;
EXTERN int ParallelParticleIO;
EXTERN int Unigrid;
EXTERN int CubeDumpEnabled;
EXTERN int PartitionNestedGrids;
EXTERN int StaticPartitionNestedGrids;
EXTERN int ExtractFieldsOnly;
EXTERN int First_Pass;
EXTERN int UnigridTranspose;
EXTERN int NumberOfRootGridTilesPerDimensionPerProcessor;
EXTERN int CosmologySimulationNumberOfInitialGrids;
EXTERN int UserDefinedRootGridLayout[3];

/* Parameters that control density dex output */

EXTERN int OutputOnDensity;
EXTERN float StartDensityOutputs;
EXTERN float CurrentDensityOutput;
EXTERN float CurrentMaximumDensity;
EXTERN float IncrementDensityOutput;

/* Parameter(s) for embedded python execution */
EXTERN int PythonTopGridSkip;
EXTERN int PythonSubcycleSkip;
EXTERN int PythonReloadScript;

/* Parameters to control inline halo finding */

EXTERN int InlineHaloFinder;
EXTERN int HaloFinderSubfind;
EXTERN int HaloFinderOutputParticleList;
EXTERN int HaloFinderMinimumSize;
EXTERN int HaloFinderCycleSkip;
EXTERN int HaloFinderRunAfterOutput;
EXTERN float HaloFinderLinkingLength;
EXTERN float HaloFinderTimestep;
EXTERN FLOAT HaloFinderLastTime;

/************************************************/
/* Global data for specific problems or methods */
/************************************************/

/* For CellFlaggingMethod = 1,
   The minimum relative slope (da/dx over a) required for refinement. */

EXTERN float MinimumSlopeForRefinement[MAX_FLAGGING_METHODS];
EXTERN int SlopeFlaggingFields[MAX_FLAGGING_METHODS];

/* For CellFlaggingMethod = 2,
   The minimum refined mass for the ByMass refining scheme
   (Usually, user sets OverDensity and code sets MinimumMass but this can be
    overridden by directely setting MinimumMass). 
   The LevelExponent is used to change the minimum mass with level,
   the formula is MinimumMassForRefinement*pow(RefineBy, level*LevelExponent)*/

EXTERN float MinimumOverDensityForRefinement[MAX_FLAGGING_METHODS];
EXTERN float MinimumMassForRefinement[MAX_FLAGGING_METHODS];
EXTERN float MinimumMassForRefinementLevelExponent[MAX_FLAGGING_METHODS];
EXTERN float DepositPositionsParticleSmoothRadius;

/* For CellFlaggingMethod = 3,
   The minimum pressure jump required to be a shock.
   The minimum internal/total energy ratio for a shock. */

EXTERN float MinimumPressureJumpForRefinement, MinimumEnergyRatioForRefinement;

/* For CellFlaggingMethod = 6,
   The number of cells by which the Jeans length should be resolved. */

EXTERN float RefineByJeansLengthSafetyFactor;

/* If > 0, this will be used instead of the temperature at all locations */

EXTERN float JeansRefinementColdTemperature;

/* For CellFlaggingMethod = 8,
   The level to which the must refine particles apply */

EXTERN int   MustRefineParticlesRefineToLevel;

/* For CellFlaggingMethod = 8,
   The physical length (in pc) to which the must refine particles apply 
   The above parameter will be automatically adjusted to match this length */

EXTERN int   MustRefineParticlesRefineToLevelAutoAdjust;

/* For CellFlaggingMethod = 8,
   For new particle system only refine around particles above the minimum mass */

EXTERN float MustRefineParticlesMinimumMass;

/* For CellFlaggingMethod = 9,   
   The minimum shear (roughly, dv accross two zones) required for 
   refinement.    */

EXTERN float MinimumShearForRefinement;

/* For CellFlaggingMethod = 9,   
   Whether to use the old method for calculating shear refinement.    */

EXTERN int OldShearMethod;

/* For CellFlaggingMethod = 11,
   The number of cells by which the Resistive length abs(B)/abs(curl(B)) 
   should be resolved. */

EXTERN float RefineByResistiveLengthSafetyFactor;

/* For CellFlaggingMethod = 14,   
   Minimum mach number required for refinement.    */

EXTERN float ShockwaveRefinementMinMach;
EXTERN float ShockwaveRefinementMinVelocity;
EXTERN int ShockwaveRefinementMaxLevel;

/* For CellFlaggingMethod = 15,   
   Minimum second derivative required for refinement.    */
EXTERN float MinimumSecondDerivativeForRefinement[MAX_FLAGGING_METHODS];
EXTERN int SecondDerivativeFlaggingFields[MAX_FLAGGING_METHODS];
EXTERN float SecondDerivativeEpsilon;

/* Noh problem switch: Upper-Right quadrant or full domain */

EXTERN int NohProblemFullBox;

/* A boolean flag indicating if we are using coordinate comoving with the
   expansion of the universe. */

EXTERN int   ComovingCoordinates;

/* A flag indicating if we are using star particles. */

EXTERN int   StarParticleCreation;
EXTERN int   StarParticleFeedback;
EXTERN int   NumberOfParticleAttributes;
EXTERN int   AddParticleAttributes;
EXTERN int   BigStarFormation;
EXTERN int   BigStarFormationDone;
EXTERN float BigStarSeparation;
EXTERN double SimpleQ;
EXTERN float SimpleRampTime;

/* Set this flag to allow star formation only once per root grid time
   step (at the beginning) and with a SFR proportional to the full
   root grid time step (as in Kravtsov 2004, for example). Currently
   only implemented for H2REG_STAR. */
EXTERN int   StarFormationOncePerRootGridTimeStep;


/* Parameters governing certain time or redshift-dependent actions. */

EXTERN int   TimeActionType[MAX_TIME_ACTIONS];
EXTERN FLOAT TimeActionTime[MAX_TIME_ACTIONS];
EXTERN FLOAT TimeActionRedshift[MAX_TIME_ACTIONS];
EXTERN float TimeActionParameter[MAX_TIME_ACTIONS];

/* Parameters for direct unigrid dumps of entire top grid */

EXTERN char  *CubeDumps[MAX_CUBE_DUMPS];

/* Parameters governing whether tracer particles are on or off. */

EXTERN int   TracerParticleOn;
EXTERN int   TracerParticleOutputVelocity;
EXTERN FLOAT TracerParticleCreationSpacing;
EXTERN FLOAT TracerParticleCreationLeftEdge[MAX_DIMENSION];
EXTERN FLOAT TracerParticleCreationRightEdge[MAX_DIMENSION];

EXTERN int   ParticleTypeInFile;
EXTERN int   OutputParticleTypeGrouping;

EXTERN int   ExternalBoundaryIO;
EXTERN int   ExternalBoundaryTypeIO;
EXTERN int   ExternalBoundaryValueIO;
EXTERN int   ExternalBoundaryField;
EXTERN int   SimpleConstantBoundary;

EXTERN Eint64 TaskMemory[MAX_NUMBER_OF_TASKS];
EXTERN int    TaskMap[MAX_NUMBER_OF_TASKS];

EXTERN double NodeMem[MAX_NUMBER_OF_NODES];
EXTERN int    NodeMap[MAX_NUMBER_OF_NODES];

/* Boolean to control reading when fields and particles get read in. */

EXTERN int LoadGridDataAtStart;

/* Storing the parameter file name for rebuilding the */
/* cpu and grid file names */
EXTERN char PrevParameterFileName[MAX_NAME_LENGTH];

/* MetaData identifier string */
EXTERN char *MetaDataIdentifier;

/* Zhiling Lan's modified code */

#ifdef MPI_INSTRUMENTATION
EXTERN double GlobalCommunication;
EXTERN double RecvComm;
EXTERN double WaitComm;
EXTERN double timer[MAX_COUNTERS];
EXTERN int counter[MAX_COUNTERS];
EXTERN FILE *filePtr;
EXTERN char tracename[MAX_NAME_LENGTH];
EXTERN double starttime, endtime;
EXTERN double Start_Wall_Time, End_Wall_Time, WallTime;
EXTERN int flagging_count, in_count, out_count, moving_count;
EXTERN float flagging_pct, moving_pct;
#endif /* MPI_INSTRUMENTATION */
EXTERN char name[MAX_NAME_LENGTH];
EXTERN FILE *tracePtr;
EXTERN int traceMPI;
#ifdef MEM_TRACE
EXTERN FILE *memtracePtr;
EXTERN int traceMEM;
EXTERN char memtracename[MAX_NAME_LENGTH];
#endif

/* New Movie Data */

EXTERN int MovieDataField[MAX_MOVIE_FIELDS];
EXTERN int MovieSkipTimestep;
EXTERN int Movie3DVolumes;
EXTERN int MovieVertexCentered;
EXTERN char *NewMovieName;
EXTERN int NewMovieDumpNumber;
EXTERN int NewMovieParticleOn;
EXTERN FLOAT *StarParticlesOnProcOnLvl_Position[128][3]; 
EXTERN float *StarParticlesOnProcOnLvl_Velocity[128][3], *StarParticlesOnProcOnLvl_Mass[128];
EXTERN float *StarParticlesOnProcOnLvl_Attr[128][MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
EXTERN int *StarParticlesOnProcOnLvl_Type[128];
EXTERN PINT *StarParticlesOnProcOnLvl_Number[128];

/* Stanford Hydro Solver variables */

/* Hydro parameters */

EXTERN int UseHydro;
EXTERN int Coordinate;
EXTERN int NSpecies;
EXTERN int NColor;
EXTERN float Theta_Limiter;
EXTERN int RKOrder;
EXTERN int UsePhysicalUnit;
EXTERN int iden;
EXTERN int ietot;
EXTERN int ivx;
EXTERN int ivy;
EXTERN int ivz;
EXTERN int iBx;
EXTERN int iBy;
EXTERN int iBz;
EXTERN int iPhi;
EXTERN int ieint;
EXTERN int iD;
EXTERN int iEtot;
EXTERN int iS1;
EXTERN int iS2;
EXTERN int iS3;
EXTERN int iEint;
EXTERN float SmallRho;
EXTERN float SmallP;
EXTERN float SmallEint;
EXTERN float SmallT;
EXTERN float MaximumAlvenSpeed;
EXTERN int NEQ_HYDRO;
EXTERN int NEQ_MHD;
EXTERN int ReconstructionMethod;
EXTERN int PositiveReconstruction;
EXTERN int RiemannSolverFallback;
EXTERN int RiemannSolver;
EXTERN int ConservativeReconstruction;
EXTERN int EOSType;
EXTERN float EOSSoundSpeed;
EXTERN float EOSCriticalDensity;
EXTERN float EOSGamma;
EXTERN float C_h;
EXTERN float C_p;
EXTERN float DivBDampingLength;
EXTERN int UseConstantAcceleration;
EXTERN float ConstantAcceleration[3];
EXTERN float Mu;
EXTERN int ExternalGravity;
EXTERN float StringKick;
EXTERN int StringKickDimension;
EXTERN int UseFloor;
EXTERN int UseViscosity;
EXTERN float ViscosityCoefficient;
EXTERN int UseAmbipolarDiffusion;
EXTERN int UseResistivity;

/* Chemistry & cooling parameters */

EXTERN float CoolingCutOffDensity1;
EXTERN float CoolingCutOffDensity2;
EXTERN float CoolingPowerCutOffDensity1;
EXTERN float CoolingPowerCutOffDensity2;
EXTERN float CoolingCutOffTemperature;

/* Gravity parameters */

EXTERN double HaloMass;
EXTERN float HaloConcentration;
EXTERN float HaloRedshift;
EXTERN double HaloCentralDensity;
EXTERN double HaloVirialRadius;
EXTERN float ExternalGravityConstant;
EXTERN float ExternalGravityDensity;
EXTERN FLOAT ExternalGravityPosition[MAX_DIMENSION];
EXTERN double ExternalGravityRadius;
EXTERN FLOAT ExternalGravityOrientation[MAX_DIMENSION];

/* Poisson Clean */

EXTERN int UseDivergenceCleaning;
EXTERN int DivergenceCleaningBoundaryBuffer;
EXTERN float DivergenceCleaningThreshold;
EXTERN float PoissonApproximationThreshold;
EXTERN int PoissonBoundaryType;



/* Star Particle paramters */

EXTERN int ShiningParticleID;
EXTERN float SinkMergeDistance;
EXTERN float SinkMergeMass;
EXTERN float TotalSinkMass;
EXTERN int StellarWindFeedback;
EXTERN float StellarWindTurnOnMass;
EXTERN float MSStellarWindTurnOnMass;
EXTERN int NBodyDirectSummation;

/* Turbulence simulation parameters */
EXTERN int UseDrivingField;
EXTERN float DrivingEfficiency;

/* Parameters to use CUDA extensions */ 
EXTERN int UseCUDA;

/* End of Stanford block */


/* ran1 initialization flag for star_maker5 */

EXTERN int ran1_init;

/* random number initialization flag */

EXTERN int rand_init;

/* test problem stuff */
EXTERN TestProblemDataType TestProblemData;

/* Memory Limit */

EXTERN long_int MemoryLimit;

/* Staged input */

#ifdef STAGE_INPUT
EXTERN int StageInput;
EXTERN char LocalPath[MAX_LINE_LENGTH];
EXTERN char GlobalPath[MAX_LINE_LENGTH];
#endif

#ifdef USE_PYTHON
EXTERN int NumberOfPythonCalls;
EXTERN int NumberOfPythonTopGridCalls;
EXTERN int NumberOfPythonSubcycleCalls;
EXTERN PyObject *grid_dictionary;
EXTERN PyObject *old_grid_dictionary;
EXTERN PyObject *hierarchy_information;
EXTERN PyObject *yt_parameter_file;
EXTERN PyObject *conversion_factors;
EXTERN PyObject *my_processor;
#endif
/* Multi-species rate equation flag and associated data. */

EXTERN int MetalCooling;
EXTERN char *MetalCoolingTable;
EXTERN int CIECooling;
EXTERN int H2OpticalDepthApproximation;

//   1 - Adaptive ray tracing transfer
//   0 - none
EXTERN int RadiativeTransfer;
EXTERN int RadiativeTransferHydrogenOnly;
#ifdef TRANSFER
EXTERN long *pix2x;
EXTERN long *pix2y;
EXTERN int  *x2pix;
EXTERN int  *y2pix;
EXTERN FLOAT PhotonTime;
EXTERN float dtPhoton;
#include "RadiationSource.h"
#include "RadiativeTransferParameters.h"
EXTERN RadiationSourceEntry *GlobalRadiationSources;
EXTERN SuperSourceEntry *SourceClusteringTree;
EXTERN SuperSourceEntry *OldSourceClusteringTree;
#ifdef MEMORY_POOL
EXTERN MPool::MemoryPool *PhotonMemoryPool;
#endif

/* [0]: Emitted photons
   [1]: escaped past 0.5 RadiativeTransferPhotonEscapeRadius
   [2]:              1.0           -"-
   [3]:              2.0           -"-
*/
EXTERN double EscapedPhotonCount[4];  
EXTERN double TotalEscapedPhotonCount[4];
EXTERN char *PhotonEscapeFilename;
EXTERN int FieldsToInterpolate[MAX_NUMBER_OF_BARYON_FIELDS];

#include "RadiativeTransferSpectrumTable.h"
EXTERN RadiativeTransferSpectrumTableType RadiativeTransferSpectrumTable;

#endif /* TRANSFER  */

EXTERN int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];
EXTERN int LevelSubCycleCount[MAX_DEPTH_OF_HIERARCHY];
EXTERN float dtRebuildHierarchy[MAX_DEPTH_OF_HIERARCHY];
EXTERN float TimeSinceRebuildHierarchy[MAX_DEPTH_OF_HIERARCHY];
EXTERN float dtThisLevelSoFar[MAX_DEPTH_OF_HIERARCHY];
EXTERN float dtThisLevel[MAX_DEPTH_OF_HIERARCHY];

/* RebuildHierarchy on this level every N cycles. */
EXTERN int RebuildHierarchyCycleSkip[MAX_DEPTH_OF_HIERARCHY];
EXTERN int ConductionDynamicRebuildHierarchy;
EXTERN int ConductionDynamicRebuildMinLevel;

/* Coupled radiative transfer, cooling, and rate solver */
EXTERN int RadiativeTransferCoupledRateSolver;


//   2 - FLD radiation transfer only (no ray-tracing at all)
//   1 - FLD radiation transfer (for optically-thin LW radiation)
//   0 - none
EXTERN int RadiativeTransferFLD;


/* Implicit problem decision flag (only 0 through 3 work for now)
      0 => do not use any implicit solver
      1 => use the gFLDProblem module for single-group coupled FLD
      2 => use the FSProb module for free-streaming FLD radiation 
      3 => use the gFLDSplit module for single-group split FLD
      4 => use the MFProb, multi-frequency fully implicit module
      5 => use the MFSplit, multi-frequency split implicit module
*/
EXTERN int ImplicitProblem;

/* Star-Maker emissivity field generator and uv_param used in calculating Geoffrey's Emissivity0 baryon field */

EXTERN int StarMakerEmissivityField;
EXTERN float uv_param;

/* Shearing Boundary Conditions */

EXTERN float AngularVelocity;
EXTERN float VelocityGradient;
EXTERN int ShearingBoundaryDirection;
EXTERN int ShearingVelocityDirection;
EXTERN int ShearingOtherDirection;
EXTERN int UseMHD;
EXTERN FLOAT TopGridDx[MAX_DIMENSION];
EXTERN int ShearingBoxProblemType; // 0 = advecting sphere; 1 = shearing box; 2 = vortex wave ; 3 = stratified

EXTERN float IsothermalSoundSpeed;
EXTERN int RefineByJeansLengthUnits;



EXTERN int MoveParticlesBetweenSiblings;

/* Particle Splitter */

EXTERN int ParticleSplitterIterations;
EXTERN float ParticleSplitterChildrenParticleSeparation;
EXTERN int ParticleSplitterRandomSeed;

/* Magnetic Field Resetter */

EXTERN int ResetMagneticField;
EXTERN float ResetMagneticFieldAmplitude[MAX_DIMENSION];

/* Star Class MBH Particle IO (PARTICLE_TYPE_MBH) */

EXTERN int MBHParticleIO;
EXTERN char *MBHParticleIOFilename;
EXTERN double MBHParticleIOTemp[30][5+MAX_DIMENSION];
EXTERN char *MBHInsertLocationFilename;
EXTERN int OutputWhenJetsHaveNotEjected;

/* Vorticity Calculations */

EXTERN int VelAnyl;
EXTERN int BAnyl;

/* Write Out External Acceleratopm Field */
EXTERN int WriteExternalAccel;

/* Gas drag */
EXTERN int UseGasDrag;
EXTERN float GasDragCoefficient;

EXTERN char current_error[255];

/* Thermal conduction */

EXTERN int IsotropicConduction;  // TRUE OR FALSE
EXTERN int AnisotropicConduction;  // TRUE OR FALSE
EXTERN float IsotropicConductionSpitzerFraction;  // f_Spitzer
EXTERN float AnisotropicConductionSpitzerFraction;  // f_Spitzer
EXTERN float ConductionCourantSafetyNumber;
EXTERN int SpeedOfLightTimeStepLimit; // TRUE OR FALSE

/* SMBH Feedback in galaxy clusters*/
EXTERN int ClusterSMBHFeedback;  // TRUE OR FALSE
EXTERN float ClusterSMBHJetMdot;  // JetMdot in SolarMass/yr 
EXTERN float ClusterSMBHJetVelocity;  // JetVelocity in km/s 
EXTERN float ClusterSMBHJetRadius;  // JetRadius in cellwidth 
EXTERN float ClusterSMBHJetLaunchOffset;  //in cellwidth
EXTERN float ClusterSMBHStartTime;  // in codeunits, usually is InitialTime of restart 
EXTERN float ClusterSMBHTramp;  // in Myr
EXTERN float ClusterSMBHJetOpenAngleRadius;  // in cellwidth 
EXTERN float ClusterSMBHFastJetRadius;  // FastJetRadius in cellwidth 
EXTERN float ClusterSMBHFastJetVelocity;  // FastJetVelocity in km/s 
EXTERN float ClusterSMBHJetEdot;  // Total feedback Edot in 10^44 ergs/s 
EXTERN float ClusterSMBHKineticFraction;  // fraction of kinetic feedback (0-1)
EXTERN float ClusterSMBHJetAngleTheta;  // from 0 to 1/2, in pi
EXTERN float ClusterSMBHJetAnglePhi;  // from 0 to 2, in pi
EXTERN float ClusterSMBHJetPrecessionPeriod;  //in Myr
EXTERN int ClusterSMBHCalculateGasMass;  // TRUE OR FALSE
EXTERN int ClusterSMBHFeedbackSwitch;  // TRUE OR FALSE
EXTERN float ClusterSMBHEnoughColdGas;  // To turn jet on, in SolarMass 
EXTERN float ClusterSMBHAccretionTime;  // Used only when CalculateGasMass=2
EXTERN int ClusterSMBHJetDim;  // Jet dimension
EXTERN float ClusterSMBHAccretionEpsilon;  // Edot=epsilon*Mdot(accreted/removed)*c^2

EXTERN int MHDCT_debug_flag;
EXTERN int MHDCTSlopeLimiter;
EXTERN int MHDCTDualEnergyMethod;
EXTERN int MHDCTPowellSource;
EXTERN int MHDCTUseSpecificEnergy;
EXTERN float FixedTimestep;
EXTERN int WriteBoundary;
EXTERN int WriteAcceleration;
EXTERN int TracerParticlesAddToRestart;// forces addition of tracer particles to already initialized simulations
EXTERN int MHD_ProjectThisFace[3]; //Used for determining face projection/communication needs for 
                                   //face centered fields
EXTERN float CT_AthenaDissipation;
EXTERN int MHD_WriteElectric;
EXTERN float tiny_pressure;
EXTERN int MHD_CT_Method;
EXTERN int MHD_ProjectB;// Should always be FALSE for the evoloution. May be used in initialization.
EXTERN int MHD_ProjectE;// Should always be TRUE for the evoloution
EXTERN int UseMHDCT;
EXTERN int EquationOfState;
EXTERN char *MHDLabel[3];
EXTERN char *MHDcLabel[3];
EXTERN char *MHDUnits[3];
EXTERN char *MHDeLabel[3];
EXTERN char *MHDeUnits[3];

/* For the database */
EXTERN char *DatabaseLocation;
EXTERN int ExtraOutputs[MAX_EXTRA_OUTPUTS];
EXTERN int CorrectParentBoundaryFlux;

/* For EnzoTiming Behavior */
EXTERN int TimingCycleSkip; // Frequency of timing data dumps.

/* For the shock pool boundary method */
EXTERN float ShockPoolAngle;
EXTERN float ShockPoolShockSpeed;
EXTERN float ShockPoolDelay;

EXTERN float ShockPoolDensity;
EXTERN float ShockPoolTotalEnergy;
EXTERN float ShockPoolVelocity[MAX_DIMENSION];

EXTERN float ShockPoolShockDensity;
EXTERN float ShockPoolShockTotalEnergy;
EXTERN float ShockPoolShockVelocity[MAX_DIMENSION];

/* For the galaxy simulation boundary method */
EXTERN int GalaxySimulationRPSWind;
/* GalaxySimulationRPSWind
 *   0 - OFF 
 *   1 - Simple Shock w/ angle and delay 
 *   2 - Lookup table of density and velocity
 */
EXTERN float GalaxySimulationRPSWindShockSpeed;
EXTERN float GalaxySimulationRPSWindDelay;

EXTERN float GalaxySimulationRPSWindDensity;
EXTERN float GalaxySimulationRPSWindTotalEnergy;
EXTERN float GalaxySimulationRPSWindVelocity[MAX_DIMENSION];
EXTERN float GalaxySimulationRPSWindPressure;

EXTERN float GalaxySimulationPreWindDensity;
EXTERN float GalaxySimulationPreWindTotalEnergy; 
EXTERN float GalaxySimulationPreWindVelocity[MAX_DIMENSION];
 
#endif
