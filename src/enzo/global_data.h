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
#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* Load Balancing.  Currently only memory count method implemented
                          0 = off
                          1 = Equalize processor memory count
                         2 = Load balance only on a node
*/
EXTERN int LoadBalancing;
EXTERN int LoadBalancingCycleSkip;
EXTERN int ResetLoadBalancing;
EXTERN int CoresPerNode;
EXTERN int PreviousMaxTask;

/* FileDirectedOutput checks for file existence: 
   stopNow (writes, stops),   outputNow, subgridcycleCount */
EXTERN int FileDirectedOutput;


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
EXTERN int ProblemType;

/* Hydrodynamics method:
       0 - PPM_DE      1 - PPM_LR (not working)    2 - ZEUS        */

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
                         12 = FlagCellsToBeRefinedByMustRefineRegion
			 13 = FlagCellsToBeRefinedByMetallicity
 */

EXTERN int CellFlaggingMethod[MAX_FLAGGING_METHODS];

/* left and right boundaries of the 'must refine region'
   for CellFlaggingMethod = 10 */

EXTERN FLOAT MustRefineRegionLeftEdge[MAX_DIMENSION];  // left edge

EXTERN FLOAT MustRefineRegionRightEdge[MAX_DIMENSION];  // right edge

/* specifies the level to which FlagCellsToBeRefinedByMustRefineRegion
   will refine up to (does not prevent refinement to higher levels) */

EXTERN int MustRefineRegionMinRefinementLevel;

/* specifies the level to which FlagGridCellsToBeRefinedByMetallicity
   will refine up to (does not prevent refinement to higher levels) */
EXTERN int MetallicityRefinementMinLevel;

/* threshold metallicity for FlagGridCellsToBeRefinedByMetallicity */
EXTERN float MetallicityRefinementMinMetallicity;


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

/* Uniform gravity: on/off flag, direction, and strength. */

EXTERN int UniformGravity, UniformGravityDirection;
EXTERN float UniformGravityConstant;

/* point source gravity: on/off flag position, and strength. */

EXTERN int PointSourceGravity;
EXTERN FLOAT PointSourceGravityPosition[MAX_DIMENSION];
EXTERN float PointSourceGravityConstant;
EXTERN float PointSourceGravityCoreRadius;

/* SelfGravity (TRUE or FALSE) */

EXTERN int SelfGravity;

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

/* Cloudy cooling parameters and data. */

EXTERN CloudyCoolingDataType CloudyCoolingData;

/* Gadget Equilibrium cooling on/off flag */

EXTERN int GadgetEquilibriumCooling;

/* Random Forcing on/off flag and associated data. */ //AK

EXTERN int     RandomForcing;
EXTERN FLOAT   RandomForcingEdot;
EXTERN FLOAT   RandomForcingMachNumber;
EXTERN fpos_t  BaryonFileNamePosition;

/* Multi-species rate equation flag and associated data. */

EXTERN int MultiSpecies;
EXTERN RateDataType RateData;

/* Glover chemistry/cooling network flags */
EXTERN int GloverChemistryModel;  // 0 is off, on is 1-7, excluding 6
EXTERN int GloverRadiationBackground; // 1: low Z, 2: ISM
EXTERN int GloverOpticalDepth; // 0: opticaly thin, 1: single-cell

/* Multi-element metallicity field flag and count. */

EXTERN int MultiMetals;

/* Type of radiation field. 
   0 - none,                    1 - Haardt & Madau alpha=-1.5
   2 - H&M alpha = -1.8       
   10 - homogenous internal radiation field (a la Renyue's work) */

EXTERN int RadiationFieldType;
EXTERN int AdjustUVBackground; 
EXTERN float SetUVBAmplitude;
EXTERN float SetHeIIHeatingScale;
EXTERN RadiationFieldDataType RadiationData;
EXTERN int RadiationFieldLevelRecompute;
EXTERN int RadiationXRaySecondaryIon;

/* Output cooling time with grid data. */

EXTERN int OutputCoolingTime;

/* Output temperature with grid data. */

EXTERN int OutputTemperature;

/* Output smoothed dark matter fields. */

EXTERN int OutputSmoothedDarkMatter;
EXTERN int SmoothedDarkMatterNeighbors;

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
EXTERN int ExtractFieldsOnly;
EXTERN int First_Pass;
EXTERN int UnigridTranspose;
EXTERN int NumberOfRootGridTilesPerDimensionPerProcessor;

/* Parameter(s) for embedded python execution */
EXTERN int PythonSubcycleSkip;

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

/* For CellFlaggingMethod = 8,
   The level to which the must refine particles apply */

EXTERN int   MustRefineParticlesRefineToLevel;

/* For CellFlaggingMethod = 9,   
   The minimum shear (roughly, dv accross two zones) required for 
   refinement.    */



EXTERN float MinimumShearForRefinement;

/* For CellFlaggingMethod = 11,
   The number of cells by which the Resistive length abs(B)/abs(curl(B)) 
   should be resolved. */

EXTERN float RefineByResistiveLengthSafetyFactor;


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

/* Parameters governing certain time or redshift-dependent actions. */

EXTERN int   TimeActionType[MAX_TIME_ACTIONS];
EXTERN FLOAT TimeActionTime[MAX_TIME_ACTIONS];
EXTERN FLOAT TimeActionRedshift[MAX_TIME_ACTIONS];
EXTERN float TimeActionParameter[MAX_TIME_ACTIONS];

/* Parameters for direct unigrid dumps of entire top grid */

EXTERN char  *CubeDumps[MAX_CUBE_DUMPS];

/* Parameters governing whether tracer particles are on or off. */

EXTERN int   TracerParticleOn;
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

/* Zhiling Lan's modified code */

#ifdef MPI_INSTRUMENTATION
EXTERN double GlobalCommunication;
EXTERN double RecvComm;
EXTERN double WaitComm;
EXTERN double timer[MAX_COUNTERS];
EXTERN int counter[MAX_COUNTERS];
EXTERN FILE *filePtr;
EXTERN char tracename[MAX_NAME_LENGTH];
EXTERN char memtracename[MAX_NAME_LENGTH];
EXTERN FILE *memtracePtr;
EXTERN int traceMEM;
EXTERN double starttime, endtime;
EXTERN double Start_Wall_Time, End_Wall_Time, WallTime;
EXTERN int flagging_count, in_count, out_count, moving_count;
EXTERN float flagging_pct, moving_pct;
#endif /* MPI_INSTRUMENTATION */
EXTERN char name[MAX_NAME_LENGTH];
EXTERN FILE *tracePtr;
EXTERN int traceMPI;

/* New Movie Data */

EXTERN int MovieDataField[MAX_MOVIE_FIELDS];
EXTERN int MovieSkipTimestep;
EXTERN int Movie3DVolumes;
EXTERN int MovieVertexCentered;
EXTERN char *NewMovieName;
EXTERN int NewMovieDumpNumber;
EXTERN int NewMovieParticleOn;

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
EXTERN int RiemannSolver;
EXTERN int EOSType;
EXTERN float EOSSoundSpeed;
EXTERN float EOSCriticalDensity;
EXTERN float EOSGamma;
EXTERN float C_h;
EXTERN float C_p;
EXTERN int UseConstantAcceleration;
EXTERN float ConstantAcceleration[3];
EXTERN float Mu;
EXTERN int ExternalGravity;
EXTERN int StringKick;
EXTERN int UseFloor;
EXTERN int UseViscosity;
EXTERN int UseAmbipolarDiffusion;
EXTERN int UseResistivity;

/* Chemistry & cooling parameters */

EXTERN int UseH2OnDust;
EXTERN double PhotoelectricHeating;
EXTERN float CoolingCutOffDensity1;
EXTERN float CoolingCutOffDensity2;
EXTERN float CoolingPowerCutOffDensity1;
EXTERN float CoolingPowerCutOffDensity2;
EXTERN float CoolingCutOffTemperature;
EXTERN int CoolingModel;

/* Gravity parameters */

EXTERN double HaloMass;
EXTERN float HaloConcentration;
EXTERN float HaloRedshift;
EXTERN double HaloCentralDensity;
EXTERN double HaloVirialRadius;
EXTERN float ExternalGravityDensity;
EXTERN double ExternalGravityRadius;

/* Poisson Clean */

EXTERN int UseDivergenceCleaning;
EXTERN int DivergenceCleaningBoundaryBuffer;
EXTERN float DivergenceCleaningThreshold;
EXTERN float PoissonApproximationThreshold;



/* Star Particle paramters */

EXTERN int ShiningParticleID;
EXTERN double SinkMergeDistance;
EXTERN float SinkMergeMass;
EXTERN float TotalSinkMass;
EXTERN int StellarWindFeedback;
EXTERN float StellarWindTurnOnMass;
EXTERN int NBodyDirectSummation;

/* Turbulence simulation parameters */
EXTERN int UseDrivingField;
EXTERN float DrivingEfficiency;
/* Parameters to use CUDA extensions */ 
EXTERN int UseCUDA;

/* End of Stanford block */


/* ran1 initialization flag for star_maker5 */

EXTERN int ran1_init;

/* test problem stuff */
EXTERN TestProblemDataType TestProblemData;

/* Memory Limit */

EXTERN int MemoryLimit;

/* Staged input */

#ifdef STAGE_INPUT
EXTERN int StageInput;
EXTERN char LocalPath[MAX_LINE_LENGTH];
EXTERN char GlobalPath[MAX_LINE_LENGTH];
#endif

#ifdef USE_PYTHON
EXTERN int NumberOfPythonCalls;
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

//   1 - Adaptive ray tacing transfer
//   0 - none
EXTERN int RadiativeTransfer;
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

/* [0]: Emitted photons
   [1]: escaped past 0.5 RadiativeTransferPhotonEscapeRadius
   [2]:              1.0           -"-
   [3]:              2.0           -"-
*/
EXTERN double EscapedPhotonCount[4];  
EXTERN double TotalEscapedPhotonCount[4];
EXTERN char *PhotonEscapeFilename;
EXTERN int FieldsToInterpolate[MAX_NUMBER_OF_BARYON_FIELDS];

#endif /* TRANSFER  */

/* Coupled radiative transfer, cooling, and rate solver */

EXTERN int RadiativeTransferCoupledRateSolver;



/* Shearing Boundary Conditions */

EXTERN float AngularVelocity;
EXTERN float VelocityGradient;
EXTERN int ShearingBoundaryDirection;
EXTERN int ShearingVelocityDirection;
EXTERN int ShearingOtherDirection;
EXTERN int useMHD;
EXTERN FLOAT TopGridDx[MAX_DIMENSION];
EXTERN int ShearingBoxProblemType; // 0 = advecting sphere; 1 = shearing box; 2 = vortex wave ; 3 = stratified

EXTERN float IsothermalSoundSpeed;
EXTERN int RefineByJeansLengthUnits;

EXTERN int MoveParticlesBetweenSiblings;

#endif
