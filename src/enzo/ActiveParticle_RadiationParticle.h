/***********************************************************************
/
/  Header File for the Radiation Particle. 
/  Based on template provided by Nathan Goldbaum - ActiveParticle_Skeleton.h
/
/  Author: John Regan
/
/  Date: July 2014
************************************************************************/

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EventHooks.h"
#include "ActiveParticle.h"

#define DEBUG 0

/* Structures */
struct InitData {   //LL used to store Radiation Particle Data
  FLOAT Position[3];
  float Redshift;
  float Redshift_end;
  bool Alive;
  bool Create;
  bool Deleteme;
  InitData *next;
};
/* Function prototypes */
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
bool CheckForParticleAction(InitData *Root);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);



class ActiveParticleType_RadiationParticle : public ActiveParticleType
{
public:
  // void constructor (no arguments)
  ActiveParticleType_RadiationParticle(void) : ActiveParticleType()
    {
      
    };

  // copy constructor
  ActiveParticleType_RadiationParticle(ActiveParticleType_RadiationParticle* part) :
    ActiveParticleType(static_cast<ActiveParticleType*>(part))
    {
      // Copy values of instance members using data from the particle instance
      // that is passed as an argument to this function
    };

  // Needed to Create a copy of a particle when only a pointer to the base
  // class is available.
  ActiveParticleType* clone() 
  {
    return static_cast<ActiveParticleType*>(
	   new ActiveParticleType_RadiationParticle(this)
					   );
  };
  /*
   * Run an algorithm to determine whether a particle forms in a grid.
   *
   * This is called inside the grid::ActiveParticleHandler function, which is in
   *   turn called just after the grids take a hydro timestep in EvolveLevel.
   */
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);


  /*
   * Run an algorithm to do feedback local to the grid the particle lives in.
   *
   * This is called inside the grid::ActiveParticleHandler function, which is in
   * turn called just after the grids take a hydro timestep in EvolveLevel.
   */
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);


  /*
   * Set boolean flags to determine which supplemental data is passed to
   * EvaluateFormation.
   *
   * This is called inside the grid::ActiveParticleHandler function, which is in
   * turn called just after the grids take a hydro timestep in EvolveLevel.
   */
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);

  /*
   * Perform distributed operations on active particles.  This function can
   * access the entire grid hierarchy and can in principle make arbitrary
   * modifications to any grid.
   *
   * This function is called in the ActiveParticleInitialize function, which is
   * in turn called before the loop that advances the grids by a timestep
   * happens in EvolveLevel.
   *
   * Note that this is a template function so it must be implemented in a header
   * file.
   */
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
				 int ThisLevel, bool CallEvolvePhotons,
				 int TotalActiveParticleCountPrevious[],
				 int RadiationParticleID);

  /*
   * Perform distributed operations on the hierarchy.  This function can
   * access the entire grid hierarchy and can in principle make arbitrary
   * modifications to any grid.
   *
   * This function is called in the ActiveParticleFinalize function, which is
   * in turn called after the loop that advances the grids by a timestep
   * happens in EvolveLevel.
   *
   * Note that this is a template function so it must be implemented in a header
   * file.
   */
  template <class active_particle_class>
  static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
				int ThisLevel, int TotalActiveParticleCountPrevious[],
				int RadiationParticleID)
    {

       if(DEBUG)
	 {
	   FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
	   float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
	     VelocityUnits;
	   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		    &TimeUnits, &VelocityUnits, Time);
	   int nParticles = 0, ipart = 0;
           ActiveParticleList<ActiveParticleType> RadiationParticleList;
	   ActiveParticleFindAll(LevelArray, &nParticles, RadiationParticleID, RadiationParticleList);
	   
	   ActiveParticleType_RadiationParticle *ThisParticle;
	   for (ipart = 0; ipart < nParticles; ipart++) {
	     if (RadiationParticleList[ipart]->IsARadiationSource(Time)) {
	       ThisParticle =
		 static_cast<ActiveParticleType_RadiationParticle*>(RadiationParticleList[ipart]);
	       ThisParticle->PrintInfo();
	     }
	   }
	 }
       
       return SUCCESS;
    };

  /*
   * This function allows fine-grained control over how an active particle is
   * deposited onto the GravitatingMassField.  If your particle represents an
   * extended object, this might be a better way to handle its self-gravity
   * instead of treating the particle as a point mass.
   *
   * Note that the particle mass will still be deposited in
   * grid::DepositParticlePositionsLocal.
   *
   * This is called inside the DepositActiveParticleMassFlaggingField function,
   * which is itself called inside the RebuildHierarchy function.
   */
  template <class active_particle_class>
    static int DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
				int ThisLevel, int GalaxyParticleID) {return SUCCESS; };

  /*
   * This function allows fine-grained control over grid refinement.  If for
   * some reason your particle type needs to refine to a certain level in its
   * vicinity, this is the place to control that.
   *
   * This is called inside the DepositActiveParticleMassFlaggingField function
   * which is itself called inside RebuildHierarchy.
   */
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level,
      int TopGridDims[], int ActiveParticleID);

  /*
   * Register class-level metadata about your particle type.  This includes
   * reading in parameter values, any class-level activities that only happen
   * once at the beginning of the simulation, and registering instance members
   * with the AttributeHandler.
   *
   * This is called inside the EnableActiveParticleType function which is itself
   * called inside ReadParameterFile.
   */
  static int InitializeParticleType();
  bool IsARadiationSource(FLOAT Time);
  static void SetRadiationDefaults();
  static int* GetGridIndices(FLOAT *position, FLOAT *GridLeftEdge, float CellWidth);
  static int ReadRadiationParameterFile();

  static int ResetAcceleration(float *ActiveParticleAcceleration);
  static int CreateParticle(grid *thisgrid_orig, ActiveParticleFormationData &supp_data,
			    int particle_index);
  /*
   * This must appear unchanged in all particle types.  It is a macro that
   * returns a unique ID for the particle type.  The ID may be different from
   * simulation to simulation.
   */
  ENABLED_PARTICLE_ID_ACCESSOR;

  /*
   * Static variables should be defined here.  Since they are static, there is
   * only one copy of these variables for all instances of the class.  This is
   * useful for storing the value of a runtime parameter, for example.
   */
  static int RadiationNumSources, RadiationSEDNumberOfBins;
  static float *RadiationEnergyBins, *RadiationSED, RadiationLifetime; 
  static float RadiationPhotonsPerSecond;
  static InitData *Root;
  static bool FixedInSpace;
  /*
   * Instance variables should be defined here. Each instance of the particle
   * class will have its own version of these variables.
   */
  
  /*
   * The AttributeHandler is used to save active particles to output files and
   * communitcate active particles over MPI.  Instance variables need to be
   * registered with the attribute handler.
   */
  static std::vector<ParticleAttributeHandler *> AttributeHandlers;
};
