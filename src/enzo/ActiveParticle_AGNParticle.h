/***********************************************************************
/
/  A "skeleton" active particle that compiles but doesn't do much.
/
************************************************************************/
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "EventHooks.h"
#include "ActiveParticle.h"
#include "TopGridData.h"
#include "communication.h"
#include "CommunicationUtilities.h"
#include "phys_constants.h"

class ActiveParticleType_AGNParticle : public ActiveParticleType
{
public:

  // Default constructor (no arguments)
  // All AGNParticle attributes are set to their default values.
  ActiveParticleType_AGNParticle(void) : ActiveParticleType()
    {
      // Initialize any instance (i.e. not static) member variables here
      CoolingRadius        = 1.0;
      FeedbackRadius       = 1.0;
      CondensationFraction = 0.1;
      FeedbackEfficiency   = 1.0e-3;
      TimescaleThreshold   = 10.0;
      KineticFraction      = 1.0;
      JetAngle             = (10.0/360.0) * 2.0 * M_PI;
      JetPhi               = (20.0/360.0) * 2.0 * M_PI;
      JetTheta             = 2.0 * M_PI;
      StoredEnergy         = 0.0;
      StoredMass           = 0.0;
      FixedJetIsOn         = 0;
      TimeOfLastShock      = 0.0;
    };

  // copy constructor
  ActiveParticleType_AGNParticle(ActiveParticleType_AGNParticle* part) :
    ActiveParticleType(static_cast<ActiveParticleType*>(part))
    {
      CondensationFraction = part->CondensationFraction;
      FeedbackEfficiency   = part -> FeedbackEfficiency;
      JetAngle             = part -> JetAngle;
      TimescaleThreshold   = part -> TimescaleThreshold;
      KineticFraction      = part -> KineticFraction;
      CoolingRadius        = part -> CoolingRadius;
      FeedbackRadius       = part -> FeedbackRadius;
      JetPhi               = part -> JetPhi;
      JetTheta             = part -> JetTheta;
      Edot                 = part -> Edot;
      StoredEnergy         = part -> StoredEnergy;
      StoredMass           = part -> StoredMass;
      FixedJetIsOn         = part -> FixedJetIsOn;
      TimeOfLastShock      = part -> TimeOfLastShock;

      pos[0] = part->pos[0];
      pos[1] = part->pos[1];
      pos[2] = part->pos[2];

      vel[0] = part->vel[0];
      vel[1] = part->vel[1];
      vel[2] = part->vel[2];

      Mass = part->Mass;
      BirthTime = part->BirthTime;
      level = part->level;
    };
  
  // Needed to Create a copy of a particle when only a pointer to the base
  // class is available.
  ActiveParticleType* clone() 
    {
      return static_cast<ActiveParticleType*>(
          new ActiveParticleType_AGNParticle(this)
        );
    }
   // AGNParticle data members
   float CondensationFraction;
   float FeedbackEfficiency;
   float JetAngle;
   float TimescaleThreshold;
   float KineticFraction;
   float CoolingRadius;
   float FeedbackRadius;
   float JetPhi;
   float JetTheta;
   float Edot;

   int FixedJetIsOn;

   float TimeOfLastShock;

   float StoredEnergy;
   float StoredMass;

   float AccretionHistoryTime[256];
   float AccretionHistoryRate[256];

   // These two need to be static as they are used to construct the feedback
   // zone
   static float static_cooling_radius;
   static float static_feedback_radius;


  static bool AGNInsideGrid(grid*);
  static int CreateAGN(grid*, ActiveParticleFormationData&);
  /*
   * Force an AGN particle to spawn in the highest refined grid at the location
   * and time specified in the parameter file.
   * This is called inside the grid::ActiveParticleHandler function, which is in
   *   turn called just after the grids take a hydro timestep in EvolveLevel.
   */
  static int InsertAGN(grid *thisgrid_orig, HierarchyEntry* SubgridPointer, ActiveParticleFormationData &data);

  /*
   * Run an algorithm to determine whether a particle forms in a grid.
   *
   * This is called inside the grid::ActiveParticleHandler function, which is in
   *   turn called just after the grids take a hydro timestep in EvolveLevel.
   */
  static int EvaluateFormation(grid *thisgrid_orig, TopGridData *MetaData, ActiveParticleFormationData &data);


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
   * This funcion and the corresponding AfterEvolveLevel function are the best
   * place to implement any feedback or subgrid physics algorithm that needs
   * information outside of the grid a particle lives on.  These functions
   * have access to the full AMR hierarchy and can arbitrarily change the
   * state of the simulation.
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
      int AGNParticleID);

  /*
   * Perform distributed operations on the hierarchy.  This function can
   * access the entire grid hierarchy and can in principle make arbitrary
   * modifications to any grid.
   *
   * This function and the corresponding BeforeEvolveLevel function are the best
   * place to implement any feedback or subgrid physics algorithm that needs
   * information outside of the grid a particle lives on.  These functions
   * have access to the full AMR hierarchy and can arbitrarily change the
   * state of the simulation.
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
				int AGNParticleID);
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
				int ThisLevel, int AGNParticleID) {return SUCCESS; };

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
   * This function allows for the particle acceleration to be set to zero.
   * Useful for fixing the position of an active particle.
   */
  static int ResetAcceleration(float *ActiveParticleAcceleration);

  /*
   * This function is a brute force method to create this particle type.
   */
  static int CreateParticle(grid *thisgrid_orig, ActiveParticleFormationData &supp_data,
			    int particle_index);
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
  /*
   * This must appear unchanged in all particle types.  It is a macro that
   * returns a unique ID for the particle type.  The ID may be different from
   * simulation to simulation.
   */
   ENABLED_PARTICLE_ID_ACCESSOR

   static int Handle_AGN_Feedback(int nParticles, ActiveParticleList<ActiveParticleType>& ParticleList,
                      FLOAT dx, LevelHierarchyEntry *LevelArray[], int ThisLevel);
  /*
   * The AttributeHandler is used to save active particles to output files and
   * communitcate active particles over MPI.  Instance variables need to be
   * registered with the attribute handler.
   */
  static std::vector<ParticleAttributeHandler *> AttributeHandlers;
};

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
                      HierarchyEntry **Grids[]);
int AssignActiveParticlesToGrids(
    ActiveParticleList<ActiveParticleType>& ParticleList, int nParticles,
    LevelHierarchyEntry *LevelArray[]);
/* Before EvolveLevel
 * Gets called before Evolve level.
 * AGN feedback is being here.
 * */
template <class active_particle_class>
int ActiveParticleType_AGNParticle::BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
        int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
        int ThisLevel, bool CallEvolvePhotons, int TotalActiveParticleCountPrevious[],
        int AGNParticleID)
{
   if (debug)
      printf("Entering BeforeEvolveLevel [%"ISYM"]\n", MyProcessorNumber);

   int i, nParticles, nMergedParticles;
   ActiveParticleList<ActiveParticleType> ParticleList;

   ActiveParticleFindAll(LevelArray, &nParticles, AGNParticleID, ParticleList);
    TotalAGNParticlesCreated = nParticles;
    
    if (nParticles == 0)
        return SUCCESS;

   if (ThisLevel == MaximumRefinementLevel)
   {  
       //ActiveParticleList<ActiveParticleType> ParticleList; 
       // Get a list of AGN particles
       

       // Calculate the cell size
       // This assumes a cubic box and may not work for simulations with MinimumMassForRefinementLevelExponent
       FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
                  (MetaData->TopGridDims[0] * POW(FLOAT(RefineBy), FLOAT(MaximumRefinementLevel)));
       // Calls Do AGN feedback
       if (Handle_AGN_Feedback(nParticles, ParticleList, dx, LevelArray, ThisLevel) == FAIL)
           ENZO_FAIL("AGN Particle feedback failed. \n");

       // Get rid of particles that are no longer on this processor
       for (i = 0; i < nParticles; i++)
           if (ParticleList[i] != NULL) {
              if (ParticleList[i]->ReturnCurrentGrid()->ReturnProcessorNumber() != MyProcessorNumber) {
                  delete ParticleList[i];
                  ParticleList[i] = NULL;
               }
	   }

      //delete [] ParticleList;
      //ParticleList.clear();
   }
   ParticleList.clear();
   if (debug)
      printf("Leaving BeforeEvolveLevel [%"ISYM"]\n", MyProcessorNumber);

   return SUCCESS;
}

/* After EvolveLevel
 * Gets called after EvolveLevel
 * Currently no action is being taken.
 *    */
template <class active_particle_class>
int ActiveParticleType_AGNParticle::AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
        int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
        int ThisLevel, int TotalStarParticleCountPrevious[],
        int AGNParticleID)
{
   //if (debug)
   //   printf("Entering AfterEvolveLevel [%"ISYM"]\n", MyProcessorNumber);


   // All finished
   //if (debug)
   //   printf("Leaving AfterEvolveLevel [%"ISYM"]\n", MyProcessorNumber);
   
   return SUCCESS;
   }
