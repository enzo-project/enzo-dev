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
#include "phys_constants.h"

class ActiveParticleType_Kravtsov;

class ActiveParticleType_Kravtsov : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_Kravtsov(void) : ActiveParticleType() {};
  ActiveParticleType_Kravtsov(ActiveParticleType_Kravtsov* part) :
    ActiveParticleType(static_cast<ActiveParticleType*>(part)) {};

  ActiveParticleType* clone() {
    return static_cast<ActiveParticleType*>(
        new ActiveParticleType_Kravtsov(this)
      );
  }
  // Static members
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &supp_data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, bool CallEvolvePhotons,
				 int TotalStarParticleCountPrevious[],
				 int SampleParticleID) { return SUCCESS; };
  template <class active_particle_class>
    static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int SampleParticleID) { return SUCCESS; };
  template <class active_particle_class>
    static int DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int GalaxyParticleID) {return SUCCESS; };
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID);
  static int ResetAcceleration(float *ActiveParticleAcceleration);
  static int CreateParticle(grid *thisgrid_orig, ActiveParticleFormationData &supp_data,
			    int particle_index);
  static int InitializeParticleType(void);
  ENABLED_PARTICLE_ID_ACCESSOR

  static float DensityThreshold, StarFormationTimeConstant, MinimumStarMass;

  static std::vector<ParticleAttributeHandler *> AttributeHandlers;

};
