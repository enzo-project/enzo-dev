/***********************************************************************
/
/  Cen & Ostriker star formation
/
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

class ActiveParticleType_CenOstriker : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_CenOstriker(void) : ActiveParticleType() {};
  ActiveParticleType_CenOstriker(ActiveParticleType* part) : 
    ActiveParticleType(static_cast<ActiveParticleType*>(part)) {};
  ActiveParticleType* clone() {
    return static_cast<ActiveParticleType*>(
        new ActiveParticleType_CenOstriker(this)
      );
  }
  static int EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static int EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, bool CallEvolvePhotons, 
				 int TotalActiveParticleCountPrevious[],
				 int CenOstrikerID) { return SUCCESS; };
  template <class active_particle_class>
    static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalActiveParticleCountPrevious[],
				int CenOstrikerID) { return SUCCESS; };
  template <class active_particle_class>
    static int DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int GalaxyParticleID) {return SUCCESS; };
  static int SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, int TopGridDims[], int ActiveParticleID);
  static int ResetAcceleration(float *ActiveParticleAcceleration);
  static int CreateParticle(grid *thisgrid_orig, ActiveParticleFormationData &supp_data,
			    int particle_index);
  static int InitializeParticleType();
  ENABLED_PARTICLE_ID_ACCESSOR
  
  static float OverdensityThreshold, MassEfficiency, MinimumDynamicalTime, 
    MinimumStarMass, MassEjectionFraction, EnergyToThermalFeedback, MetalYield;

  static int FeedbackDistTotalCells, FeedbackDistRadius, FeedbackDistCellStep;

  static bool JeansMassCriterion, StochasticStarFormation, UnigridVelocities, 
    PhysicalOverdensity, dtDependence;

  static std::vector<ParticleAttributeHandler *> AttributeHandlers;
};

