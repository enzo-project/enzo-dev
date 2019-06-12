/***********************************************************************
/
/ Smart Star Particle
/ Author: John Regan
/         Based on previous active particle code from many authors
/ Date: Early 2017
/
************************************************************************/

#ifdef USE_MPI
#endif 

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
#include "CommunicationUtilities.h"
#include "communication.h"
#include "phys_constants.h"
#include "FofLib.h"
#define DEBUG 0
/* Every how many times will the accretion rate be updated */
#define FREQUENCY 100
#define MAXACCRETIONRADIUS  128 /* Times the minimum cell width */
#define ACCRETIONRADIUS  4
#define NUMRADIATIONBINS 5
#define CRITICAL_ACCRETION_RATE 0.04 //Msolar/yr
#define TIMEGAP            900   //yrs
/* Prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
static FLOAT Dist(FLOAT *tempPos, FLOAT *tempPos1);
static float R_ISCO(float a);
//Particle Classes
#define POPIII 0
#define SMS    1
#define BH     2

//Accretion Modes
#define SPHERICAL_BONDI_HOYLE_FORMALISM 1
#define ANGULAR_MOMENTUM_LIMITED_ACCRETION 4
#define SPHERICAL_BONDI_HOYLE_FORMALISM_WITH_VORTICITY 5
#define VISCOUS_ANGULAR_MOMENTUM_TRANSPORT 6
#define ALPHA_DISK_CEN_2012 7
#define CONVERGING_MASS_FLOW 8

class ActiveParticleType_SmartStar : public ActiveParticleType
{
public:
  // Constructors
  ActiveParticleType_SmartStar(void) : ActiveParticleType() {
    AccretionRadius = -1;
    ParticleClass = 0;
    NotEjectedMass = 0;
    MassToBeEjected = 0;
    eta_disk = 1 - sqrt(1 - 2.0/(3.0*R_ISCO(SmartStarSpin)));
    mass_in_accretion_sphere = 0;
    epsilon_deltat = 1.0;
    beta_jet = 0.0;
    for(int i = 0; i < 3; i++)
      {
	Accreted_angmom[i] = 0.0;
      }
    TimeIndex = -1;
    oldmass = -1;
    for(int i = 0; i < NTIMES; i++)
      {
	AccretionRateTime[i] = -222222;
	AccretionRate[i] = -222222;
      }
  };
  ActiveParticleType_SmartStar(ActiveParticleType_SmartStar* part) :
    ActiveParticleType(static_cast<ActiveParticleType*>(part)) {   
    AccretionRadius = part->AccretionRadius;
    ParticleClass = part->ParticleClass;
    for(int i = 0; i < 3; i++) {
      Accreted_angmom[i] = 0.0;
    }
   
    TimeIndex = part->TimeIndex;
    oldmass = part->oldmass;
    for(int i = 0; i < NTIMES; i++) {
      AccretionRateTime[i] = part->AccretionRateTime[i];
      AccretionRate[i] = part->AccretionRate[i];
    }
    NotEjectedMass = part->NotEjectedMass;
    MassToBeEjected = part->MassToBeEjected;
    eta_disk = 1 - sqrt(1 - 2.0/(3.0*R_ISCO(SmartStarSpin)));
    epsilon_deltat = part->epsilon_deltat;
    beta_jet = part->beta_jet;
    mass_in_accretion_sphere = part->mass_in_accretion_sphere;
    };
  ActiveParticleType* clone() {
    return static_cast<ActiveParticleType*>(
        new ActiveParticleType_SmartStar(this)
      );
  };
  static int EvaluateFormation(
      grid *thisgrid_orig, ActiveParticleFormationData &data);
  static void DescribeSupplementalData(ActiveParticleFormationDataFlags &flags);
  static int EvaluateFeedback(
      grid *thisgrid_orig, ActiveParticleFormationData &data);
  template <class active_particle_class>
    static int BeforeEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				 int ThisLevel, bool CallEvolvePhotons,
				 int TotalStarParticleCountPrevious[],
				 int SmartStarID);
  template <class active_particle_class>
    static int AfterEvolveLevel(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int TotalStarParticleCountPrevious[],
				int SmartStarID);
  template <class active_particle_class>
    static int DepositMass(HierarchyEntry *Grids[], TopGridData *MetaData,
				int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
				int ThisLevel, int GalaxyParticleID) {return SUCCESS; };
  static int SetFlaggingField(
      LevelHierarchyEntry *LevelArray[], 
      int level, 
      int TopGridDims[], 
      int ActiveParticleID);

  static int SmartStarParticleFeedback(
             int nParticles, 
             ActiveParticleList<ActiveParticleType>& ParticleList,
		     FLOAT dx, LevelHierarchyEntry *LevelArray[], int ThisLevel);
  

  static int ResetAcceleration(float *ActiveParticleAcceleration);
  static int CreateParticle(grid *thisgrid_orig, ActiveParticleFormationData &supp_data,
			    int particle_index);
  static int InitializeParticleType();
  void   SmartMerge(ActiveParticleType_SmartStar *a);
  int CalculateAccretedAngularMomentum();
  int SmartStarAddFeedbackSphere();
  ENABLED_PARTICLE_ID_ACCESSOR
  bool IsARadiationSource(FLOAT Time);
  
  // sink helper routines

  template <class active_particle_class>
  static void MergeSmartStars(
      int *nParticles, ActiveParticleList<ActiveParticleType>& ParticleList, 
      int *ngroups, 
      LevelHierarchyEntry *LevelArray[], int ThisLevel, 
      ActiveParticleList<active_particle_class>& MergedParticles);

  static int Accrete(int nParticles, 
      ActiveParticleList<ActiveParticleType>& ParticleList,
      FLOAT AccretionRadius, FLOAT dx, 
      LevelHierarchyEntry *LevelArray[], int ThisLevel);

  static int UpdateAccretionRateStats(int nParticles,
				      ActiveParticleList<ActiveParticleType>& ParticleList,
				      FLOAT dx, 
				      LevelHierarchyEntry *LevelArray[], int ThisLevel);

 
  static float EjectedMassThreshold;
  FLOAT AccretionRadius;   // in units of CellWidth on the maximum refinement level
  static int RadiationParticle;
  static double LuminosityPerSolarMass;
  static int RadiationSEDNumberOfBins;
  static float* RadiationEnergyBins;
  static float* RadiationSED;
  static float RadiationLifetime;
  //float acc[3];
  int ParticleClass;

  float AccretionRate[NTIMES];
  float AccretionRateTime[NTIMES];
  int TimeIndex;
  float oldmass; //To calculate accmass do accmass = mass - oldmass; oldmass = mass;

  static int FeedbackDistTotalCells, FeedbackDistRadius, FeedbackDistCellStep;
  float NotEjectedMass, eta_disk, mass_in_accretion_sphere, MassToBeEjected;
  float beta_jet, epsilon_deltat;
  float Accreted_angmom[MAX_DIMENSION];
  static std::vector<ParticleAttributeHandler *> AttributeHandlers;
};

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

int CalculateAccretedAngularMomentum();
template <class active_particle_class>
void ActiveParticleType_SmartStar::MergeSmartStars(
    int *nParticles, ActiveParticleList<ActiveParticleType>& ParticleList, 
    int *ngroups, LevelHierarchyEntry *LevelArray[], 
    int ThisLevel, ActiveParticleList<active_particle_class>& MergedParticles)
{
  int i,j;
  int dim;
  int GroupNumberAssignment[*nParticles];
  FLOAT* tempPos = NULL;
  int *groupsize = NULL;
  int **grouplist = NULL;

  HierarchyEntry** LevelGrids = NULL;

  int NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &LevelGrids);
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  /* Construct list of sink particle positions to pass to Foflist */
  FLOAT ParticleCoordinates[3*(*nParticles)];

  fflush(stdout);
  /* Particles merge once they come within 4 cells of one another */
  /* This is OK since we only want to merge particles when they come really
   * close to one another. Their accretion zones overlapping does 
   * not necessarily constitute a merger. 
   */
  FLOAT MergingRadius = LevelArray[ThisLevel]->GridData->GetCellWidth(0,0)*ACCRETIONRADIUS; 
  MergingRadius = MergingRadius*3.0;
  for (i=0; i<(*nParticles); i++) {
    tempPos = ParticleList[i]->ReturnPosition();
    for (dim=0; dim<3; dim++)
      ParticleCoordinates[3*i+dim] = tempPos[dim];
  }

  /* Find mergeable groups using an FOF search */
  *ngroups = FofList((*nParticles), ParticleCoordinates, MergingRadius, 
      GroupNumberAssignment, &groupsize, &grouplist);
  
  MergedParticles.reserve(*ngroups);
  
  /* Merge the mergeable groups */
 
  for (i=0; i<*ngroups; i++) {
    MergedParticles.copy_and_insert(
        *(static_cast<active_particle_class*>(ParticleList[grouplist[i][0]])));
    if (groupsize[i] != 1) {
      for (j=1; j<groupsize[i]; j++) {
    MergedParticles[i]->SmartMerge(
          static_cast<active_particle_class*>(ParticleList[grouplist[i][j]]));
        if (ParticleList[grouplist[i][j]]->DisableParticle(
                LevelArray, 
                MergedParticles[i]->ReturnCurrentGrid()->ReturnProcessorNumber()
              ) == FAIL)
        {
          ENZO_FAIL("MergeSmartStars: DisableParticle failed!\n");
        }
      }
    }
  }
  
  ParticleList.clear();

  free(groupsize);
  groupsize = NULL;
  for (i=0; i<*ngroups; i++)
    free(grouplist[i]);
  free(grouplist);
  grouplist = NULL;
  
  /* Loop over the grids and check if any of the merged particles have
     moved. If so, disable the particle on the current grid and assign
     it to the new grid*/

  int NewGrid = -1;

  for (i = 0; i < *ngroups; i++) {
    if (MergedParticles[i]->ReturnCurrentGrid()->PointInGrid(
            MergedParticles[i]->ReturnPosition()) == false) {
      // Find the grid to transfer to 
      for (j = 0; j < NumberOfGrids; j++) {
        if (LevelGrids[j]->GridData->PointInGrid(
                MergedParticles[i]->ReturnPosition())) {
          NewGrid = j;
          break;
        }
      }
      if (NewGrid == -1)
        ENZO_FAIL("Cannot assign particle to grid after merging!\n");
      int OldProc = MergedParticles[i]->CurrentGrid->ReturnProcessorNumber();
      active_particle_class *temp = static_cast<active_particle_class*>(
          MergedParticles[i]->clone());
      MergedParticles[i]->DisableParticle(LevelArray, 
          LevelGrids[NewGrid]->GridData->ReturnProcessorNumber()); 
      if (LevelGrids[NewGrid]->GridData->AddActiveParticle(
              static_cast<ActiveParticleType*>(temp)) == FAIL)
      	ENZO_FAIL("Active particle grid assignment failed!\n");
      if (MyProcessorNumber == OldProc) {
        MergedParticles.erase(i);
        MergedParticles.insert(*temp);
      }
      else if (MyProcessorNumber != temp->CurrentGrid->ReturnProcessorNumber())
        delete temp;
      MergedParticles[i]->AssignCurrentGrid(LevelGrids[NewGrid]->GridData);
    }
  }

  delete [] LevelGrids;

  *nParticles = *ngroups;
}


int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
 
int AssignActiveParticlesToGrids(
    ActiveParticleList<ActiveParticleType>& ParticleList, int nParticles, 
    LevelHierarchyEntry *LevelArray[]); 

template <class active_particle_class>
int ActiveParticleType_SmartStar::AfterEvolveLevel(
    HierarchyEntry *Grids[], TopGridData *MetaData,
    int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
    int ThisLevel, int TotalStarParticleCountPrevious[],
    int SmartStarID)
{

  /* SmartStar particles live on the maximum refinement level.  If we are on a lower level, this does not concern us */

  if (ThisLevel == MaximumRefinementLevel)
    {

      /* Generate a list of all sink particles in the simulation box */
      int i = 0, nParticles = 0, NumberOfMergedParticles = 0;
      ActiveParticleList<ActiveParticleType> ParticleList;
      FLOAT accradius = -10.0; //dummy
      
      ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID,
        ParticleList);

      /* Return if there are no smartstar particles */
      
      if (nParticles == 0)
        return SUCCESS;

      /* Calculate CellWidth on maximum refinement level */

      // This assumes a cubic box and may not work for simulations with MinimumMassForRefinementLevelExponent
      FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
        (MetaData->TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));

      /* Do Merging   */

      ActiveParticleList<active_particle_class> MergedParticles;
      
      /* Generate new merged list of sink particles */
      
      MergeSmartStars<active_particle_class>(
          &nParticles, ParticleList, 
          &NumberOfMergedParticles, LevelArray, ThisLevel, MergedParticles);

      ParticleList.clear();

      // Do merging twice to catch pathological cases where merging
      // leaves multiple sinks inside the same accretion zone.

      for (i = 0; i<NumberOfMergedParticles; i++)
      {
        ParticleList.copy_and_insert(
            *(static_cast<ActiveParticleType*>(MergedParticles[i])));
      }

      MergedParticles.clear();

      MergeSmartStars<active_particle_class>(
          &nParticles, ParticleList,
          &NumberOfMergedParticles, LevelArray, ThisLevel, MergedParticles);
      
      ParticleList.clear();
      
      if (debug)
        printf("Number of particles after merging: %"ISYM"\n",NumberOfMergedParticles);
      
      /* Assign local particles to grids */
      
      ParticleList.reserve(NumberOfMergedParticles);
      
      // need to use a bit of redirection because C++ pointer arrays have
      // trouble with polymorphism
      for (i = 0; i<NumberOfMergedParticles; i++)
      {
        ParticleList.copy_and_insert(
            *(static_cast<ActiveParticleType*>(MergedParticles[i])));
      }

      if (AssignActiveParticlesToGrids(ParticleList, NumberOfMergedParticles, 
              LevelArray) == FAIL)
        return FAIL;
      
      ParticleList.clear();
      MergedParticles.clear();
      
      /* Regenerate the global active particle list */
      
      ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID, 
        ParticleList);

      /* Do accretion */
     
      if (Accrete(nParticles, ParticleList, accradius, dx, LevelArray, 
              ThisLevel) == FAIL)
        ENZO_FAIL("SmartStar Particle accretion failed. \n");

      if(UpdateAccretionRateStats(nParticles, ParticleList, dx, LevelArray, ThisLevel) == FAIL)
	ENZO_FAIL("Failed to update accretion rate stats. \n");

      /* Apply feedback */
      if (SmartStarParticleFeedback(nParticles, ParticleList,
        dx, LevelArray, ThisLevel) == FAIL)
	ENZO_FAIL("SmartStar Particle Feedback failed. \n");
      /* This applies all of the updates made above */
      if (AssignActiveParticlesToGrids(ParticleList, NumberOfMergedParticles, 
              LevelArray) == FAIL)
        return FAIL;      

      ParticleList.clear();

    }

  return SUCCESS;
}

static FLOAT Dist(FLOAT *tempPos, FLOAT *tempPos1)
{
  return sqrt((tempPos[0] - tempPos1[0])*(tempPos[0] - tempPos1[0]) +
	      (tempPos[1] - tempPos1[1])*(tempPos[1] - tempPos1[1]) +
	      (tempPos[2] - tempPos1[2])*(tempPos[2] - tempPos1[2]));
}

/* Calculate the Innermost stable orbit (Abramowicz & Fragile (2013) Eqn 32) */
float R_ISCO(float a)
{
  float z1 = 1 + POW((1-a*a), 1.0/3.0)*(POW((1+a), 1.0/3.0) + POW((1-a), 1.0/3.0));
  float z2 = sqrt(3*a*a + z1*z1);
  float Z = 3 + z2 - sqrt((3-z1)*(3+z1+2*z2)); 
  return Z;
}
float MadauFit(float a, float mdot, float medddot);
