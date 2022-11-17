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
#define DEBUG 0 // SG. Changing back to debug 0.
/* Every how many times will the accretion rate be updated */
//#define FREQUENCY 100
#define MAXACCRETIONRADIUS  128 /* Times the minimum cell width */
#define ACCRETIONRADIUS  4 // SG. Change to something arbitrarily small for testing. Changing it back to 4 for formation purposes.
#define NUMRADIATIONBINS 5
#define CRITICAL_ACCRETION_RATE 0.001 //Msolar/yr (Haemerlee et al (2018))
#define TIMEGAP            100   // * timestep SG CHANGED FROM 100 TO 1 for testing
#define POPIII_RESOLUTION  0.001 //pc
#define SMS_RESOLUTION     0.1   //pc
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
#define POPII  3


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
    RadiationLifetime = 0.0;
    StellarAge = 0.0;
    ParticleClass = 0;
    NotEjectedMass = 0;
    MassToBeEjected = 0;
    eta_disk = 1 - sqrt(1 - 2.0/(3.0*R_ISCO(SmartStarSpin)));
    mass_in_accretion_sphere = 0;
    epsilon_deltat = 1.0;
    beta_jet = 0.0;
    InfluenceRadius = 0.0;
    for(int i = 0; i < 3; i++)
      {
	Accreted_angmom[i] = 0.0;
      }
    TimeIndex = -1;
    //oldmass = -1;
    for(int i = 0; i < NTIMES; i++)
      {
	AccretionRateTime[i] = -222222;
	AccretionRate[i] = -222222;
      }
  };
  ActiveParticleType_SmartStar(ActiveParticleType_SmartStar* part) :
    ActiveParticleType(static_cast<ActiveParticleType*>(part)) {   
    AccretionRadius = part->AccretionRadius;
    RadiationLifetime = part->RadiationLifetime;
    StellarAge = part->StellarAge;
    ParticleClass = part->ParticleClass;
    for(int i = 0; i < 3; i++) {
      Accreted_angmom[i] = 0.0;
    }
    TimeIndex = part->TimeIndex;
    
    //oldmass = part->oldmass;
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
    InfluenceRadius = part->InfluenceRadius;
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
  template <class active_particle_class> active_particle_class *copy(void);
    void  IncreaseLevel() { level++; AccretionRadius *= 2; fprintf(stderr, "%s: (void) AccretionRadius increased.\n", __FUNCTION__);}; // SG. Update AccretionRadius with each increase in level.
    void  ReduceLevel() { level--; AccretionRadius /= 2; fprintf(stderr, "%s: (void) AccretionRadius reduced.\n"), __FUNCTION__;}; // SG. Update AccretionRadius with each increase in level.
    void  ReduceLevel(int x) { level -= x; AccretionRadius /= 2; fprintf(stderr, "%s: (int x) AccretionRadius reduced.\n", __FUNCTION__);};
    void  IncreaseLevel(int x) { level += x; AccretionRadius *= 2; fprintf(stderr, "%s: (int x) AccretionRadius increased.\n", __FUNCTION__);};

  
  
  static int SmartStarParticleFeedback(
             int nParticles, 
             ActiveParticleList<ActiveParticleType>& ParticleList,
		     FLOAT dx, LevelHierarchyEntry *LevelArray[], int ThisLevel, int SmartStarID);
  

  static int ResetAcceleration(float *ActiveParticleAcceleration);
  static int CreateParticle(grid *thisgrid_orig, ActiveParticleFormationData &supp_data,
			    int particle_index);
  static int InitializeParticleType();

  void SmartMerge(ActiveParticleType_SmartStar *a);
  void AssignMassFromIMF();
  int CalculateAccretedAngularMomentum();
  int SmartStarAddFeedbackSphere();
  int DetermineSEDParameters(FLOAT Time, FLOAT dx);
  ENABLED_PARTICLE_ID_ACCESSOR
  bool IsARadiationSource(FLOAT Time);
  // SG. New Function.
  int FindAccretionSphere(LevelHierarchyEntry *LevelArray[], int level, 
           float StarLevelCellWidth, float &Radius, float TargetSphereMass, float &MassEnclosed, 
           float &Metallicity2, float &Metallicity3, float &ColdGasMass, float &ColdGasFraction,
			     int &SphereContained, bool &MarkedSubgrids);
  
  // sink helper routines

  template <class active_particle_class>
  static void MergeSmartStars(
      int *nParticles, ActiveParticleList<ActiveParticleType>& ParticleList, 
      int *ngroups, 
      LevelHierarchyEntry *LevelArray[], int ThisLevel, 
      ActiveParticleList<active_particle_class>& MergedParticles);
      
// SG. Took dx out of argument declaration in the following 3 functions + added SmartStarID to accrete.
  static int Accrete(int nParticles, 
      ActiveParticleList<ActiveParticleType>& ParticleList,
      FLOAT &AccretionRadius,
      LevelHierarchyEntry *LevelArray[], int ThisLevel, int SmartStarID);

  static int UpdateAccretionRateStats(int nParticles,
				      ActiveParticleList<ActiveParticleType>& ParticleList,
				      LevelHierarchyEntry *LevelArray[], int ThisLevel);

  static int UpdateRadiationLifetimes(int nParticles,
				      ActiveParticleList<ActiveParticleType>& ParticleList,
				      LevelHierarchyEntry *LevelArray[], int ThisLevel);

  static int  RemoveMassFromGridAfterFormation(int nParticles, 
              ActiveParticleList<ActiveParticleType>& ParticleList,
              LevelHierarchyEntry *LevelArray[], int ThisLevel);

  // SG. New Func within RemoveMassFromGridAfterFormation.
  static int  PopIIIFormationFromSphere(ActiveParticleType_SmartStar* SS, 
              grid* APGrid, int ThisProcessorNum, FLOAT StarLevelCellWidth, 
              FLOAT CellVolumeStarLevel,FLOAT Time, LevelHierarchyEntry *LevelArray[], 
              LevelHierarchyEntry *Temp, int ThisLevel);

  static float EjectedMassThreshold;
  FLOAT AccretionRadius;   // in units of CellWidth on the maximum refinement level
 
  static double LuminosityPerSolarMass;
  static int RadiationSEDNumberOfBins;
  static float* RadiationEnergyBins;
  static float* RadiationSED;
  double RadiationLifetime, StellarAge;
  //float acc[3];
  int ParticleClass;
  //float oldmass; //To calculate accmass do accmass = mass - oldmass; oldmass = mass;
  float AccretionRate[NTIMES];
  float AccretionRateTime[NTIMES];
  int TimeIndex;
  static int FeedbackDistTotalCells, FeedbackDistRadius, FeedbackDistCellStep;
  float NotEjectedMass, eta_disk, mass_in_accretion_sphere, MassToBeEjected;
  float beta_jet, epsilon_deltat;
  float Accreted_angmom[MAX_DIMENSION];
  float InfluenceRadius;
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

  // SG. Trying to get rid of merging error for POPIII particles.
  // for (i=0; i<(*nParticles); i++) {
  //   if(ParticleList[i]->ReturnType() == POPIII){
  //     continue;
  //   }
  // }

  int NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &LevelGrids);
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  /* Construct list of sink particle positions to pass to Foflist */
  FLOAT ParticleCoordinates[3*(*nParticles)];
  /* Particles merge once they come within 1 accretion radii of one another */
  //fprintf(stderr, "%s: no. particles = %"ISYM"\n.", __FUNCTION__, *nParticles);
  FLOAT MergingRadius = 1.5*(LevelArray[ThisLevel]->GridData->GetCellWidth(0,0))*ACCRETIONRADIUS; // SG. Making merging radius = 0. Unmaking it 0. How does merging radius relate to feedback radius?
  //fprintf(stderr, "%s: merging radius = %e\n.", __FUNCTION__, MergingRadius);

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
          ENZO_FAIL("MergeSmartStars: DisableAParticle failed!\n");
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
              // // SG. Debugging.
              // if(ParticleList[i]->ReturnType() == POPIII){
              //   continue;
              //  }
      // Find the grid to transfer to
      for (j = 0; j < NumberOfGrids; j++) {
        if (LevelGrids[j]->GridData->PointInGrid(
                MergedParticles[i]->ReturnPosition())) {
          NewGrid = j;
          fprintf(stderr,"%s: new grid = %"ISYM".\n", __FUNCTION__, NewGrid);
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
// SG. Working with all grids on ThisLevel. Maybe try to grab level the SS is on?


      /* SmartStar particles live on the maximum refinement level.  If we are on a lower level, this does not concern us */
      int i,j,k,index;

      /* Generate a list of all sink particles in the simulation box */
      int nParticles = 0, NumberOfMergedParticles = 0;
      ActiveParticleList<ActiveParticleType> ParticleList;
      FLOAT accradius = -10.0; //dummy
      // SG. Units for dx to pc conversion.
      float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
      FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
      double MassUnits;
      GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);

      //   // SG. Return if ThisLevel != APGrid level -> doesn't work as particle traverses levels.
      // for (int i = 0; i < nParticles; i++) {
      //   grid* APGrid = ParticleList[i]->ReturnCurrentGrid();
			// 	int MyLevel = APGrid->GridLevel;
			// 	if (ThisLevel != MyLevel){
      //     return SUCCESS;
      //   }
      // }

      ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID,
        ParticleList);

      /* Return if there are no smartstar particles */

      if (nParticles == 0){
        return SUCCESS;
      } // end if

      LevelHierarchyEntry *Temp = NULL;
	    HierarchyEntry *Temp2 = NULL;
      Temp = LevelArray[ThisLevel];

      while (Temp != NULL) {

        // // SG. Check we're on the maximum LOCAL refinement level from the get-go.
        // for (k = Temp->GridStartIndex[2]; k <= Temp->GridEndIndex[2]; k++) {
        //   for (j = Temp->GridStartIndex[1]; j <= Temp->GridEndIndex[1]; j++) {
        //   index = GRIDINDEX_NOGHOST(Temp->GridStartIndex[0], j, k);
        //     for (i = Temp->GridStartIndex[0]; i <= Temp->GridEndIndex[0]; i++, index++) {
        // if (Temp->BaryonField[Temp->NumberOfBaryonFields][index] != 0.0){
        //   return SUCCESS;
        // } else{
        
        /* Zero under subgrid field */
      
          Temp->GridData->
      ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
          Temp2 = Temp->GridHierarchyEntry->NextGridNextLevel;
          while (Temp2 != NULL) { // SG. this is doing the check 1 or 0 in baryon refinement field
      Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
                ZERO_UNDER_SUBGRID_FIELD);
      Temp2 = Temp2->NextGridThisLevel;
          }
        
        Temp = Temp->NextGridThisLevel; // how we loop over all grids on the level.

      } // END: Grids


      ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID,
        ParticleList);
    

      /* Calculate CellWidth on maximum refinement level */
      // SG. May need to fix this.
      FLOAT dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
        (MetaData->TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(14))); // SG. Replaced MaximumRefinementLevel with ThisLevel.
      // //fprintf(stderr,"%s: CellWidth dx = %e and ThisLevel = %"ISYM" (but dx calculated for level 14).\n", __FUNCTION__, dx*LengthUnits/pc_cm, ThisLevel); fflush(stdout);

      /* Remove mass from grid from newly formed particles */
      RemoveMassFromGridAfterFormation(nParticles, ParticleList, 
				       LevelArray, ThisLevel);
 
      /* Clean any particles marked for deletion */
      for (i = 0; i<nParticles; i++) {
	if(ParticleList[i]->ShouldDelete() == true) {
	  printf("%s: Delete SS %d following RemoveMassFromGridAfterFormation\n", __FUNCTION__,
		 static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->ReturnID());
	  fflush(stdout);
	  ParticleList.erase(i);
	  i = -1;
	  nParticles--;
	}
      }
      ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID,
      			    ParticleList);
      if (AssignActiveParticlesToGrids(ParticleList, nParticles, 
      				       LevelArray) == FAIL)
	return FAIL;
      if (debug)
	fprintf(stderr,"%s: Number of particles before merging: %"ISYM"\n",__FUNCTION__, nParticles);fflush(stdout);
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
        fprintf(stderr,"%s: Number of particles after merging: %"ISYM"\n",__FUNCTION__, NumberOfMergedParticles);fflush(stdout);
      
      /* Assign local particles to grids */
      /* Do Merging   */
     
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

      if (Accrete(nParticles, ParticleList, accradius, LevelArray, // SG. Error is being thrown here. AFTER AssignActiveParticlesToGrid was called before it.
              ThisLevel, SmartStarID) == FAIL)
        ENZO_FAIL("SmartStar Particle accretion failed. \n");

      if(UpdateAccretionRateStats(nParticles, ParticleList, LevelArray, ThisLevel) == FAIL)
	      ENZO_FAIL("Failed to update accretion rate stats. \n");

    //  // SG. For debugging.
    //   fprintf(stderr, "%s: edge 0 = %e, edge 1 = %e and edge 3 = %e.\n", __FUNCTION__, APGrid->GetGridLeftEdge(0), APGrid->GetGridLeftEdge(1), APGrid->GetGridLeftEdge(2));
		// 	fprintf(stderr, "%s: edge 0 = %e, edge 1 = %e and edge 3 = %e.\n", __FUNCTION__, APGrid->GetGridRightEdge(0), APGrid->GetGridRightEdge(1), APGrid->GetGridRightEdge(2));
      if(UpdateRadiationLifetimes(nParticles, ParticleList, LevelArray, ThisLevel) == FAIL)
	ENZO_FAIL("Failed to update radiation lifetimes. \n");

    //  // SG. For debugging.
    //   fprintf(stderr, "%s: edge 0 = %e, edge 1 = %e and edge 3 = %e.\n", __FUNCTION__, APGrid->GetGridLeftEdge(0), APGrid->GetGridLeftEdge(1), APGrid->GetGridLeftEdge(2));
		// 	fprintf(stderr, "%s: edge 0 = %e, edge 1 = %e and edge 3 = %e.\n", __FUNCTION__, APGrid->GetGridRightEdge(0), APGrid->GetGridRightEdge(1), APGrid->GetGridRightEdge(2));
    
      /* Apply feedback */
      if (SmartStarParticleFeedback(nParticles, ParticleList,
        dx, LevelArray, ThisLevel, SmartStarID) == FAIL)
	ENZO_FAIL("SmartStar Particle Feedback failed. \n");

      
      /* Clean any particles marked for deletion. 
       * After each deletion I need to reloop and check it again. 
       */
      for (i = 0; i<nParticles; i++) {
	if(ParticleList[i]->ShouldDelete() == true) { 
	  printf("%s: Delete SS %d following Feedback\n", __FUNCTION__,
		 static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->ReturnID());
	  fflush(stdout);
	  ParticleList.erase(i);
	  i = -1; //reset counter
	  nParticles--;
	}
      }
      ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID, 
       ParticleList);
      /* This applies all of the updates made above */
      if (AssignActiveParticlesToGrids(ParticleList, nParticles, 
              LevelArray) == FAIL)
        return FAIL;      
      ParticleList.clear();

      //       }
      //     }
      //   }
      // }
  return SUCCESS;
} // End AfterEvolveLevel function

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
