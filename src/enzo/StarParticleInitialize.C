/***********************************************************************
/
/  STAR PARTICLE INITIALIZATION
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE: Contains all routines to initialize the star particles.
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "performance.h"
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

int StarParticlePopIII_IMFInitialize(void);
int StarParticleFindAll(LevelHierarchyEntry *LevelArray[], Star *&AllStars);
int StarParticleMergeNew(LevelHierarchyEntry *LevelArray[], Star *&AllStars);
int StarParticleMergeMBH(LevelHierarchyEntry *LevelArray[], Star *&AllStars);
int FindTotalNumberOfParticles(LevelHierarchyEntry *LevelArray[]);
void RecordTotalStarParticleCount(HierarchyEntry *Grids[], int NumberOfGrids,
				  int TotalStarParticleCountPrevious[]);


int StarParticleIndividual_IMFInitialize(void);
int IndividualStarProperties_Initialize(TopGridData &MetaData);
int IndividualStarRadiationProperties_Initialize(void);
int InitializeStellarYields(const float & time);

int StarParticleInitialize(HierarchyEntry *Grids[], TopGridData *MetaData,
			   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			   int ThisLevel, Star *&AllStars,
			   int TotalStarParticleCountPrevious[]
#ifdef INDIVIDUALSTAR
                           , int SkipFeedbackFlag = 0
#endif
                           )


{

  /* Return if this does not concern us */
  if (!(StarParticleCreation || StarParticleFeedback)) 
    return SUCCESS;

  LCAPERF_START("StarParticleInitialize");
  TIMER_START("StarParticleInitialize");

  /* Set MetaData->NumberOfParticles and prepare TotalStarParticleCountPrevious
     these are to be used in CommunicationUpdateStarParticleCount 
     in StarParticleFinalize */  

  MetaData->NumberOfParticles = FindTotalNumberOfParticles(LevelArray);
  NumberOfOtherParticles = MetaData->NumberOfParticles - NumberOfStarParticles;
  RecordTotalStarParticleCount(Grids, NumberOfGrids, 
			       TotalStarParticleCountPrevious);

  /* Initialize the IMF lookup table if requested and not defined */

  if (PopIIIInitialMassFunction && STARMAKE_METHOD(INDIVIDUAL_STAR) == FALSE)
    StarParticlePopIII_IMFInitialize();

  /* Initialize IMF lookup table if needed and radiation table if needed */
  if(STARMAKE_METHOD(INDIVIDUAL_STAR)){
    StarParticleIndividual_IMFInitialize();

    /* Initialize individual star properties (L, T, R) */
    IndividualStarProperties_Initialize(*MetaData);

    /* Initialize radiation data table */
    if((RadiativeTransfer && IndividualStarBlackBodyOnly == FALSE) || IndividualStarFUVHeating){
      IndividualStarRadiationProperties_Initialize();
    }

    /* StellarYields */
    InitializeStellarYields(MetaData->Time);

  }

  int level, grids;
  Star *cstar;
  LevelHierarchyEntry *Temp;
  grid *GridPointer[MAX_NUMBER_OF_SUBGRIDS];
  FLOAT TimeNow = LevelArray[ThisLevel]->GridData->ReturnTime();

  /* Initialize all star particles if this is a restart */

  if (MetaData->FirstTimestepAfterRestart)
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++)
      for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
	if (Temp->GridData->FindAllStarParticles(level) == FAIL) {
	  	  ENZO_FAIL("Error in grid::FindAllStarParticles.");
	}

  /* Create a master list of all star particles */

  if (StarParticleFindAll(LevelArray, AllStars) == FAIL) {
        ENZO_FAIL("Error in StarParticleFindAll.");
  }

  if (MetaData->FirstTimestepAfterRestart == FALSE) {

    /* Merge any newly created, clustered particles */

    /* If individual stars, no merging */
    if (StarParticleMergeNew(LevelArray, AllStars) == FAIL) {
        ENZO_FAIL("Error in StarParticleMergeNew.");
    }

  /* Merge MBH particles that are close enough.  Ji-hoon Kim, Sep.2009 */

    if (StarParticleMergeMBH(LevelArray, AllStars) == FAIL) {
      ENZO_FAIL("Error in StarParticleMergeMBH.\n");
    }

  } // ENDIF !restart

  /* 
     Set feedback flags.  

     Sync all star and normal particles that are stored in the grids
     to the global list (AllStars) so these changes are reflected
     there.
  */

//  if (MyProcessorNumber == ROOT_PROCESSOR) {

//    for (cstar = AllStars; cstar; cstar = cstar->NextStar)
//      cstar->PrintInfo();
//  }

#ifdef INDIVIDUALSTAR
  if (SkipFeedbackFlag){
#endif
  AllStars->MakeStarsMap();
  for (cstar = AllStars; cstar; cstar = cstar->NextStar) {
    float dtForThisStar   = LevelArray[ThisLevel]->GridData->ReturnTimeStep();

    cstar->SetFeedbackFlag(TimeNow, dtForThisStar);
    cstar->CopyToGrid();
    cstar->MirrorToParticle();

    if(STARMAKE_METHOD(INDIVIDUAL_STAR))
      cstar->AssertInterpolationPositions(); // should be set at init, but double check
//      cstar->AssignInterpolationTablePositions();


  }
#ifdef INDIVIDUALSTAR
  }
#endif


//  fprintf(stdout, "\nin StarParticleInitialize.C \n", MetaData->NumberOfParticles); 
//  fprintf(stdout, "MetaData->NumberOfParticles = %d\n", MetaData->NumberOfParticles); 
//  fprintf(stdout, "NumberOfStarParticles now = %d\n", NumberOfStarParticles);
//  fprintf(stdout, "NumberOfOtherParticles now = %d\n", NumberOfOtherParticles);


  LCAPERF_STOP("StarParticleInitialize");
  TIMER_STOP("StarParticleInitialize");
  return SUCCESS;

}
