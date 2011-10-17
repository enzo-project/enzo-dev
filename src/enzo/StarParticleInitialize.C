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

int StarParticleInitialize(HierarchyEntry *Grids[], TopGridData *MetaData,
			   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			   int ThisLevel, Star *&AllStars,
			   int TotalStarParticleCountPrevious[])
{

  /* Return if this does not concern us */
  if (!(StarParticleCreation || StarParticleFeedback)) 
    return SUCCESS;

  LCAPERF_START("StarParticleInitialize");

  /* Set MetaData->NumberOfParticles and prepare TotalStarParticleCountPrevious
     these are to be used in CommunicationUpdateStarParticleCount 
     in StarParticleFinalize */  

  MetaData->NumberOfParticles = FindTotalNumberOfParticles(LevelArray);
  NumberOfOtherParticles = MetaData->NumberOfParticles - NumberOfStarParticles;
  RecordTotalStarParticleCount(Grids, NumberOfGrids, 
			       TotalStarParticleCountPrevious);

  /* Initialize the IMF lookup table if requested and not defined */

  if (PopIIIInitialMassFunction)
    StarParticlePopIII_IMFInitialize();

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

  for (cstar = AllStars; cstar; cstar = cstar->NextStar) {
    cstar->SetFeedbackFlag(TimeNow);
    cstar->CopyToGrid();
    cstar->MirrorToParticle();
  }

//  fprintf(stdout, "\nin StarParticleInitialize.C \n", MetaData->NumberOfParticles); 
//  fprintf(stdout, "MetaData->NumberOfParticles = %d\n", MetaData->NumberOfParticles); 
//  fprintf(stdout, "NumberOfStarParticles now = %d\n", NumberOfStarParticles);
//  fprintf(stdout, "NumberOfOtherParticles now = %d\n", NumberOfOtherParticles);


  LCAPERF_STOP("StarParticleInitialize");
  return SUCCESS;

}
