/***********************************************************************
/
/  PARTICLE SPLITTER
/
/  written by: Ji-hoon Kim
/  date:       October, 2009
/  modified1:
/
/  PURPOSE: Contains routines to split particles.  Ideal for setting up 
/           an initial condition with a relatively low computational cost,
/           and then restarting for an extremely high-resolution 
/           re-simulation. See Grid_ParticleSplitter.C for details.
/
/           Currently it implicitly assumes that only DM and conventional 
/           star particles get split.  Other particles - which usually 
/           become Star class particles - seem to have no reason to be 
/           split.  (as of Oct.2009)
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
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
#include "StarParticleData.h"

#define DEBUG_PS

int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids);

int ParticleSplitter(LevelHierarchyEntry *LevelArray[], int ThisLevel,
		     TopGridData *MetaData)
{

  /* Return if this does not concern us */

  if (ParticleSplitterIterations <= 0 || !(MetaData->FirstTimestepAfterRestart)) 
    return SUCCESS;

  int level, i, grid1;
  HierarchyEntry **Grids;
  int NumberOfGrids;

  /* Set MetaData->NumberOfParticles; this is needed in 
     CommunicationUpdateStarParticleCount below */

  MetaData->NumberOfParticles = 0;

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
      NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) 
	MetaData->NumberOfParticles += Grids[grid1]->GridData->ReturnNumberOfParticles();
  }
	
#ifdef DEBUG_PS
  fprintf(stdout, "MetaData->NumberOfParticles = %d\n", MetaData->NumberOfParticles);
#endif

  /* Initialize all star particles if this is a restart */

  for (i = 0; i < ParticleSplitterIterations; i++) {

    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {

#ifdef DEBUG_PS
      fprintf(stdout, "ParticleSplitter [level=%d] starts. \n", level);
#endif
      NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

#ifdef DEBUG_PS
	fprintf(stdout, "ParticleSplitter [grid->NumberOfParticles=%d] starts. \n", 
		Grids[grid1]->GridData->ReturnNumberOfParticles());
#endif
	
	if (Grids[grid1]->GridData->ParticleSplitter(level) == FAIL) {
	  fprintf(stderr, "Error in grid::ParticleSplitter.\n");
	  ENZO_FAIL("");
	}

      }  // loop for grid1

      /* Update the star particle counters using the same routine 
      in StarParticleFinalize */
      
      if (CommunicationUpdateStarParticleCount(Grids, MetaData,
					       NumberOfGrids) == FAIL) {
	fprintf(stderr, "Error in CommunicationUpdateStarParticleCount.\n");
	ENZO_FAIL("");
      }

    }  // loop for level

    /* Redistribute the particles as the newly created particles 
       might have crossed the grid boundaries */
  
    RebuildHierarchy(MetaData, LevelArray, 0);
    
    
#ifdef DEBUG_PS
    fprintf(stdout, "NumberOfStarParticles now = %d\n", NumberOfStarParticles);
#endif

  }  // loop for i

  return SUCCESS;

}
