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

int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids);

int ParticleSplitter(HierarchyEntry *Grids[], LevelHierarchyEntry *LevelArray[], int ThisLevel,
		     TopGridData *MetaData, int NumberOfGrids)
{

  /* Return if this does not concern us */
  if (ParticleSplitterIterations == 0 || !(MetaData->FirstTimestepAfterRestart)) 
    return SUCCESS;

  int level, i;
  LevelHierarchyEntry *Temp;
  grid *GridPointer[MAX_NUMBER_OF_SUBGRIDS];
  FLOAT TimeNow = LevelArray[ThisLevel]->GridData->ReturnTime();

  /* Initialize all star particles if this is a restart */

  for (i = 0; i < ParticleSplitterIterations; i++)
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
      fprintf(stdout, "ParticleSplitter [level=%d] starts. \n", level);
      for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
	fprintf(stdout, "ParticleSplitter [grid->NumberOfParticles=%d] starts. \n", 
		Temp->GridData->ReturnNumberOfParticles());
	if (Temp->GridData->ParticleSplitter(level) == FAIL) {
	  fprintf(stderr, "Error in grid::ParticleSplitter.\n");
	  ENZO_FAIL("");
	}
      }
    }

  /* Now redistribute the particles as the newly created particles 
     might have crossed the grid boundaries */
  
  RebuildHierarchy(MetaData, LevelArray, 0);

  /* Update the star particle counters using the same routine in StarParticleFinalize */

  if (CommunicationUpdateStarParticleCount(Grids, MetaData,
					   NumberOfGrids) == FAIL) {
    fprintf(stderr, "Error in CommunicationUpdateStarParticleCount.\n");
    ENZO_FAIL("");
  }

  return SUCCESS;

}
