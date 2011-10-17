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

#ifdef USE_MPI
#include "mpi.h"
#endif

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
#include "CosmologyParameters.h"
#include "CommunicationUtilities.h"

#define NO_DEBUG_PS 

int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids,
					 int TotalStarParticleCountPrevious[]);
int FindTotalNumberOfParticles(LevelHierarchyEntry *LevelArray[]);
void RecordTotalStarParticleCount(HierarchyEntry *Grids[], int NumberOfGrids,
				  int TotalStarParticleCountPrevious[]);

int ParticleSplitter(LevelHierarchyEntry *LevelArray[], int ThisLevel,
		     TopGridData *MetaData)
{

  int level, i, grid1;
  HierarchyEntry **Grids;
  int NumberOfGrids;

  /* First, rebuild the hierarchy */
  
  RebuildHierarchy(MetaData, LevelArray, 0);  

  /* Return if this does not concern us */

  if (ParticleSplitterIterations <= 0 || 
      ParticleSplitterChildrenParticleSeparation <=0) 
    return SUCCESS;

  /* Find total NumberOfParticles in all grids; this is needed in 
     CommunicationUpdateStarParticleCount below */

  MetaData->NumberOfParticles = FindTotalNumberOfParticles(LevelArray);
  NumberOfOtherParticles = MetaData->NumberOfParticles - NumberOfStarParticles;

  /* Initialize all star particles if this is a restart */

  if (ParticleSplitterIterations > 1)
    fprintf(stderr, "WARNING: ParticleSplitterIterations > 1 is not properly tested yet.\n");

  for (i = 0; i < ParticleSplitterIterations; i++) {

    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {

#ifdef DEBUG_PS
      fprintf(stdout, "MetaData->NumberOfParticles when ParticleSplitter [level=%d] starts = %d\n", 
	      level, MetaData->NumberOfParticles);
#endif

      NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
      int *TotalStarParticleCountPrevious = new int[NumberOfGrids];

      /* Record the star particle number on each grid,
	 will be used in CommunicationUpdateStarParticleCount */

      RecordTotalStarParticleCount(Grids, NumberOfGrids, 
				   TotalStarParticleCountPrevious);

      //      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) 
      //	fprintf(stdout, "TotalStarParticleCountPrevious[grid=%d] = %d\n", grid1, 
      //		TotalStarParticleCountPrevious[grid1]);
      
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

#ifdef DEBUG_PS
	fprintf(stdout, "ParticleSplitter [grid->NumberOfParticles=%d] starts. \n", 
		Grids[grid1]->GridData->ReturnNumberOfParticles());
#endif
	
	if (Grids[grid1]->GridData->ParticleSplitter(level) == FAIL) {
	  ENZO_FAIL("Error in grid::ParticleSplitter.\n");
	}

      }  // loop for grid1

      CommunicationBarrier(); 

      /* Assign indices for star particles and update the star particle counters 
	 using the same routine in StarParticleFinalize */
      
      if (CommunicationUpdateStarParticleCount(Grids, MetaData,
					       NumberOfGrids, 
					       TotalStarParticleCountPrevious) == FAIL) {
	ENZO_FAIL("Error in CommunicationUpdateStarParticleCount.\n");
      }

    }  // loop for level

    /* Redistribute the particles as the newly created particles 
       might have crossed the grid boundaries */
  
    RebuildHierarchy(MetaData, LevelArray, 0);  
    
    /* Set MetaData->NumberOfParticles again; might be needed somewhere */

    MetaData->NumberOfParticles = FindTotalNumberOfParticles(LevelArray);

#ifdef DEBUG_PS
    fprintf(stdout, "ParticleSplitter 1 cycle done!\n");
    fprintf(stdout, "MetaData->NumberOfParticles = %d\n", MetaData->NumberOfParticles);
    fprintf(stdout, "NumberOfStarParticles now = %d\n", NumberOfStarParticles);
#endif

  }  // loop for i

  /* Set splitter parameter zero; otherwise it will split particles again at next restart */

  ParticleSplitterIterations = 0;

  return SUCCESS;

}










/* Find the total number of particles in all grids
   to update MetaData->NumberOfParticles */

int FindTotalNumberOfParticles(LevelHierarchyEntry *LevelArray[])
{

  int level, TotalNumberOfParticles = 0;
  LevelHierarchyEntry *Temp;

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) 
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      TotalNumberOfParticles += Temp->GridData->ReturnNumberOfParticles();

  return TotalNumberOfParticles;

}









/* Record the star particle number on each grid; 
   mostly used in CommunicationUpdateStarParticleCount */

void RecordTotalStarParticleCount(HierarchyEntry *Grids[], int NumberOfGrids,
				  int TotalStarParticleCountPrevious[])
{

  int grid, *PartialStarParticleCountPrevious = new int[NumberOfGrids];

  for (grid = 0; grid < NumberOfGrids; grid++) {
    TotalStarParticleCountPrevious[grid] = 0;
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      PartialStarParticleCountPrevious[grid] =
	Grids[grid]->GridData->ReturnNumberOfStarParticles();
    }
    else {
      PartialStarParticleCountPrevious[grid] = 0;
    }
  }
 
#ifdef USE_MPI
  /* Get counts from each processor to get total list of new particles. */
 
  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;

  MPI_Arg GridCount = NumberOfGrids;
   
  MPI_Allreduce(PartialStarParticleCountPrevious, TotalStarParticleCountPrevious, GridCount,
		DataTypeInt, MPI_SUM, MPI_COMM_WORLD);
#endif

  delete [] PartialStarParticleCountPrevious;

}
