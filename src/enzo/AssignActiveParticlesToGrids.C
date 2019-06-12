/***********************************************************************
/
/  AssignActiveParticlesToGrids:
/  Assigning a list of active particles to the appropriate grid.
/
/  written by: Nathan Goldbaum
/  date:       June, 2012
/  modified1:
/
/  PURPOSE: Given a list of active particles, this function will 
/           traverse the hierarchy and assign them to the appropriate 
/           grid.  If a particle with the same ID as the one under 
/           consideration is already assigned to the wrong grid, that 
/           particle is deleted.  If the particle is already assigned 
/           to the correct grid, the mirrored particle is updated.
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
#include "phys_constants.h"

#define NO_DEBUG_ASSIGN

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],int NumberOfGrids);

int AssignActiveParticlesToGrids(
    ActiveParticleList<ActiveParticleType>& ParticleList, 
    int nParticles, 
    LevelHierarchyEntry *LevelArray[]) 
{
  int LevelMax, GlobalLevelMax = -1, SavedGrid, NumberOfGrids = -1, i, level, NumberOfLevelGrids = -1, gridnum;
  int dim, SavedGridOffProc = -1, ProcessorNumber = -1;
  HierarchyEntry **LevelGrids = NULL;
  FLOAT* pos, pos1[3];
  float mass;
  
  FLOAT period[3];
  for (dim = 0; dim < 3; dim++) {
    period[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
  }
  
  for (i = 0; i<nParticles; i++) {
    pos = ParticleList[i]->ReturnPosition();
    for (dim = 0; dim < 3; dim++) {
      pos1[dim] = fmod(pos[dim], period[dim]);
      if (pos1[dim] < 0) {
        pos1[dim] += period[dim];
      }
    }
    // Find the grid and processor this particle lives on
    GlobalLevelMax = LevelMax = SavedGrid = SavedGridOffProc = -1;
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      NumberOfLevelGrids = GenerateGridArray(LevelArray, level, &LevelGrids);     
      for (gridnum = 0; gridnum < NumberOfLevelGrids; gridnum++) {
        if (LevelGrids[gridnum]->GridData->PointInGrid(pos1) == true) {
	      SavedGridOffProc = gridnum;
	      if (LevelGrids[gridnum]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
	        if (LevelGrids[gridnum]->GridData->isLocal() == true) { 
	          SavedGrid = gridnum;
	          LevelMax = level;
	        }
	      }
          GlobalLevelMax = level;
	    } // if MyProc
	  } // for gridnum
      delete [] LevelGrids;
      LevelGrids = NULL;
    }
    
    /* Assign the merged particles to grids. */
    
    if (NumberOfProcessors == 1) {
      grid* OldGrid = ParticleList[i]->ReturnCurrentGrid();
      int ID = ParticleList[i]->ReturnID();
      int foundAP = FALSE;
      NumberOfGrids = GenerateGridArray(LevelArray, LevelMax, &LevelGrids); 

      if (SavedGrid != -1) {
        foundAP = OldGrid->RemoveActiveParticle(ID,LevelGrids[SavedGrid]->
            GridData->ReturnProcessorNumber());
      }
      if (LevelGrids[SavedGrid]->GridData->AddActiveParticle(
              static_cast<ActiveParticleType*>(ParticleList[i])) == FAIL)
        ENZO_FAIL("Active particle grid assignment failed!\n");
    }
    else {
#ifdef USE_MPI
      /* Find the processor which has the maximum value of
	  LevelMax and assign the particle to the
	  SavedGrid on that processor.  */
      struct { Eint32 value; Eint32 rank; } sendbuf, recvbuf;
      MPI_Comm_rank(MPI_COMM_WORLD, &sendbuf.rank);
      NumberOfGrids = GenerateGridArray(LevelArray, GlobalLevelMax, &LevelGrids); 
      // We're moving it, make sure that the particle position is fixed (if required).
      ParticleList[i]->SetPositionPeriod(period);
      // Below this effectively removes the particle from the sending proc.
      grid* OldGrid = ParticleList[i]->ReturnCurrentGrid();
      int ID = ParticleList[i]->ReturnID();
      OldGrid->RemoveActiveParticle(ID,LevelGrids[SavedGridOffProc]->GridData->ReturnProcessorNumber());
      // if this is the receiving proc, add it.
      if (LevelMax == GlobalLevelMax) {
	    if (LevelGrids[SavedGrid]->GridData->AddActiveParticle(
                static_cast<ActiveParticleType*>(ParticleList[i])) == FAIL) {
	      ENZO_FAIL("Active particle grid assignment failed"); 
	    } 
      }
      LevelMax = GlobalLevelMax;
#endif // endif parallel
    } // end else
    
    /* Sync the updated particle counts accross all proccessors */
    
    CommunicationSyncNumberOfParticles(LevelGrids, NumberOfGrids);
    
    delete [] LevelGrids;
    
  }

  return SUCCESS;
}
