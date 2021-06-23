/***********************************************************************
/
/  COMMUNICATION ROUTINE: UPDATE STAR PARTICLE COUNT
/
/  written by: Greg Bryan
/  date:       June, 1999
/  modified1:  Robert Harkness
/  date:       March, 2005
/
/  modified2:  Ji-hoon Kim
/  date:       November, 2009
/
/              Now this new routines can correctly assigns indices to new 
/              particles (whether they are star or DM) and keep track 
/              of NumberOfStarParticles and NumberOfOtherParticles.
/              This function has to have the arguments 
/              Total(Star)ParticleCountPrevious, which can be achieved
/              by RecordTotalParticleCount. (check StarParticleInitialize)
/  modified3:  John Wise
/  date:       December, 2011 -- active particles
/
/  PURPOSE:
/    This routines keeps track of the number of new active particles.
/
************************************************************************/
 
#ifdef USE_MPI
#endif /* USE_MPI */
 
#include "preincludes.h"
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"

#define NO_DEBUG

int CommunicationUpdateActiveParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids, 
					 int NumberOfNewActiveParticles[])
{

  LCAPERF_START("UpdateActiveParticleCount");

  int grid, *TotalNewActiveParticleCount = new int[NumberOfGrids],
      *PartialNewActiveParticleCount = new int[NumberOfGrids];
 
  /* Set ParticleCount to zero and record number of particles for grids
     on this processor. */
 
  for (grid = 0; grid < NumberOfGrids; grid++) {
    TotalNewActiveParticleCount[grid] = 0;
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      PartialNewActiveParticleCount[grid] = NumberOfNewActiveParticles[grid];	
    }
    else {
      PartialNewActiveParticleCount[grid] = 0;
    }
}
 
#ifdef USE_MPI
 
  /* Get counts from each processor to get total list of new particles. */
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;

  MPI_Arg GridCount = NumberOfGrids;
   
  MPI_Allreduce(PartialNewActiveParticleCount, TotalNewActiveParticleCount, GridCount,
		DataTypeInt, MPI_SUM, MPI_COMM_WORLD);

#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[11] += endtime-starttime;
  counter[11] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
#endif /* USE_MPI */
 
  /* Set new particle count for each grid. */
 
  for (grid = 0; grid < NumberOfGrids; grid++) {

    NumberOfActiveParticles += TotalNewActiveParticleCount[grid];
    
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) 
      /* If this grid is on this processor, then call routine to set the
	 particle index numbers.  */
        Grids[grid]->GridData->SetNewActiveParticleIndex(NextActiveParticleID);
    else {
 
      /* If not on this processor, then keep track of the number of new
	 star particles to ensure the active particle IDs are unique */

      NextActiveParticleID += TotalNewActiveParticleCount[grid];

    }

  }

  // Update ActiveParticleType counts as well

  int i, j, idx, nap;
  int stride = EnabledActiveParticlesCount;
  int *buffer = new int[NumberOfGrids * stride];

  for (i = 0, idx = 0; i < NumberOfGrids; i++, idx += stride)
    if (Grids[i]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      nap = Grids[i]->GridData->ReturnNumberOfActiveParticles();
	  for (j = 0; j < EnabledActiveParticlesCount; j++) {
	    buffer[idx+j] = Grids[i]->GridData->ReturnNumberOfActiveParticlesOfThisType(j);
	  }
    }
    else {
      for (j = 0; j < EnabledActiveParticlesCount; j++)
        buffer[idx + j] = 0;
    }
  
#ifdef USE_MPI
  CommunicationAllReduceValues(buffer, NumberOfGrids * stride, MPI_SUM);
#endif

  for (i = 0, idx = 0; i < NumberOfGrids; i++, idx += stride)
    for (j = 0; j < EnabledActiveParticlesCount; j++)
      Grids[i]->GridData->SetActiveParticleTypeCounts(j, buffer[idx+j]);
  
  /* Clean up. */
 
  delete [] TotalNewActiveParticleCount;
  delete [] PartialNewActiveParticleCount;
  delete [] buffer;

  LCAPERF_STOP("UpdateActiveParticleCount");
 
  return SUCCESS;
}
