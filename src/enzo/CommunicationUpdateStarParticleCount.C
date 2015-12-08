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
/
/  PURPOSE:
/    This routines keeps track of the number of new star particles.
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
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

int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids, 
					 int TotalStarParticleCountPrevious[])
{

  LCAPERF_START("UpdateStarParticleCount");

  int grid, *TotalParticleCount = new int[NumberOfGrids],
          *PartialParticleCount = new int[NumberOfGrids],
        *TotalStarParticleCount = new int[NumberOfGrids],
      *PartialStarParticleCount = new int[NumberOfGrids];
  int *buffer = new int[2*NumberOfGrids];
  int *rbuffer = new int[2*NumberOfGrids];
 
  /* Set ParticleCount to zero and record number of particles for grids
     on this processor. */
 
  for (grid = 0; grid < NumberOfGrids; grid++) {
    TotalParticleCount[grid] = 0;
    TotalStarParticleCount[grid] = 0;
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      PartialParticleCount[grid] =
	Grids[grid]->GridData->ReturnNumberOfParticles();
      PartialStarParticleCount[grid] =
	Grids[grid]->GridData->ReturnNumberOfStarParticles();
    }
    else {
      PartialParticleCount[grid] = 0;
      PartialStarParticleCount[grid] = 0;
    }
}
 
#ifdef USE_MPI
 
  /* Get counts from each processor to get total list of new particles. */
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;

  MPI_Arg GridCount = NumberOfGrids;
   
  int index;
  for (grid = 0, index = 0; grid < NumberOfGrids; index += 2, grid++) {
    buffer[index] = PartialParticleCount[grid];
    buffer[index+1] = PartialStarParticleCount[grid];
  }
  MPI_Allreduce(buffer, rbuffer, 2*GridCount,
		DataTypeInt, MPI_SUM, MPI_COMM_WORLD);
  for (grid = 0, index = 0; grid < NumberOfGrids; index += 2, grid++) {
    TotalParticleCount[grid] = rbuffer[index];
    TotalStarParticleCount[grid] = rbuffer[index+1];
  }

#ifdef UNUSED
  if (MyProcessorNumber == ROOT_PROCESSOR)
    for (grid = 0; grid < NumberOfGrids; grid++) {
      fprintf(stdout, "PartialParticleCount[%d] = %d\n", grid, PartialParticleCount[grid]); 
      fprintf(stdout, "TotalParticleCount[%d]   = %d\n", grid, TotalParticleCount[grid]);
      fprintf(stdout, "PartialStarParticleCount[%d] = %d\n", grid, PartialStarParticleCount[grid]); 
      fprintf(stdout, "TotalStarParticleCount[%d]   = %d\n", grid, TotalStarParticleCount[grid]);
      //fprintf(stdout, "TotalParticleCountPrevious[%d]   = %d\n", grid, TotalParticleCountPrevious[grid]);
      fprintf(stdout, "TotalStarParticleCountPrevious[%d]   = %d\n\n", grid, TotalStarParticleCountPrevious[grid]);
    }
#endif

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

    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) 

      /* If this grid is on this processor, then call routine to set the
	 particle index numbers.  This also updates NumberOfStarParticles. */

      Grids[grid]->GridData->SetNewParticleIndex(NumberOfStarParticles,
						 NumberOfOtherParticles);

    else {
 
      /* If not on this processor, then keep track of the number of new
	 star particles (which is the difference between the number 
	 got from the communication and what is currently stored).
	 Finally, correct the number of particles in our record. */

      NumberOfStarParticles += TotalStarParticleCount[grid] - 
	                       TotalStarParticleCountPrevious[grid];
      NumberOfOtherParticles += (TotalParticleCount[grid] - TotalStarParticleCount[grid]) - 
	                        (Grids[grid]->GridData->ReturnNumberOfParticles()
	                         - TotalStarParticleCountPrevious[grid]);
      Grids[grid]->GridData->SetNumberOfParticles(TotalParticleCount[grid]);

    }

    //printf("NumberOfStarParticles = %"ISYM"\n", NumberOfStarParticles); 

  }

#ifdef UNUSED
  fprintf(stdout, "\nin CUSPC.C \n", MetaData->NumberOfParticles); 
  fprintf(stdout, "MetaData->NumberOfParticles = %d\n", MetaData->NumberOfParticles); 
  fprintf(stdout, "NumberOfStarParticles now = %d\n", NumberOfStarParticles);
  fprintf(stdout, "NumberOfOtherParticles now = %d\n", NumberOfOtherParticles);
#endif

  /* Clean up. */
 
  delete [] TotalParticleCount;
  delete [] PartialParticleCount;
  delete [] TotalStarParticleCount;
  delete [] PartialStarParticleCount;
  delete [] buffer;
  delete [] rbuffer;
 
  LCAPERF_STOP("UpdateStarParticleCount");
 
  return SUCCESS;
}










/* Old version: not only it doesn't work when the new particles are not stars, 
   but it has assigned indices incorrectly */

#ifdef UNUSED
int CommunicationUpdateStarParticleCountOld(HierarchyEntry *Grids[],
					    TopGridData *MetaData,
					    int NumberOfGrids)
{
 
  int grid, *TotalParticleCount = new int[NumberOfGrids],
          *PartialParticleCount = new int[NumberOfGrids],
        *TotalStarParticleCount = new int[NumberOfGrids],
      *PartialStarParticleCount = new int[NumberOfGrids];
 
  /* Set ParticleCount to zero and record number of particles for grids
     on this processor. */
 
  for (grid = 0; grid < NumberOfGrids; grid++) {
    TotalParticleCount[grid] = 0;
    TotalStarParticleCount[grid] = 0;
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      PartialParticleCount[grid] =
	Grids[grid]->GridData->ReturnNumberOfParticles();
      PartialStarParticleCount[grid] =
	Grids[grid]->GridData->ReturnNumberOfStarParticles();
    }
    else {
      PartialParticleCount[grid] = 0;
      PartialStarParticleCount[grid] = 0;
    }
}
 
#ifdef USE_MPI
 
  /* Get counts from each processor to get total list of new particles. */
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;

  MPI_Arg GridCount = NumberOfGrids;
   
  MPI_Allreduce(PartialParticleCount, TotalParticleCount, GridCount,
		DataTypeInt, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(PartialStarParticleCount, TotalStarParticleCount, GridCount,
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

    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) 

      /* If this grid is on this processor, then call routine to set the
	 particle index numbers.  This also updates NumberOfStarParticles. */

      Grids[grid]->GridData->SetNewParticleIndexOld(NumberOfStarParticles,
						    MetaData->NumberOfParticles);

    else {
 
      /* If not on this processor, then keep track of the number of new
	 star particles (which is the difference between the number 
	 got from the communication and what is currently stored).
	 Finally, correct the number of particles in our record. */

      NumberOfStarParticles += TotalParticleCount[grid] -
                       Grids[grid]->GridData->ReturnNumberOfParticles();
      Grids[grid]->GridData->SetNumberOfParticles(TotalParticleCount[grid]);

    }

    //  printf("NumberOfStarParticles = %"ISYM"\n", NumberOfStarParticles); 

  }

  /* Clean up. */
 
  delete [] TotalParticleCount;
  delete [] PartialParticleCount;
  delete [] TotalStarParticleCount;
  delete [] PartialStarParticleCount;

  return SUCCESS;
}
#endif //UNUSED
