/***********************************************************************
/
/  RECORD TOTAL NUMBER OF ACTIVE PARTICLES
/
/  written by: Ji-hoon Kim
/  date:       October, 2009
/  modified1:  December, 2011 by John Wise -- active particle counts
/
/  PURPOSE: 
/
************************************************************************/

#ifdef USE_MPI
#endif

#include "preincludes.h"
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

void RecordTotalActiveParticleCount(HierarchyEntry *Grids[], int NumberOfGrids,
				    int TotalActiveParticleCountPrevious[])
{

  int grid, *PartialActiveParticleCountPrevious = new int[NumberOfGrids];

  for (grid = 0; grid < NumberOfGrids; grid++) {
    TotalActiveParticleCountPrevious[grid] = 0;
    if (Grids[grid]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      PartialActiveParticleCountPrevious[grid] =
	Grids[grid]->GridData->ReturnNumberOfActiveParticles();
    }
    else {
      PartialActiveParticleCountPrevious[grid] = 0;
    }
  }
 
#ifdef USE_MPI
  /* Get counts from each processor to get total list of new particles. */
 
  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;

  MPI_Arg GridCount = NumberOfGrids;
   
  MPI_Allreduce(PartialActiveParticleCountPrevious, 
		TotalActiveParticleCountPrevious, GridCount,
		DataTypeInt, MPI_SUM, MPI_COMM_WORLD);
#endif

  delete [] PartialActiveParticleCountPrevious;

}
