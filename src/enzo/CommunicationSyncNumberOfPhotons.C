/***********************************************************************
/
/  COMMUNICATION ROUTINE: SYNCHRONIZE NUMBER OF PHOTONS OVER ALL PROCS
/
/  written by: John Wise
/  date:       September, 2009
/  modified:   
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
 
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

int CommunicationSyncNumberOfPhotons(LevelHierarchyEntry *LevelArray[])
{

  if (NumberOfProcessors == 1)
    return SUCCESS;

  LevelHierarchyEntry *Temp;
  int i, level, NumberOfGrids = 0;
  int *NumberOfPhotons;

  /* Count number of grids for allocation purposes */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      NumberOfGrids++;

  /* Allocate array */

  NumberOfPhotons = new int[NumberOfGrids];

  /* Get number of photons in each grid */

  i = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel, i++)
      NumberOfPhotons[i] = Temp->GridData->CountPhotonNumber();

  /* Sum up numbers */
  
  CommunicationAllSumValues(NumberOfPhotons, NumberOfGrids);

#ifdef UNUSED
  int total = 0;
  for (i = 0; i < NumberOfGrids; i++)
    total += NumberOfPhotons[i];
  printf("P%d: SyncPhotons, total NumberOfPhotons = %d\n",
	 MyProcessorNumber, total);
#endif

  /* Put back into grids */

  i = 0;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel, i++)
      Temp->GridData->SetNumberOfPhotonPackages(NumberOfPhotons[i]);

  delete [] NumberOfPhotons;

  return SUCCESS;
}

