/***********************************************************************
/
/  COMMUNICATION ROUTINE: SYNCHRONIZE NUMBER OF PARTICLES OVER ALL PROCS
/
/  written by: John Wise
/  date:       May, 2009
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

int CommunicationAllSumIntegerValues(int *Values, int Number);
 
int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],
				       int NumberOfGrids)
{

  int i;
  int *AllNumberOfParticles = new int[NumberOfGrids];

  for (i = 0; i < NumberOfGrids; i++)
    if (GridHierarchyPointer[i]->GridData->ReturnProcessorNumber() == MyProcessorNumber)
      AllNumberOfParticles[i] = GridHierarchyPointer[i]->GridData->
	ReturnNumberOfParticles();
    else
      AllNumberOfParticles[i] = 0;

  CommunicationAllSumIntegerValues(AllNumberOfParticles, NumberOfGrids);

  for (i = 0; i < NumberOfGrids; i++)
    GridHierarchyPointer[i]->GridData->SetNumberOfParticles(AllNumberOfParticles[i]);

  delete [] AllNumberOfParticles;

  return SUCCESS;
}

