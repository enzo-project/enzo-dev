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

int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],
				       int NumberOfGrids)
{

  int i, idx;
  //int *AllNumberOfParticles = new int[NumberOfGrids];
  //int *AllNumberOfStars = new int[NumberOfGrids];
  int *buffer = new int[2*NumberOfGrids];

  for (i = 0, idx = 0; i < NumberOfGrids; i++, idx += 2)
    if (GridHierarchyPointer[i]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      buffer[idx] = GridHierarchyPointer[i]->GridData->
	ReturnNumberOfParticles();
      buffer[idx+1] = GridHierarchyPointer[i]->GridData->ReturnNumberOfStars();
    } else {
      buffer[idx] = 0;
      buffer[idx+1] = 0;
    }

#ifdef USE_MPI
  CommunicationAllReduceValues(buffer, 2*NumberOfGrids, MPI_SUM);
  //CommunicationAllReduceValues(AllNumberOfParticles, NumberOfGrids, MPI_SUM);
  //CommunicationAllReduceValues(AllNumberOfStars, NumberOfGrids, MPI_SUM);
#endif

  for (i = 0, idx = 0; i < NumberOfGrids; i++, idx += 2) {
    GridHierarchyPointer[i]->GridData->SetNumberOfParticles(buffer[idx]);
    GridHierarchyPointer[i]->GridData->SetNumberOfStars(buffer[idx+1]);
  }

  delete [] buffer;
  //delete [] AllNumberOfParticles;
  //delete [] AllNumberOfStars;

  return SUCCESS;
}

