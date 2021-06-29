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

  int i, j, idx;
  int stride = NUMBER_ENZO_PARTICLE_TYPES + MAX_ACTIVE_PARTICLE_TYPES;
  int *buffer = new int[NumberOfGrids * stride];

  for (i = 0, idx = 0; i < NumberOfGrids; i++, idx += stride)
    if (GridHierarchyPointer[i]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
      buffer[idx] = GridHierarchyPointer[i]->GridData->ReturnNumberOfParticles();
      buffer[idx+1] = GridHierarchyPointer[i]->GridData->ReturnNumberOfStars();
      buffer[idx+2] = GridHierarchyPointer[i]->GridData->ReturnNumberOfActiveParticles();
      for (j = 0; j < MAX_ACTIVE_PARTICLE_TYPES; j++) {
        if (j < EnabledActiveParticlesCount) {
          buffer[idx+3+j] = GridHierarchyPointer[i]->GridData->
            ReturnNumberOfActiveParticlesOfThisType(j);
        } else {
          buffer[idx+3+j] = 0.;
        }
      }
    } else {
      buffer[idx] = 0;
      buffer[idx+1] = 0;
      buffer[idx+2] = 0;
      for (j = 0; j < MAX_ACTIVE_PARTICLE_TYPES; j++) {
        buffer[idx+3+j] = 0;
      }
    }

#ifdef USE_MPI
  CommunicationAllReduceValues(buffer, NumberOfGrids * stride, MPI_SUM);
#endif

  for (i = 0, idx = 0; i < NumberOfGrids; i++, idx += stride) {
    GridHierarchyPointer[i]->GridData->SetNumberOfParticles(buffer[idx]);
    GridHierarchyPointer[i]->GridData->SetNumberOfStars(buffer[idx+1]);
    GridHierarchyPointer[i]->GridData->SetNumberOfActiveParticles(buffer[idx+2]);
    for (j = 0; j < MAX_ACTIVE_PARTICLE_TYPES; j++) {
      GridHierarchyPointer[i]->GridData->SetActiveParticleTypeCounts(j,
        buffer[idx+3+j]);
    }
  }

  delete [] buffer;

  return SUCCESS;
}

