/***********************************************************************
/
/  FIND SINK PARTICLES
/
/  written by: Peng Wang
/  date:       July, 2008
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
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

#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_SHINE;
#endif

int FindSinkParticles(LevelHierarchyEntry *LevelArray[]) 
{

  int i, grid, level, nShine;
  LevelHierarchyEntry *Temp;

  nShine = 0;

  /* Initialize BirthTime to -1 to indicate that there is no particle
     stored there. */

  for (i = 0; i < MAX_SHINE_PARTICLES; i++) 
    ShiningParticles[i].BirthTime = -1;

  minStarLifetime = 1e20;

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    Temp = LevelArray[level];

    while (Temp != NULL) {

      if (Temp->GridData->GetShiningParticles(ShiningParticles, nShine)
	  == FAIL) {
	fprintf(stderr, "Error in GetShiningParticles\n");
	return FAIL;
      }

      Temp = Temp->NextGridThisLevel;

    } // END: grid
  }  // END: level

  /* Gather all shining particles on all processors */

#ifdef USE_MPI

  if (NumberOfProcessors > 1) {

    if (FirstTimeCalled) {
      MPI_Type_contiguous(sizeof(ShineParticle), MPI_BYTE, &MPI_SHINE);
      MPI_Type_commit(&MPI_SHINE);
      FirstTimeCalled = FALSE;
    }

    /* Gather a list of particle counts on each processor */

    int totalShine = 0;
    int *nCount = new int[NumberOfProcessors];
    int *displace = new int[NumberOfProcessors];
    ShineParticle *sendBuffer, *recvBuffer;

    MPI_Allgather(&nShine, 1, MPI_INT, nCount, 1, MPI_INT, MPI_COMM_WORLD);

    /* Generate displacement list. */

    for (i = 0; i < NumberOfProcessors; i++) {
      displace[i] = totalShine;
      totalShine += nCount[i];
    }

    /* If any, gather all shining particles */

    if (totalShine > 0) {

      recvBuffer = new ShineParticle[totalShine];
      sendBuffer = new ShineParticle[nShine];
      for (i = 0; i < nShine; i++)
	sendBuffer[i] = ShiningParticles[i];

      /* Share all data with all processors */

      MPI_Allgatherv(sendBuffer, nShine, MPI_SHINE,
		     recvBuffer, nCount, displace, MPI_SHINE,
		     MPI_COMM_WORLD);
      
      /* Transfer dynamic array to the static one used elsewhere */
      
      for (i = 0; i < totalShine; i++)
	ShiningParticles[i] = recvBuffer[i];

      delete [] recvBuffer;
      delete [] sendBuffer;

      NumberOfShineParticles = totalShine;

    } /* ENDIF totalShine > 0 */
  }  /* ENDIF NumberOfProcessors > 1 */
  else
    NumberOfShineParticles = nShine;

#else /* USE_MPI */
  
  NumberOfShineParticles = nShine;

#endif /* USE_MPI */

  if (debug)
    printf("FindShiningParticles: found %"ISYM" shining particles\n", 
	   NumberOfShineParticles);
  
  return SUCCESS;

}
