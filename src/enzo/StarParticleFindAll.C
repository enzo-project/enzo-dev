/***********************************************************************
/
/  FIND ALL STAR PARTICLES OVER ALL PROCESSORS
/
/  written by: John Wise
/  date:       March, 2009
/  modified1: August 2021 by Ka Hou Leong 
              (fixed MPI issue and lack of memory)
/
/  PURPOSE: First synchronizes particle information in the normal and 
/           star particles.  Then we make a global particle list, which
/           simplifies adding the feedback to the all of the grids.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "CommunicationUtilities.h"

#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_STAR;
#endif

/* Global communication buffers to avoid reallocations and 
   this memory fragmentation */

StarBuffer *recvBuffer = NULL, *sendBuffer = NULL;
int recvBufferSize = 0, sendBufferSize = 0;

void InsertStarAfter(Star * &Node, Star * &NewNode);
void DeleteStarList(Star * &Node);
Star* StarBufferToList(StarBuffer *buffer, int n);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

int StarParticleFindAll(LevelHierarchyEntry *LevelArray[], Star *&AllStars)
{

  int i, level, GridNum, TotalNumberOfStars, LocalNumberOfStars;
  int SavedP3IMFCalls;
  Star *LocalStars = NULL, *GridStars = NULL, *cstar = NULL, *lstar = NULL;
  HierarchyEntry **Grids;
  int NumberOfGrids, *NumberOfStarsInGrids;

  minStarLifetime = 1e20;
  TotalNumberOfStars = 0;
  LocalNumberOfStars = 0;

  if (AllStars != NULL)
    DeleteStarList(AllStars);

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
    NumberOfStarsInGrids = new int[NumberOfGrids];

    for (GridNum = 0; GridNum < NumberOfGrids; GridNum++) {

      // First update any existing star particles (e.g. position,
      // velocity)
      if (Grids[GridNum]->GridData->UpdateStarParticles(level) == FAIL) {
		ENZO_FAIL("Error in grid::UpdateStarParticles.");
      }

      // Then find any newly created star particles
      if (Grids[GridNum]->GridData->FindNewStarParticles(level) == FAIL) {
		ENZO_FAIL("Error in grid::FindNewStarParticles.");
      }

      // Now copy any stars into the local linked list
      NumberOfStarsInGrids[GridNum] = 0;
      GridStars = Grids[GridNum]->GridData->ReturnStarPointer();
      while (GridStars != NULL) {
	cstar = GridStars->copy();
	InsertStarAfter(LocalStars, cstar);
	GridStars = GridStars->NextStar;
	NumberOfStarsInGrids[GridNum]++;
      } // ENDWHILE stars

      LocalNumberOfStars += NumberOfStarsInGrids[GridNum];

    } // ENDFOR grids

    /* Synchronize number of stars across processors */

#ifdef USE_MPI
//    CommunicationAllReduceValues(NumberOfStarsInGrids, NumberOfGrids, MPI_MAX);
//    for (GridNum = 0; GridNum < NumberOfGrids; GridNum++)
//      Grids[GridNum]->GridData->SetNumberOfStars(NumberOfStarsInGrids[GridNum]);
#endif

    delete [] Grids;
    delete [] NumberOfStarsInGrids;

  } // ENDFOR level

  /***********************************************/
  /*                                             */
  /* Gather all star particles on all processors */
  /*                                             */
  /***********************************************/

  if (NumberOfProcessors > 1) {

#ifdef USE_MPI
    if (FirstTimeCalled) {
      MPI_Type_contiguous(sizeof(StarBuffer), MPI_BYTE, &MPI_STAR);
      MPI_Type_commit(&MPI_STAR);
      FirstTimeCalled = FALSE;
    }

    /* Gather a list of particle counts on each processor */

    Eint32 *nCount = new Eint32[NumberOfProcessors];
    Eint32 *displace = new Eint32[NumberOfProcessors];

    MPI_Allgather(&LocalNumberOfStars, 1, MPI_INT, nCount, 1, MPI_INT, 
		  MPI_COMM_WORLD);

    /* Generate displacement list. */

    for (i = 0; i < NumberOfProcessors; i++) {
      displace[i] = TotalNumberOfStars;
      TotalNumberOfStars += nCount[i];
    }

    /* If any, gather all shining particles */

    if (TotalNumberOfStars > 0) {
      if (recvBufferSize > 2 * ceil_log2(TotalNumberOfStars))
      {
       // Avoiding recvBuffer occurs memoeries which exceed 2 times of the powers of buffer space which has minimum space to contain TotalNumberOfStars.
        delete [] recvBuffer;
        recvBufferSize = 0;
      }
      if (TotalNumberOfStars > recvBufferSize) {
        if (recvBufferSize > 0)
        {
            recvBufferSize = ceil_log2(TotalNumberOfStars);
            delete [] recvBuffer;
        }
        else recvBufferSize = ceil_log2(TotalNumberOfStars);
        recvBuffer = new StarBuffer[recvBufferSize];
      }
      if ((LocalNumberOfStars > 0) && (sendBufferSize > 2 * ceil_log2(LocalNumberOfStars)))
      {
       // Avoiding sendBuffer occurs memoeries which exceed 2 times of the powers of buffer space which has minimum space to contain LocalNumberOfStars.
        delete [] sendBuffer;
        sendBufferSize = 0;
      }
      if (LocalNumberOfStars > sendBufferSize) 
      { 
        if(sendBufferSize > 0)
        {
            sendBufferSize = ceil_log2(LocalNumberOfStars);
            delete [] sendBuffer;
        }
        else sendBufferSize = ceil_log2(LocalNumberOfStars); 
        sendBuffer = new StarBuffer[sendBufferSize];
      }
      if (LocalNumberOfStars > 0)
      {
          LocalStars->StarListToBuffer(sendBuffer, LocalNumberOfStars);
      }
      else
      { 
        // Due to No local star, reinitialise sendbuffer.
        // release memories and reset sendBufferSize.
        if(sendBufferSize > 0)
        {
            delete [] sendBuffer;
            sendBufferSize = 0;
        }
        sendBuffer = NULL;
      }

      /* Share all data with all processors */

      MPI_Allgatherv(sendBuffer, LocalNumberOfStars, MPI_STAR,
		     recvBuffer, nCount, displace, MPI_STAR,
		     MPI_COMM_WORLD);

      AllStars = StarBufferToList(recvBuffer, TotalNumberOfStars);

      /* Re-assign CurrentGrid pointers to local particles */

      int i0, i1;
      cstar = AllStars;
      lstar = LocalStars;
      for (i = 0; i < TotalNumberOfStars; i++) {
	i0 = displace[MyProcessorNumber];
	i1 = (MyProcessorNumber < NumberOfProcessors-1) ? 
	  displace[MyProcessorNumber+1] : TotalNumberOfStars;

	// local processors from Allgatherv
	if (i >= i0 && i < i1) {
	  cstar->AssignCurrentGrid(lstar->ReturnCurrentGrid());
	  lstar = lstar->NextStar;
	} // ENDIF local
	else
	  cstar->AssignCurrentGrid(NULL);

	cstar = cstar->NextStar;

      } // ENDFOR stars

      DeleteStarList(LocalStars);

    } /* ENDIF TotalNumberOfStars > 0 */

    delete [] nCount;
    delete [] displace;
#endif /* USE_MPI */
  }  /* ENDIF NumberOfProcessors > 1 */
  else {
    TotalNumberOfStars = LocalNumberOfStars;
    AllStars = LocalStars;
  }

  /* Find minimum stellar lifetime */
  
  for (cstar = AllStars; cstar; cstar = cstar->NextStar)
    if (cstar->ReturnMass() > 1e-9)
      minStarLifetime = min(minStarLifetime, cstar->ReturnLifetime());

  /* Store in global variable */
  
  G_TotalNumberOfStars = TotalNumberOfStars;

  return SUCCESS;

}
