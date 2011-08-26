/***********************************************************************
/
/  COMMUNICATION ROUTINE: COLLECT PARTICLES AFTER REBUILDING HIERARCHY
/
/  written by: John Wise
/  date:       May, 2009
/  modified:   July, 2009 by John Wise to collect stars as well
/
/  PURPOSE:
/
/  NOTE: communication modeled after the optimized version of 
/        CommunicationTransferParticles.
/  NOTE1: If MoveStars is false, the temp. memory is still allocated
/         and deallocated but nothing is moved.  This is a bit
/         inefficient but reduced the number of if-statements.
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
 
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
#include "SortCompareFunctions.h"
void my_exit(int status);
 
// function prototypes
 
Eint32 compare_grid(const void *a, const void *b);
Eint32 compare_star_grid(const void *a, const void *b);
int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],
				       int NumberOfGrids);
int CommunicationShareParticles(int *NumberToMove, particle_data* &SendList,
				int &NumberOfReceives,
				particle_data* &SharedList);
int CommunicationShareStars(int *NumberToMove, star_data* &SendList,
			    int &NumberOfReceives, star_data* &SharedList);

#define NO_DEBUG_CCP
#define GRIDS_PER_LOOP 100000
#define PARTICLES_PER_LOOP 100000000
 
int CommunicationCollectParticles(LevelHierarchyEntry *LevelArray[],
				  int level, bool ParticlesAreLocal, 
				  bool SyncNumberOfParticles, 
				  bool MoveStars, int CollectMode)
{

  /* Create pointer arrays and count grids */

  int NumberOfGrids = 0, NumberOfSubgrids = 0;
  HierarchyEntry *SubgridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
  HierarchyEntry *GridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
  HierarchyEntry *Subgrid;
  LevelHierarchyEntry *Temp;

  Temp = LevelArray[level];
  while (Temp != NULL) {
    GridHierarchyPointer[NumberOfGrids++] = Temp->GridHierarchyEntry;
    Temp = Temp->NextGridThisLevel;
  }

  if (level < MAX_DEPTH_OF_HIERARCHY-1) {
    Temp = LevelArray[level+1];
    while (Temp != NULL) {
      Temp->GridData->SetGridID(NumberOfSubgrids);
      SubgridHierarchyPointer[NumberOfSubgrids++] = Temp->GridHierarchyEntry;
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* 
     TWO MODES: 
     1) Move particles to subgrids by either keeping them local or sending 
        them to the home processor.
     2) Move the rest from their "local" processor (i.e. the processor 
        of the previous grid it existed on) to their "home" processor 
	(grid->ProcessorNumber).
  */

  particle_data *SendList = NULL;
  particle_data *SharedList = NULL;
  star_data *StarSendList = NULL;
  star_data *StarSharedList = NULL;

  int NumberOfReceives, StarNumberOfReceives, TotalNumber;
  int *NumberToMove = new int[NumberOfProcessors];
  int *StarsToMove = new int[NumberOfProcessors];

  int proc, i, j, k, jstart, jend, ThisID;
  int particle_data_size, star_data_size;
  int Zero = 0;

  /*********************************************************************/
  // First move all particles that exist in new subgrids

  if (CollectMode == SUBGRIDS_LOCAL || CollectMode == SUBGRIDS_GLOBAL ||
      CollectMode == ALLGRIDS) {

    // Nothing to do, but still sync number of particles in grids.
    if (NumberOfSubgrids == 0) {
      CommunicationSyncNumberOfParticles(GridHierarchyPointer, NumberOfGrids);
      delete [] NumberToMove;
      delete [] StarsToMove;
      return SUCCESS;
    }

    grid** SubgridPointers = new grid*[NumberOfSubgrids];
    bool KeepLocal = (CollectMode == SUBGRIDS_LOCAL);

#ifdef DEBUG_CCP
    for (i = 0; i < NumberOfGrids; i++)
      printf("CCP[P%"ISYM"B]: grid %"ISYM", %"ISYM" proc, %"ISYM" particles, "
	     "%"ISYM" stars\n",
	     MyProcessorNumber, i,
	     GridHierarchyPointer[i]->GridData->ReturnProcessorNumber(),
	     GridHierarchyPointer[i]->GridData->ReturnNumberOfParticles(),
	     GridHierarchyPointer[i]->GridData->ReturnNumberOfStars());
    for (i = 0; i < NumberOfSubgrids; i++)
      printf("CCP[P%"ISYM"B]: subgrid %"ISYM", %"ISYM" proc, %"ISYM" particles, "
	     "%"ISYM" stars\n",
	     MyProcessorNumber, i,
	     SubgridHierarchyPointer[i]->GridData->ReturnProcessorNumber(),
	     SubgridHierarchyPointer[i]->GridData->ReturnNumberOfParticles(),
	     SubgridHierarchyPointer[i]->GridData->ReturnNumberOfStars());
#endif /* DEBUG_CCP */

    for (i = 0; i < NumberOfProcessors; i++) {
      NumberToMove[i] = 0;
      StarsToMove[i] = 0;
    }

    for (j = 0; j < NumberOfSubgrids; j++)
      SubgridPointers[j] = SubgridHierarchyPointer[j]->GridData;

    /* In each grid, mark subgrid field and copy out particles in those
       subgrids. */

    int ZeroOnAllProcs = (ParticlesAreLocal) ? FALSE : TRUE;

    /* Count number of particles to move first to allocate memory */

    for (j = 0; j < NumberOfGrids; j++)
      if (GridHierarchyPointer[j]->NextGridNextLevel != NULL) {

	if (GridHierarchyPointer[j]->GridData->ReturnNumberOfParticles() == 0 &&
	    GridHierarchyPointer[j]->GridData->ReturnNumberOfStars() == 0)
	  continue;

	GridHierarchyPointer[j]->GridData->
	  ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD, 1.0, 
				   ZeroOnAllProcs);
 
	for (Subgrid = GridHierarchyPointer[j]->NextGridNextLevel;
	     Subgrid; Subgrid = Subgrid->NextGridThisLevel) {
	  ThisID = Subgrid->GridData->GetGridID();
	  GridHierarchyPointer[j]->GridData->ZeroSolutionUnderSubgrid
	    (Subgrid->GridData, ZERO_UNDER_SUBGRID_FIELD, float(ThisID+1),
	     ZeroOnAllProcs);
	}

	if (MoveStars)
	  GridHierarchyPointer[j]->GridData->TransferSubgridStars
	    (SubgridPointers, NumberOfSubgrids, StarsToMove, Zero, Zero, 
	     StarSendList, KeepLocal, ParticlesAreLocal, COPY_OUT); 

	GridHierarchyPointer[j]->GridData->TransferSubgridParticles
	    (SubgridPointers, NumberOfSubgrids, NumberToMove, Zero, Zero, 
	     SendList, KeepLocal, ParticlesAreLocal, COPY_OUT, FALSE, TRUE);
 
      } // ENDIF subgrids exist

    /* Now allocate the memory once and store the particles to move */

    TotalNumber = 0;
    for (j = 0; j < NumberOfProcessors; j++) {
      TotalNumber += NumberToMove[j];
      NumberToMove[j] = 0;  // Zero-out to use in the next step
    }
    SendList = new particle_data[TotalNumber];

    for (j = 0; j < NumberOfGrids; j++)
      if (GridHierarchyPointer[j]->NextGridNextLevel != NULL) {

	if (GridHierarchyPointer[j]->GridData->ReturnNumberOfParticles() == 0 &&
	    GridHierarchyPointer[j]->GridData->ReturnNumberOfStars() == 0)
	  continue;

	GridHierarchyPointer[j]->GridData->TransferSubgridParticles
	    (SubgridPointers, NumberOfSubgrids, NumberToMove, Zero, Zero, 
	     SendList, KeepLocal, ParticlesAreLocal, COPY_OUT, FALSE, FALSE);
 
      } // ENDIF subgrids exist

    /* Now we have a list of particles to move to subgrids.  If
       specified, we communicate them with all processors.  If not,
       just sort them by subgrid number. */

    if (KeepLocal) {

      /* particles first */

      SharedList = SendList;
      NumberOfReceives = NumberToMove[MyProcessorNumber];
      //particle_data_size = sizeof(particle_data);
      //qsort(SharedList, NumberOfReceives, particle_data_size, compare_grid);
      std::sort(SharedList, SharedList+NumberOfReceives, cmp_grid());

      /* stars second */

      if (MoveStars) {
	StarSharedList = StarSendList;
	StarNumberOfReceives = StarsToMove[MyProcessorNumber];
	//star_data_size = sizeof(star_data);
	//qsort(StarSharedList, StarNumberOfReceives, star_data_size, 
	//      compare_star_grid);
	std::sort(StarSharedList, StarSharedList+StarNumberOfReceives,
		  cmp_star_grid());
      }

    } // ENDIF local
    else {

      CommunicationShareParticles(NumberToMove, SendList, NumberOfReceives,
				  SharedList);
      if (MoveStars)
	CommunicationShareStars(StarsToMove, StarSendList, StarNumberOfReceives,
				StarSharedList);

    } // ENDELSE local

    /*******************************************************************/
    /****************** Copy particles back to grids. ******************/
    /*******************************************************************/

    jstart = 0;
    jend = 0;

    // Copy shared particles to grids, if any

    if (NumberOfReceives > 0)
      for (j = 0; j < NumberOfSubgrids && jend < NumberOfReceives; j++) {
	while (SharedList[jend].grid <= j) {
	  jend++;
	  if (jend == NumberOfReceives) break;
	}

	SubgridPointers[j]->TransferSubgridParticles
	  (SubgridPointers, NumberOfSubgrids, NumberToMove, jstart, jend, 
	   SharedList, KeepLocal, ParticlesAreLocal, COPY_IN);

	jstart = jend;
      } // ENDFOR grids

    /*******************************************************************/
    /******************** Copy stars back to grids. ********************/
    /*******************************************************************/

    jstart = 0;
    jend = 0;

    // Copy shared stars to grids, if any

    if (MoveStars) {

    if (StarNumberOfReceives > 0)
      for (j = 0; j < NumberOfSubgrids && jend < StarNumberOfReceives; j++) {
	while (StarSharedList[jend].grid <= j) {
	  jend++;
	  if (jend == StarNumberOfReceives) break;
	}

	SubgridPointers[j]->TransferSubgridStars
	  (SubgridPointers, NumberOfSubgrids, StarsToMove, jstart, jend, 
	   StarSharedList, KeepLocal, ParticlesAreLocal, COPY_IN);

	jstart = jend;
      } // ENDFOR grids

    } // ENDIF MoveStars

    /************************************************************************
       If the particles and stars are only on the grid's host
       processor, set number of particles so everybody agrees. 
    ************************************************************************/

    if ((!KeepLocal && NumberOfProcessors > 1) || 
	(ParticlesAreLocal && SyncNumberOfParticles)) {

      CommunicationSyncNumberOfParticles(GridHierarchyPointer, NumberOfGrids);
      CommunicationSyncNumberOfParticles(SubgridHierarchyPointer, NumberOfSubgrids);

    } // ENDIF synchronize particle counts

#ifdef DEBUG_CCP
    for (i = 0; i < NumberOfGrids; i++)
      printf("CCP[P%"ISYM"A]: grid %"ISYM", %"ISYM" proc, %"ISYM" particles, "
	     "%"ISYM" stars\n",
	     MyProcessorNumber, i,
	     GridHierarchyPointer[i]->GridData->ReturnProcessorNumber(),
	     GridHierarchyPointer[i]->GridData->ReturnNumberOfParticles(),
	     GridHierarchyPointer[i]->GridData->ReturnNumberOfStars());
    for (i = 0; i < NumberOfSubgrids; i++)
      printf("CCP[P%"ISYM"A]: subgrid %"ISYM", %"ISYM" proc, %"ISYM" particles, "
	     "%"ISYM" stars\n",
	     MyProcessorNumber, i,
	     SubgridHierarchyPointer[i]->GridData->ReturnProcessorNumber(),
	     SubgridHierarchyPointer[i]->GridData->ReturnNumberOfParticles(),
	     SubgridHierarchyPointer[i]->GridData->ReturnNumberOfStars());
#endif /* DEBUG_CCP */

    /* Cleanup. */

    if (SendList != SharedList)
      delete [] SharedList;
    delete [] SendList;
    if (StarSendList != StarSharedList)
      delete [] StarSendList;
    delete [] StarSharedList;
    delete [] SubgridPointers;

    SendList = NULL;
    SharedList = NULL;
    StarSendList = NULL;
    StarSharedList = NULL;

  } // ENDIF subgrid mode

  /*********************************************************************/
  //  Send all leftover particles and stars on the grids to their home
  //  processors.

  /* Create a list of particle moves and delete the particles not on
     the right processor. */

  if ((CollectMode == SIBLINGS_ONLY || CollectMode == ALLGRIDS) &&
      NumberOfProcessors > 1) {

    int StartGrid, EndGrid, StartNum, TotalNumberToMove, AllMovedParticles;
    int TotalStarsToMove, AllMovedStars;
    StartGrid = 0;
    EndGrid = 0;
    //for (StartGrid = 0; StartGrid < NumberOfGrids; StartGrid += GRIDS_PER_LOOP) {
    while (EndGrid < NumberOfGrids) {

      // Make a cumulative sum of particles until we're above
      // PARTICLES_PER_LOOP, then the EndGrid will be the minimum grid
      // number on all processors.  Remember that NumberOfParticles is
      // still the number of particles stored on this processor
      // because we're going to collect all of the particles onto the
      // host processor in this routine.

      StartGrid = EndGrid;
      TotalNumberToMove = 0;
      while (TotalNumberToMove < PARTICLES_PER_LOOP && 
	     EndGrid < NumberOfGrids) {
	if (GridHierarchyPointer[EndGrid]->GridData->ReturnProcessorNumber() != 
	    MyProcessorNumber) 
	  TotalNumberToMove += GridHierarchyPointer[EndGrid]->GridData->
	    ReturnNumberOfParticles();
	EndGrid++;
      }

#ifdef USE_MPI
      CommunicationAllReduceValues(&EndGrid, 1, MPI_MIN);
#endif

      /* Count particles and stars to move */

      TotalNumberToMove = 0;
      TotalStarsToMove = 0;
      for (i = StartGrid; i < EndGrid; i++)
	if (GridHierarchyPointer[i]->GridData->ReturnProcessorNumber() != 
	    MyProcessorNumber) {
	  TotalNumberToMove += GridHierarchyPointer[i]->GridData->
	    ReturnNumberOfParticles();
	  TotalStarsToMove += GridHierarchyPointer[i]->GridData->
	    ReturnNumberOfStars();
	}

      AllMovedParticles = TotalNumberToMove;
      AllMovedStars = TotalStarsToMove;
#ifdef USE_MPI
      int ibuffer[2];
      if (NumberOfProcessors > 1) {
	ibuffer[0] = AllMovedParticles;
	ibuffer[1] = AllMovedStars;
	CommunicationAllReduceValues(ibuffer, 2, MPI_SUM);
	AllMovedParticles = ibuffer[0];
	AllMovedStars = ibuffer[1];
      }
#endif

#ifdef DEBUG_CCP
      printf("CCP[%d]: Collecting a total of %"ISYM" (%"ISYM" local) "
	     "particles over grids %"ISYM"->%"ISYM".\n", 
	     MyProcessorNumber, AllMovedParticles, TotalNumberToMove,
	     StartGrid, EndGrid-1);  

    for (i = StartGrid; i < EndGrid; i++)
      printf("CCP[P%"ISYM"BB]: grid %"ISYM", %"ISYM" proc, %"ISYM" particles\n",
	     MyProcessorNumber, i,
	     GridHierarchyPointer[i]->GridData->ReturnProcessorNumber(),
	     GridHierarchyPointer[i]->GridData->ReturnNumberOfParticles());
#endif /* DEBUG_CCP */

    /* Count the number of particles needed to move */

    SendList = new particle_data[TotalNumberToMove];
    StarSendList = new star_data[TotalStarsToMove];

    for (i = 0; i < NumberOfProcessors; i++) {
      NumberToMove[i] = 0;
      StarsToMove[i] = 0;
    }

    StartNum = 0;
    for (j = StartGrid; j < EndGrid; j++)
      GridHierarchyPointer[j]->GridData->CollectParticles
	(j, NumberToMove, StartNum, Zero, SendList, COPY_OUT);

    if (MoveStars) {
      StartNum = 0;
      for (j = StartGrid; j < EndGrid; j++)
	GridHierarchyPointer[j]->GridData->CollectStars
	  (j, StarsToMove, StartNum, Zero, StarSendList, COPY_OUT);
    } // ENDIF MoveStars

    if (StarsToMove[MyProcessorNumber] > 0)
      printf("CCP: moving stars!\n");

    /* Share the particle move list */

    NumberOfReceives = 0;
    StarNumberOfReceives = 0;
    CommunicationShareParticles(NumberToMove, SendList, NumberOfReceives,
				SharedList);
    if (MoveStars)
      CommunicationShareStars(StarsToMove, StarSendList, StarNumberOfReceives,
			      StarSharedList);
  
    /*******************************************************************/
    /****************** Copy particles back to grids. ******************/
    /*******************************************************************/

    jstart = 0;
    jend = 0;

    // Copy shared particles to grids, if any
    if (NumberOfReceives > 0)
      for (j = StartGrid; j < EndGrid && jend < NumberOfReceives; j++) {
	while (SharedList[jend].grid <= j) {
	  jend++;
	  if (jend == NumberOfReceives) break;
	}

	GridHierarchyPointer[j]->GridData->CollectParticles
	  (j, NumberToMove, jstart, jend, SharedList, COPY_IN);

	jstart = jend;
      } // ENDFOR grids

#ifdef DEBUG_CCP
    for (i = StartGrid; i < EndGrid; i++)
      printf("CCP[P%"ISYM"CC]: grid %"ISYM", %"ISYM" proc, %"ISYM" particles\n",
	     MyProcessorNumber, i,
	     GridHierarchyPointer[i]->GridData->ReturnProcessorNumber(),
	     GridHierarchyPointer[i]->GridData->ReturnNumberOfParticles());
#endif /* DEBUG_CCP */

    /*******************************************************************/
    /******************** Copy stars back to grids. ********************/
    /*******************************************************************/

    jstart = 0;
    jend = 0;
    
    if (MoveStars) {

    // Copy shared stars to grids, if any
    if (StarNumberOfReceives > 0)
      for (j = StartGrid; j < EndGrid && jend < StarNumberOfReceives; j++) {
	while (StarSharedList[jend].grid <= j) {
	  jend++;
	  if (jend == StarNumberOfReceives) break;
	}

	GridHierarchyPointer[j]->GridData->CollectStars
	  (j, StarsToMove, jstart, jend, StarSharedList, COPY_IN);

	jstart = jend;
      } // ENDFOR grids

    } // ENDIF MoveStars

#ifdef DEBUG_CCP
    for (i = StartGrid; i < EndGrid; i++)
      printf("CCP[P%"ISYM"DD]: grid %"ISYM", %"ISYM" proc, %"ISYM" particles\n",
	     MyProcessorNumber, i,
	     GridHierarchyPointer[i]->GridData->ReturnProcessorNumber(),
	     GridHierarchyPointer[i]->GridData->ReturnNumberOfParticles());
#endif /* DEBUG_CCP */

    /* Cleanup. */

    if (SendList != SharedList)
      delete [] SendList;
    delete [] SharedList;

    if (StarSendList != StarSharedList)
      delete [] StarSendList;
    delete [] StarSharedList;

    } // ENDFOR grid batches

    /************************************************************************
       If the particles and stars are only on the grid's host
       processor, set number of particles so everybody agrees. 
    ************************************************************************/

    if (SyncNumberOfParticles)
      CommunicationSyncNumberOfParticles(GridHierarchyPointer, NumberOfGrids);
    else {
      for (i = 0; i < NumberOfGrids; i++)
	if (MyProcessorNumber != 
	    GridHierarchyPointer[i]->GridData->ReturnProcessorNumber()) {
	  GridHierarchyPointer[i]->GridData->SetNumberOfParticles(0);
	  if (MoveStars)
	    GridHierarchyPointer[i]->GridData->SetNumberOfStars(0);
	}
    }

  } // ENDIF sibling grids and multi-processor

  delete [] NumberToMove;
  delete [] StarsToMove;

  return SUCCESS;
}

