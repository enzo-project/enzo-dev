  
/*****************************************************************************
 *                                                                           *
 * Copyright 2005 David Collins, Rick Wagner                                 *
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  CREATE SUBLINGLIST
/
/  written by: David Collins, Rick Wagner
/  date:       May, 2005
/  modified1:  February, 2010 by JHW (parallelized)
/
/
/  PURPOSE: Create a list, for each grid, of all grids on the next finer level
/           that share a face ONLY.
/           Note that this configuration does NOT recognize grids that are strict subgrids
/           AT ALL, so if you're running with fewer than 2 processors in each direction
/           AND periodic boundary conditions, there will be grids at the domain edge that may
/           cause conservation errors (since they'll both share a face AND be subgrids.)
/           In that case, the user will need to do some work to a.) get those grids included
/           in this list AND b.) avoid correcting that grid multiple times.  dcc.
/
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

#ifdef FAST_SIB 
int CreateSUBlingList(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level,
		      SiblingGridList SiblingList[],
		      LevelHierarchyEntry ***SUBlingList)
#else
int CreateSUBlingList(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level,
		      LevelHierarchyEntry ***SUBlingList)
#endif
{


  int i, grid1, grid2, othergrid_id, count;
  int NumberOfGrids, NumberOfChildGrids, LocalNumberOfSUBlings, TotalNumberOfSUBlings;
  LevelHierarchyEntry *Temp, *NextEntry, *LastEntry;
  HierarchyEntry *NextGrid, *OtherGrid;
  HierarchyEntry **Grids, **ChildGrids;

  /* Create grid lists for level and level+1 */

  NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  NumberOfChildGrids = GenerateGridArray(LevelArray, level+1, &ChildGrids);

  if( FluxCorrection == 0 ) {
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      (*SUBlingList)[grid1] = NULL;
    return SUCCESS;
  }


  /************************************************************************
     Create a SUBling list of ONLY the subgrid IDs for grids on this 
     processor
  ************************************************************************/

  /* May have more SUBlings, but we can reallocate the memory as
     needed.  Array of coarse and SUBling grid IDs (packed together).
     We use NumberOfSUBlings[] for the displacements. */

  const int SUBlingArrayIncr = 2*NumberOfChildGrids;
  int SUBlingArraySize = SUBlingArrayIncr;
  int *SUBlingIDs = new int[SUBlingArraySize];

  int *GlobalSUBlingIDs = NULL;  // Complete list after MPI_Allgatherv
  int *NewList = NULL;  // Temporary list when reallocating

  count = 0;
  LocalNumberOfSUBlings = 0;
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
    if (Grids[grid1]->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
#ifdef FAST_SIB
      for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++) {
	// Not sure if the sibling list includes itself, so compare pointers
	if (Grids[grid1]->GridData != SiblingList[grid1].GridList[grid2]) {
	  othergrid_id = SiblingList[grid1].GridList[grid2]->GetGridID();
	  OtherGrid = Grids[othergrid_id];
#else
      for (grid2 = 0; grid2 < NumberOfGrids; grid2++) {
	if (grid1 != grid2) {
	  OtherGrid = Grids[grid2];
#endif /* FAST_SIB */
	  if (Grids[grid1]->GridData->
	      CheckForSharedFace(OtherGrid->GridData,
				 MetaData->LeftFaceBoundaryCondition,
				 MetaData->RightFaceBoundaryCondition
				 ) == TRUE) {
	    NextGrid = OtherGrid->NextGridNextLevel;
	    while( NextGrid ){
 
	      if (Grids[grid1]->GridData->
		  CheckForSharedFace(NextGrid->GridData,
				     MetaData->LeftFaceBoundaryCondition,
				     MetaData->RightFaceBoundaryCondition
				     ) == TRUE) {

		LocalNumberOfSUBlings++;
		SUBlingIDs[count++] = grid1;
		SUBlingIDs[count++] = NextGrid->GridData->GetGridID();
		
		// If SUBlingIDs[] has become too large, allocate some
		// more memory
		if (count >= SUBlingArraySize) {
		  SUBlingArraySize += SUBlingArrayIncr;
		  NewList = new int[SUBlingArraySize];
		  memcpy(NewList, SUBlingIDs, count * sizeof(int));
		  delete [] SUBlingIDs;
		  SUBlingIDs = NewList;
		} // ENDIF gridcount >= SUBlingArraySize

	      } // ENDIF SharedFace(grid/subgrid)
	      NextGrid = NextGrid->NextGridThisLevel;
	    } // ENDWHILE NextGrid
	  } // ENDIF SharedFace(othergrid/grid)
	} // ENDIF othergrid != grid
      } // ENDFOR grid2
    } // ENDIF MyProcessorNumber == ProcessorNumber
  } // ENDFOR grid1

//  if (debug)
//    for (i = 0; i < LocalNumberOfSUBlings; i++)
//      printf("%d: %d %d\n", i, SUBlingIDs[2*i], SUBlingIDs[2*i+1]);

  /* Now we have the SUBlings on this processor, we need to gather the
     remote ones.  The serial version is easy enough -- just copy the
     array. */

  if (NumberOfProcessors == 1) {
    TotalNumberOfSUBlings = LocalNumberOfSUBlings;
    GlobalSUBlingIDs = SUBlingIDs;
  }

  // Parallel
  else {
#ifdef USE_MPI

    /* First gather the number of SUBlings on each processor to all
       processors. */

    int *SharedListCount = new int[NumberOfProcessors];
    MPI_Arg *MPI_SharedListCount = new MPI_Arg[NumberOfProcessors];
    MPI_Arg *MPI_SharedListDisplacements = new MPI_Arg[NumberOfProcessors];
    MPI_Arg SendCount, RecvCount;

    SendCount = 1;
    RecvCount = 1;
    MPI_Allgather(&LocalNumberOfSUBlings, SendCount, IntDataType,
		  SharedListCount, RecvCount, IntDataType,
		  MPI_COMM_WORLD);

  /* Create a global displacement array and get all SUBling IDs across
     processors.  Remember that we're passing 2 integers (coarse grid ID
     and SUBling ID) per entry. */

    TotalNumberOfSUBlings = 0;
    for (i = 0; i < NumberOfProcessors; i++) {
      MPI_SharedListDisplacements[i] = 2*TotalNumberOfSUBlings;
      TotalNumberOfSUBlings += SharedListCount[i];
      MPI_SharedListCount[i] = 2*SharedListCount[i];
    }

    GlobalSUBlingIDs = new int[2*TotalNumberOfSUBlings];
    SendCount = 2*LocalNumberOfSUBlings;

    MPI_Allgatherv(SUBlingIDs, SendCount, IntDataType,
		   GlobalSUBlingIDs, MPI_SharedListCount, MPI_SharedListDisplacements,
		   IntDataType, MPI_COMM_WORLD);

    /* Cleanup */

    delete [] SharedListCount;
    delete [] MPI_SharedListCount;
    delete [] MPI_SharedListDisplacements;

#endif /* USE_MPI */
  } // ENDELSE (NumberOfProcessors == 1)


  /*****************************************************************/
  /* Create the SUBlingList now that we have all of the SUBlingIDs */
  /*****************************************************************/
 
  int CoarseGrid, SUBling;
  //*SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];
  NextEntry = NULL;
  LastEntry = NULL;

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    (*SUBlingList)[grid1] = NULL;

  for (i = 0; i < TotalNumberOfSUBlings; i++) {
    CoarseGrid = GlobalSUBlingIDs[2*i];
    SUBling = GlobalSUBlingIDs[2*i+1];

    // If a SUBling already exists, insert after the head node
    if ((*SUBlingList)[CoarseGrid] != NULL) {
      LastEntry = (*SUBlingList)[CoarseGrid]->NextGridThisLevel;
      (*SUBlingList)[CoarseGrid]->NextGridThisLevel = new LevelHierarchyEntry;
      NextEntry = (*SUBlingList)[CoarseGrid]->NextGridThisLevel;
    } 

    // If no SUBling, create the head node
    else {
      LastEntry = NULL;
      (*SUBlingList)[CoarseGrid] = new LevelHierarchyEntry;
      NextEntry = (*SUBlingList)[CoarseGrid];
    }

    NextEntry->GridHierarchyEntry = ChildGrids[SUBling];
    NextEntry->GridData = ChildGrids[SUBling]->GridData;
    NextEntry->NextGridThisLevel = LastEntry;
  } // ENDFOR grid1

  /* Cleanup */

  if (NumberOfProcessors > 1)
    delete [] GlobalSUBlingIDs;
  delete [] SUBlingIDs;
  delete [] Grids;
  delete [] ChildGrids;

  return SUCCESS;
 
}
