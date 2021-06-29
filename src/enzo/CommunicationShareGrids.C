/***********************************************************************
/
/  COMMUNICATION ROUTINE: SHARE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
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
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
 
/* function prototypes */
 
#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_PackedGrid;
#endif
 
  
int CommunicationShareGrids(HierarchyEntry *GridHierarchyPointer[],
			    int NumberOfGrids, int ShareParticles)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
  /* Declarations. */
 
  struct PackedGrid {
    int Rank;
    int Dimension[MAX_DIMENSION];
    FLOAT LeftEdge[MAX_DIMENSION];
    FLOAT RightEdge[MAX_DIMENSION];
    int NumberOfParticles;
    int NumberOfStars;
    int ParentNumber;
  };
  int i;
 
  /* Count the subgrids on this processor. */
 
  int GridsToSend = 0;
  HierarchyEntry *Temp;
  for (i = 0; i < NumberOfGrids; i++) {
    Temp = GridHierarchyPointer[i]->NextGridNextLevel;
    while (Temp != NULL) {
      GridsToSend++;
      Temp = Temp->NextGridThisLevel;
    }
  }
 
//  printf("ShareGrids (%"ISYM"): NumberOfGrids = %"ISYM" GridsToSend = %"ISYM"\n",
//	 MyProcessorNumber, NumberOfGrids, GridsToSend);
 
  /* Allocate an array of packed subgrids and fill it out. */
 
  int Counter = 0;
  PackedGrid *SendList = NULL;
  if (GridsToSend > 0) {
    SendList = new PackedGrid[GridsToSend];
    for (i = 0; i < NumberOfGrids; i++) {
      Temp = GridHierarchyPointer[i]->NextGridNextLevel;
      while (Temp != NULL) {
	Temp->GridData->ReturnGridInfo(&SendList[Counter].Rank,
				     SendList[Counter].Dimension,
				     SendList[Counter].LeftEdge,
				     SendList[Counter].RightEdge);
	if (ShareParticles == TRUE) {
	  SendList[Counter].NumberOfParticles =
	    Temp->GridData->ReturnNumberOfParticles();
	  SendList[Counter].NumberOfStars = 
	    Temp->GridData->ReturnNumberOfStars();
	} else {
	  SendList[Counter].NumberOfParticles = 0;
	  SendList[Counter].NumberOfStars = 0;
	}
	SendList[Counter++].ParentNumber = i;
	Temp = Temp->NextGridThisLevel;
      }
    }
  }
 
  /* Allocate the array to receive subgrids. */
 
  PackedGrid *SharedList = NULL;
  int NumberOfSharedGrids = 0;
 
#ifdef USE_MPI
 
  /* Generate a new MPI data type corresponding to the PackedGrid struct. */

  MPI_Datatype DataType = MPI_BYTE;
  MPI_Arg Count = sizeof(PackedGrid);
 
  if (FirstTimeCalled) {
    MPI_Type_contiguous(Count, DataType, &MPI_PackedGrid);
    MPI_Type_commit(&MPI_PackedGrid);
    FirstTimeCalled = FALSE;
  }

  int *SharedListCount = new int[NumberOfProcessors];

  MPI_Arg *MPI_SharedListCount = new MPI_Arg[NumberOfProcessors],
          *MPI_SharedListDisplacements = new MPI_Arg[NumberOfProcessors];
 
  /* Get counts from each processor to allocate buffers. */
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Arg Sendcount = 1;
  MPI_Arg Recvcount = 1;
  MPI_Datatype DataTypeInt1 = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  MPI_Datatype DataTypeInt2 = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
   
  MPI_Allgather(&GridsToSend, Sendcount, DataTypeInt1, SharedListCount, Recvcount, DataTypeInt2, MPI_COMM_WORLD);
 
  /* Allocate buffers and generated displacement list. */
 
  for (i = 0; i < NumberOfProcessors; i++) {
    MPI_SharedListDisplacements[i] = NumberOfSharedGrids;
    NumberOfSharedGrids += SharedListCount[i];
    MPI_SharedListCount[i] = SharedListCount[i];
  }
  SharedList = new PackedGrid[NumberOfSharedGrids];
 
  /* Perform sharing operation. */

  Sendcount = GridsToSend;
   
  MPI_Allgatherv(SendList, Sendcount, MPI_PackedGrid, SharedList,
		 MPI_SharedListCount, MPI_SharedListDisplacements, MPI_PackedGrid,
		 MPI_COMM_WORLD);

#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[13] += endtime-starttime;
  counter[13] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
 
  delete [] SharedListCount;
//  delete [] SharedListDisplacements;
  delete [] MPI_SharedListCount;
  delete [] MPI_SharedListDisplacements;

 
#endif /* USE_MPI */
 
  /* Unpack the subgrids. */
 
  HierarchyEntry *PreviousGrid = NULL, *ThisGrid, *SubgridParent;
  for (i = 0; i < NumberOfSharedGrids; i++) {
 
    /* If this subgrid has a parent not on this processor, then allocate it. */
 
    SubgridParent = GridHierarchyPointer[SharedList[i].ParentNumber];
    if (SubgridParent->GridData->ReturnProcessorNumber()
	!= MyProcessorNumber) {
 
      /* create hierarchy entry */
 
      ThisGrid = new HierarchyEntry;
 
      /* set hierarchy values */
 
      if (PreviousGrid != NULL)
	if (PreviousGrid->ParentGrid != SubgridParent)
	  PreviousGrid = NULL;
 
      if (PreviousGrid == NULL)
	SubgridParent->NextGridNextLevel = ThisGrid;
      else
	PreviousGrid->NextGridThisLevel = ThisGrid;
      ThisGrid->NextGridNextLevel = NULL;
      ThisGrid->NextGridThisLevel = NULL;
      ThisGrid->ParentGrid        = SubgridParent;
      PreviousGrid = ThisGrid;
 
      /* create new grid */
 
      ThisGrid->GridData = new grid;
 
      /* set some the new grid's properties (rank, field types, etc.)
	 based on the current grid */
 
      ThisGrid->GridData->InheritProperties(SubgridParent->GridData);
 
      /* Set the new grid's positional parameters.
         (The zero indicates there are no particles (for now). */
 
      ThisGrid->GridData->PrepareGrid(SharedList[i].Rank,
				      SharedList[i].Dimension,
				      SharedList[i].LeftEdge,
				      SharedList[i].RightEdge,
				      0);
 
      if (ShareParticles == TRUE) {
	ThisGrid->GridData->SetNumberOfParticles(SharedList[i].NumberOfParticles);
	ThisGrid->GridData->SetNumberOfStars(SharedList[i].NumberOfStars);
      }

      if (SubgridParent != NULL && ShareParticles == TRUE) {

	  SubgridParent->GridData->SetNumberOfParticles
	    (SubgridParent->GridData->ReturnNumberOfParticles() - 
	     SharedList[i].NumberOfParticles);
	  SubgridParent->GridData->SetNumberOfStars
	    (SubgridParent->GridData->ReturnNumberOfStars() - 
	     SharedList[i].NumberOfStars);

      } // ENDIF parent exists and sharing particles/stars
 
      ThisGrid->GridData->SetProcessorNumber
	(SubgridParent->GridData->ReturnProcessorNumber());
 
    }
  } // end: loop over shared grids
 
  /* CleanUp. */
 
  delete [] SendList;
  delete [] SharedList;
 
  return SUCCESS;
}
