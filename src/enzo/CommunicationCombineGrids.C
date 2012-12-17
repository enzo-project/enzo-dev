/***********************************************************************
/
/  COMMUNICATION ROUTINE: COMBINE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
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
 
int CommunicationCombineGrids(HierarchyEntry *OldHierarchy,
			      HierarchyEntry **NewHierarchyPointer,
			      FLOAT WriteTime, int RestartDump = FALSE)
{
 
  /* If there is only one proc, then just point the new one at the old one. */
 
  if (NumberOfProcessors == 1 || ParallelRootGridIO == TRUE) {
    *NewHierarchyPointer = OldHierarchy;
    return SUCCESS;
  }
 
  /* Otherwise generate a new hierarchy entry and proceed. */
 
  HierarchyEntry *NewHierarchy = new HierarchyEntry;
  *NewHierarchyPointer = NewHierarchy;

  int NumberOfInterpolatedFields = 0;
  int Rank, dim, Dims[MAX_DIMENSION], NewDims[MAX_DIMENSION],
      SendOffset[MAX_DIMENSION], TempDims[MAX_DIMENSION],
      StartIndex[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION], CellSize[MAX_DIMENSION];
 
  /* Compute dims, etc. for new grid assuming it fills entire domain. */
 
  OldHierarchy->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    if (dim < Rank)
      Dims[dim] -= 2*NumberOfGhostZones;
    CellSize[dim] = (Right[dim] - Left[dim])/FLOAT(Dims[dim]);
    NewDims[dim] = nint((DomainRightEdge[dim] - DomainLeftEdge[dim])/
			CellSize[dim])
                   + ((dim < Rank) ? 2*NumberOfGhostZones : 0);
  }
  if (debug)
    printf("CombineGrids: NewDims = %"ISYM" %"ISYM" %"ISYM"\n",
	   NewDims[0], NewDims[1], NewDims[2]);

  switch (OutputSmoothedDarkMatter) {
  case 1: NumberOfInterpolatedFields = 1; break;  // density
  case 2: NumberOfInterpolatedFields = 5; break;  // + rms velocity + 3-velocity
  }

  /* Generate a new grid. */
 
  NewHierarchy->GridData = new grid;
  NewHierarchy->GridData->InheritProperties(OldHierarchy->GridData);
  NewHierarchy->GridData->SetGravityParameters(
		       OldHierarchy->GridData->ReturnGravityBoundaryType());
  NewHierarchy->GridData->PrepareGrid(Rank, NewDims, DomainLeftEdge,
				      DomainRightEdge, 0);
 
  /* Loop over old grids and copy info. */
 
  HierarchyEntry *Temp = OldHierarchy;
  grid *NewGrid = NewHierarchy->GridData;
  while (Temp != NULL) {
 
    /* Compute region to send. */
 
    grid *OldGrid = Temp->GridData;
    OldGrid->ReturnGridInfo(&Rank, TempDims, Left, Right);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      SendOffset[dim] = (dim < Rank)? NumberOfGhostZones : 0;
      TempDims[dim] -= 2*SendOffset[dim];
      StartIndex[dim] = nint((Left[dim] - DomainLeftEdge[dim])/CellSize[dim])
	              + SendOffset[dim];
    }
 
    /* Copy grid region. */
 
    int RecvType = ((WriteTime < 0) && (RestartDump == FALSE)) ? 
                     NEW_ONLY : NEW_AND_OLD;
    int OldProc = OldGrid->ReturnProcessorNumber(),
        NewProc = NewGrid->ReturnProcessorNumber();
    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
//  printf("(%"ISYM"): %"ISYM" --> %"ISYM"\n", MyProcessorNumber, OldProc, NewProc);
    if (MyProcessorNumber == NewProc || MyProcessorNumber == OldProc)
      if (NewGrid->CommunicationReceiveRegion(OldGrid, OldProc, ALL_FIELDS,
			      ((WriteTime < 0) ? NEW_ONLY : NEW_AND_OLD),
			      StartIndex, TempDims, FALSE) == FAIL) {
	ENZO_FAIL("Error in grid->CommunicationReceiveRegion.\n");
      }

    /* Copy interpolated (particle) regions */

    if ((MyProcessorNumber == NewProc || MyProcessorNumber == OldProc) &&
	NumberOfInterpolatedFields > 0) {
      OldGrid->ReturnGridInfo(&Rank, TempDims, Left, Right);
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	SendOffset[dim] = 0;
	TempDims[dim] -= 2*NumberOfGhostZones;
	StartIndex[dim] = nint((Left[dim] - DomainLeftEdge[dim])/CellSize[dim])
	  + SendOffset[dim];
      }
      if (NewGrid->CommunicationReceiveRegion(OldGrid, OldProc, INTERPOLATED_FIELDS,
			      NEW_ONLY, StartIndex, TempDims, FALSE) == FAIL) {
	ENZO_FAIL("Error in grid->CommunicationReceiveRegion.\n");
      }
    }
 
    /* Copy particles. */
 
//    if (MyProcessorNumber == NewProc || MyProcessorNumber == OldProc)
      if (OldGrid->CommunicationSendParticles(NewGrid, NewProc, 0,
			    OldGrid->ReturnNumberOfParticles(), -1) == FAIL) {
	ENZO_FAIL("Error in grid->CommunicationSendParticles.\n");

      }
 
    /* Next Grid */
 
    Temp = Temp->NextGridThisLevel;
  }
//  printf("(%"ISYM"): done\n", MyProcessorNumber);
 
  /* Create a new first level of hierarchy entries that are all below the
     new one.  Below that, just point back into the old hierarchy. */
 
  NewHierarchy->ParentGrid = OldHierarchy->ParentGrid;
  NewHierarchy->NextGridNextLevel = NULL;
  NewHierarchy->NextGridThisLevel = NULL;
  Temp = OldHierarchy;
  HierarchyEntry *Previous = NULL;
  while (Temp != NULL) {
 
    HierarchyEntry *Temp2 = Temp->NextGridNextLevel;
 
    while (Temp2 != NULL) {
 
      HierarchyEntry *NewEntry = new HierarchyEntry;
      NewEntry->NextGridThisLevel = Previous;
      NewEntry->NextGridNextLevel = Temp2->NextGridNextLevel;
      NewEntry->ParentGrid        = NewHierarchy;
      NewEntry->GridData = Temp2->GridData;
      NewHierarchy->NextGridNextLevel = NewEntry;
      Previous = NewEntry;
 
      Temp2 = Temp2->NextGridThisLevel;
    }
 
    Temp = Temp->NextGridThisLevel;
  }
 
  return SUCCESS;
}
