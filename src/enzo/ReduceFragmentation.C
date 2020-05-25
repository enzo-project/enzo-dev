/***********************************************************************
/
/  REDUCE FRAGMENTATION FUNCTION
/
/  written by: Greg Bryan
/  date:       February, 1999
/  modified1:
/
/  PURPOSE:
/    This routine deletes the hierarchy and rereads it, hopefully
/    reducing memory fragmentation in the process.
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdio.h>
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
#include "CosmologyParameters.h"
#include "CommunicationUtilities.h"
 
/* function prototypes */
 
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		ExternalBoundary *Exterior, float *Initialdt);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
char LastFileNameWritten[MAX_LINE_LENGTH];
int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData,
			 LevelHierarchyEntry *LevelArray[], int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
 
/*  function */
 
int ReduceFragmentation(HierarchyEntry &TopGrid, TopGridData &MetaData,
			ExternalBoundary *Exterior,
			LevelHierarchyEntry *LevelArray[])
{
 
  /* Declarations. */
 
  int level;
  LevelHierarchyEntry *Previous, *Temp;
  float dummy;
  /* Delete hierarchy, level array data and grids themselves. */
 
  fprintf(stderr, "Fragmentation reduction: deleting...");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
 
      /* Delete grid itself. */
 
      delete Temp->GridData;
 
      /* Delete hierarchy entry except the top one. */
 
      if (Temp->GridHierarchyEntry != &TopGrid)
	delete Temp->GridHierarchyEntry;
 
      Previous = Temp;
      Temp = Temp->NextGridThisLevel;
 
      /* Delete previous level hierarchy entry. */
 
      delete Previous;
 
    }
  }
 
  /* Reload data (unfortunately also reads in ExternalBoundary). */
 
  fprintf(stderr, "reading %s...", LastFileNameWritten);
  CommunicationBarrier();
  if (ReadAllData(LastFileNameWritten, &TopGrid,
		  MetaData, Exterior, &dummy) == FAIL) {
    ENZO_VFAIL("Error reloading data: %s\n", LastFileNameWritten)
  }
  AddLevel(LevelArray, &TopGrid, 0);
  fprintf(stderr, "done\n");
 
  /* Set top grid boundary conditions. */
 
  Temp = LevelArray[0];
  while (Temp != NULL) {
    if (Temp->GridData->SetExternalBoundaryValues(Exterior) == FAIL) {
      ENZO_FAIL("Error in grid->SetExternalBoundaryValues.\n");
    }
    if (CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, 0)
	== FAIL) {
      ENZO_FAIL("Error in CopyOverlappingZones.\n");
    }
    Temp = Temp->NextGridThisLevel;
  }
 
  /* Rebuild the grids from level 0. */
 
  if (RebuildHierarchy(&MetaData, LevelArray, 0) == FAIL) {
    ENZO_FAIL("Error in RebuildHierarchy.\n");

  }
 
  return SUCCESS;
 
}
