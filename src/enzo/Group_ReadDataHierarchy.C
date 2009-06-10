/***********************************************************************
/
/  READ IN THE DATA HIERARCHY (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Alexei Kritsuk, Jan 2004 now reads RandomForcing fields //AK
/  modified2:  Robert Harkness, Jan 2006
/              Read Unigrid Grid-to-MPI task map
/              This is necessary only for ultra-large grids on
/              node memory challenged systems
/  modified3:  Robert Harkness, Jan 2007
/              For in-core driver 
/
/  PURPOSE:
/
************************************************************************/
 
// This function reads in the data hierarchy (TopGrid)

#include <hdf5.h> 
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <map>

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

extern std::map<HierarchyEntry *, int> OriginalGridID;
#ifdef USE_HDF5_GROUPS
extern std::map<HierarchyEntry *, int> OriginalTaskID;
 #endif

/* function prototypes */
 
static int ReadDataGridCounter = 0;



 
int Group_ReadDataHierarchy(FILE *fptr, HierarchyEntry *Grid, int GridID,
		           HierarchyEntry *ParentGrid, hid_t file_id)
{
 
  int TestGridID, NextGridThisLevelID, NextGridNextLevelID;
  int Task;
 
  /* Read header info for this grid */
 
  if (fscanf(fptr, "\nGrid = %"ISYM"\n", &TestGridID) != 1) {
    fprintf(stderr, "Error reading Grid # in grid %"ISYM".\n", GridID);
    ENZO_FAIL("Error in: "__FILE__);
  }
  if (TestGridID != GridID) {
    fprintf(stderr, "Unexpected GridID = %"ISYM" while reading grid %"ISYM".\n",
	    TestGridID, GridID);
    ENZO_FAIL("Error in: "__FILE__);
  }

  fscanf(fptr, "Task = %"ISYM"\n", &Task);
  Task = Task % NumberOfProcessors;

//  if ( MyProcessorNumber == 0 )
//    fprintf(stderr, "Reading Grid %"ISYM" assigned to Task %"ISYM"\n", TestGridID, Task);
 
  /* Create new grid and fill out hierarchy entry. */
 
  Grid->GridData          = new grid;
  Grid->NextGridThisLevel = NULL;
  Grid->NextGridNextLevel = NULL;
  Grid->ParentGrid        = ParentGrid;

#ifdef SINGLE_HDF5_OPEN_ON_INPUT

// Explicit task mapping cannot be used if SINGLE_HDF5_OPEN_ON_INPUT is in effect

#ifdef ENABLE_TASKMAP
  fprintf(stderr, "Task map cannot be used with SINGLE_HDF5_OPEN_ON_INPUT!\n");
  ENZO_FAIL("Error in: "__FILE__);
#endif

// Use grid-to-processor mapping from last dump

//  if ( MyProcessorNumber == 0 )
//    fprintf(stderr, "Using dumped task assignment: GridID = %"ISYM"  MPI Task = %"ISYM"\n", GridID, Task);

    Grid->GridData->SetProcessorNumber(Task);

#else

#ifdef ENABLE_TASKMAP

//  if ( MyProcessorNumber == 0 )
//    fprintf(stderr, "Task map assignment: GridID = %"ISYM"  MPI Task = %"ISYM"\n", GridID, TaskMap[GridID-1]);

  Grid->GridData->SetProcessorNumber(TaskMap[GridID-1]);

#else

// Use grid-to-processor mapping from last dump (or simple cyclic map)

//  if ( MyProcessorNumber == 0 )
//    fprintf(stderr, "Using dumped task assignment: GridID = %"ISYM"  MPI Task = %"ISYM"\n", GridID, Task);

    Grid->GridData->SetProcessorNumber(Task);

#endif

#endif /* SINGLE OPEN */


  int NewProc = ReadDataGridCounter % NumberOfProcessors;
  int ProcMap = ABS(NewProc - NumberOfProcessors) % NumberOfProcessors;

#ifdef USE_CYCLIC_CPU_DISTRIBUTION

  Grid->GridData->SetProcessorNumber(NewProc);

  if ( MyProcessorNumber == 0 )
    fprintf(stderr, "TASKMAP DISABLED: Grid->Processor assignment:  GridID = %"ISYM"  MPI Task = %"ISYM"\n", GridID, NewProc);

#endif

#ifdef USE_PERMUTED_CPU_DISTRIBUTION

  Grid->GridData->SetProcessorNumber(ProcMap);

  if ( MyProcessorNumber == 0 )
    fprintf(stderr, "TASKMAP DISABLED: Grid->Processor assignment:  GridID = %"ISYM"  MPI Task = %"ISYM"\n", GridID, ProcMap);

#endif

  if(!LoadGridDataAtStart){
    // For stand alone analysis tools,
    // use the cyclic method, since there is no mapping between
    // the original data distribution, and the current
    // number of processors
    Grid->GridData->SetProcessorNumber(NewProc);
  }

  ReadDataGridCounter++;

  /* Read grid data for this grid. */

  if(LoadGridDataAtStart){ 
    if (Grid->GridData->Group_ReadGrid(fptr, GridID, file_id) == FAIL) {
      fprintf(stderr, "Error in grid->Group_ReadGrid (grid %"ISYM").\n", GridID);
      ENZO_FAIL("Error in: "__FILE__);
    }
  }else{
    if (Grid->GridData->Group_ReadGrid(fptr, GridID, file_id, TRUE, FALSE) == FAIL) {
      fprintf(stderr, "Error in grid->Group_ReadGrid (grid %"ISYM").\n", GridID);
      ENZO_FAIL("Error in: "__FILE__);
    }
    // Store grid and task id for later grid opening
    if (Grid->GridData->ReturnProcessorNumber() == MyProcessorNumber){
      OriginalGridID[Grid] = GridID;
#ifdef USE_HDF5_GROUPS
      OriginalTaskID[Grid] = Task;
#endif
    }
  }
 
  /* Read RandomForcingFields for the grid(s) on level 0. //AK */
 
  if (RandomForcing && ParentGrid == NULL && extract != TRUE 
      && LoadGridDataAtStart )
    if (Grid->GridData->ReadRandomForcingFields(fptr) == FAIL) {
      fprintf(stderr, "Error in grid->ReadRandomForcingFields (grid %"ISYM").\n",
              GridID);
      ENZO_FAIL("Error in: "__FILE__);
    }
 
  /* Read pointer information for the next grid this level. */
 
  if (fscanf(fptr, "Pointer: Grid[%"ISYM"]->NextGridThisLevel = %"ISYM"\n",
	     &TestGridID, &NextGridThisLevelID) != 2) {
    fprintf(stderr, "Error reading NextGridThisLevel pointer for grid %"ISYM".\n",
	    GridID);
    ENZO_FAIL("Error in: "__FILE__);
  }
  if (TestGridID != GridID) {
    fprintf(stderr, "GridID = %"ISYM" does not match grid(1) %"ISYM".\n",
	    TestGridID, GridID);
    ENZO_FAIL("Error in: "__FILE__);
  }
 
  /* If the pointer was non-zero, then read that grid. */
 
  if (NextGridThisLevelID != 0) {
    Grid->NextGridThisLevel = new HierarchyEntry;
    if (Group_ReadDataHierarchy(fptr, Grid->NextGridThisLevel, NextGridThisLevelID,
			  ParentGrid, file_id) == FAIL) {
      fprintf(stderr, "Error in Group_ReadDataHierarchy(1).\n");
      ENZO_FAIL("Error in: "__FILE__);
    }
  }
 
  /* Read pointer information for the next grid next level. */
 
  if (fscanf(fptr, "Pointer: Grid[%"ISYM"]->NextGridNextLevel = %"ISYM"\n",
	     &TestGridID, &NextGridNextLevelID) != 2) {
    fprintf(stderr, "Error reading NextGridNextLevel pointer for grid %"ISYM".\n",
	    GridID);
    ENZO_FAIL("Error in: "__FILE__);
  }
  if (TestGridID != GridID) {
    fprintf(stderr, "GridID = %"ISYM" does not match grid(2) %"ISYM".\n",
	    TestGridID, GridID);
    ENZO_FAIL("Error in: "__FILE__);
  }
 
  /* If the pointer was non-zero, then read that grid. */
 
  if (NextGridNextLevelID != 0) {
    Grid->NextGridNextLevel = new HierarchyEntry;
    if (Group_ReadDataHierarchy(fptr, Grid->NextGridNextLevel, NextGridNextLevelID, Grid, file_id)
	== FAIL) {
      fprintf(stderr, "Error in Group_ReadDataHierarchy(2).\n");
      ENZO_FAIL("Error in: "__FILE__);
    }
  }
 
  return SUCCESS;
}
