/***********************************************************************
/
/  READ IN THE DATA HIERARCHY (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Alexei Kritsuk, Jan 2004 now reads RandomForcing fields //AK
/  modified2:  Robert harkness, Jan 2006
/              Read Unigrid Grid-to-MPI task map
/              This is necessary only for ultra-large grids on
/              node memory challenged systems
/  modified3:  Michael Kuhlen, October 2010, HDF5 hierarchy
/
/  PURPOSE:
/
************************************************************************/
 
// This function reads in the data hierarchy (TopGrid)

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

/* function prototypes */
 
static int ReadDataGridCounter = 0;



 
int ReadDataHierarchy(FILE *fptr, hid_t Hfile_id, HierarchyEntry *Grid, int GridID, HierarchyEntry *ParentGrid, FILE *log_fptr)
{
 
  int TestGridID, NextGridThisLevelID, NextGridNextLevelID;
  int Task = 0;
 
  char DataFilename[MAX_LINE_LENGTH];

  //dcollins, August 5 2009.  Updated failsafe for old files that don't have Task defined.
  int NewProc = ReadDataGridCounter % NumberOfProcessors;
  int ProcMap = ABS(NewProc - NumberOfProcessors) % NumberOfProcessors;
  
  FILE * ptr_task_check = fptr;
  
  if (HierarchyFileInputFormat == 1) {

    /* Read header info for this grid */
    
    if (fscanf(fptr, "\nGrid = %"ISYM"\n", &TestGridID) != 1) {
      ENZO_VFAIL("Error reading Grid # in grid %"ISYM".\n", GridID)
	}
    if (TestGridID != GridID) {
      ENZO_VFAIL("Unexpected GridID = %"ISYM" while reading grid %"ISYM".\n",
		 TestGridID, GridID)
	}
        
    if( fscanf(fptr, "Task = %"ISYM"\n", &Task) != 1){
      Task = NewProc;
      fptr = ptr_task_check;
    }
  }    

  /* Create new grid and fill out hierarchy entry. */
 
  Grid->GridData          = new grid;
  Grid->NextGridThisLevel = NULL;
  Grid->NextGridNextLevel = NULL;
  Grid->ParentGrid        = ParentGrid;


  if (HierarchyFileInputFormat % 2 == 0) {
    TestGridID = GridID;
    
    Grid->GridData->ReadHierarchyInformationHDF5(Hfile_id, TestGridID, Task, NextGridThisLevelID, NextGridNextLevelID, DataFilename, log_fptr);
    
  }

  if ( MyProcessorNumber == 0 )
      fprintf(stderr, "Reading Grid %"ISYM" assigned to Task %"ISYM"\n", TestGridID, Task);
 

// If explicit task mapping is enabled (for Unigrid) then use that map

#ifdef ENABLE_TASKMAP

  if ( MyProcessorNumber == 0 )
    fprintf(stderr, "Task map assignment: GridID = %"ISYM"  MPI Task = %"ISYM"\n", GridID, TaskMap[GridID-1]);
  
  Grid->GridData->SetProcessorNumber(TaskMap[GridID-1]);
  
#else

// Use grid-to-processor mapping from last dump

  if (Task > NumberOfProcessors-1) {
    Task = NewProc;
    fptr = ptr_task_check;
  }

  if ( MyProcessorNumber == 0 )
    fprintf(stderr, "Using dumped task assignment: GridID = %"ISYM"  MPI Task = %"ISYM"\n", GridID, Task);

    Grid->GridData->SetProcessorNumber(Task);

#endif

    //dcollins: moved higher in the code.
    //int NewProc = ReadDataGridCounter % NumberOfProcessors;
    //int ProcMap = ABS(NewProc - NumberOfProcessors) % NumberOfProcessors;

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
    if (Grid->GridData->ReadGrid(fptr, GridID, DataFilename) == FAIL) {
      ENZO_VFAIL("Error in grid->ReadGrid (grid %"ISYM").\n", GridID)
    }
  }else{
    if (Grid->GridData->ReadGrid(fptr, GridID, DataFilename, TRUE, FALSE) == FAIL) {
      ENZO_VFAIL("Error in grid->ReadGrid (grid %"ISYM").\n", GridID)
    }
    // Store grid id for later grid opening
    if (Grid->GridData->ReturnProcessorNumber() == MyProcessorNumber)
      OriginalGridID[Grid] = GridID;
  }
 
  /* Read RandomForcingFields for the grid(s) on level 0. //AK */
 
  if (RandomForcing && ParentGrid == NULL && extract != TRUE
      && LoadGridDataAtStart )
    if (Grid->GridData->ReadRandomForcingFields(fptr, DataFilename) == FAIL) {
      ENZO_VFAIL("Error in grid->ReadRandomForcingFields (grid %"ISYM").\n",
              GridID)
    }
 
  if (HierarchyFileInputFormat == 1) {
    /* Read pointer information for the next grid this level. */
    
    if (fscanf(fptr, "Pointer: Grid[%"ISYM"]->NextGridThisLevel = %"ISYM"\n",
	       &TestGridID, &NextGridThisLevelID) != 2) {
      ENZO_VFAIL("Error reading NextGridThisLevel pointer for grid %"ISYM".\n",
		 GridID)
	}
    if (TestGridID != GridID) {
      ENZO_VFAIL("GridID = %"ISYM" does not match grid(1) %"ISYM".\n",
		 TestGridID, GridID)
	}
  }

  /* If the pointer was non-zero, then read that grid. */
 
  if (NextGridThisLevelID != 0) {
    Grid->NextGridThisLevel = new HierarchyEntry;
    if (ReadDataHierarchy(fptr, Hfile_id, Grid->NextGridThisLevel, NextGridThisLevelID, ParentGrid, log_fptr) == FAIL) {
      ENZO_FAIL("Error in ReadDataHierarchy(1).\n");
    }
  }
 
  if (HierarchyFileInputFormat == 1) {
    /* Read pointer information for the next grid next level. */
    
    if (fscanf(fptr, "Pointer: Grid[%"ISYM"]->NextGridNextLevel = %"ISYM"\n",
	       &TestGridID, &NextGridNextLevelID) != 2) {
      ENZO_VFAIL("Error reading NextGridNextLevel pointer for grid %"ISYM".\n",
		 GridID)
	}
    if (TestGridID != GridID) {
      ENZO_VFAIL("GridID = %"ISYM" does not match grid(2) %"ISYM".\n",
		 TestGridID, GridID)
	}
  }

  /* If the pointer was non-zero, then read that grid. */
 
  if (NextGridNextLevelID != 0) {
    Grid->NextGridNextLevel = new HierarchyEntry;
    if (ReadDataHierarchy(fptr, Hfile_id, Grid->NextGridNextLevel, NextGridNextLevelID,Grid, log_fptr)
	== FAIL) {
      ENZO_FAIL("Error in ReadDataHierarchy(2).\n");

    }
  }
 
  return SUCCESS;
}
