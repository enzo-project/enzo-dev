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
/  modified4:  Michael Kuhlen, October 2010, HDF5 hierarchy
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


int Group_ReadDataHierarchy(FILE *fptr, hid_t Hfile_id, HierarchyEntry *Grid,
                            TopGridData &MetaData, int GridID,
                            HierarchyEntry *ParentGrid, hid_t file_id,
                            int NumberOfRootGrids, int *RootGridProcessors,
                            bool ReadParticlesOnly, FILE *log_fptr)
{
 
  int TestGridID, NextGridThisLevelID, NextGridNextLevelID;
  int Task;
 
  char DataFilename[MAX_LINE_LENGTH];
  for(int i=0;i<MAX_LINE_LENGTH;i++) DataFilename[i] = 0;

  if (HierarchyFileInputFormat == 1) {

    /* Read header info for this grid */
    
    if (fscanf(fptr, "\nGrid = %"ISYM"\n", &TestGridID) != 1) {
      ENZO_VFAIL("Error reading Grid # in grid %"ISYM".\n", GridID)
	}
    if (TestGridID != GridID) {
      ENZO_VFAIL("Unexpected GridID = %"ISYM" while reading grid %"ISYM".\n",
		 TestGridID, GridID)
	}
    
    fscanf(fptr, "Task = %"ISYM"\n", &Task);
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
  
  Task = Task % NumberOfProcessors;
  //  if ( MyProcessorNumber == 0 )
  //    fprintf(stderr, "Reading Grid %"ISYM" assigned to Task %"ISYM"\n", TestGridID, Task);


#ifdef SINGLE_HDF5_OPEN_ON_INPUT

// Explicit task mapping cannot be used if SINGLE_HDF5_OPEN_ON_INPUT is in effect

#ifdef ENABLE_TASKMAP
  ENZO_FAIL("Task map cannot be used with SINGLE_HDF5_OPEN_ON_INPUT!\n");
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

  /* If requested, reset the grid processors. */

  if (ResetLoadBalancing) {
    if (GridID <= NumberOfRootGrids) {
      if (LoadBalancing == 2 && PreviousMaxTask < NumberOfProcessors-1) {
        Task = (GridID-1) * NumberOfProcessors / (PreviousMaxTask+1);
      } else {
        Task = (GridID-1) % NumberOfProcessors;
      } 
    } else {
      Task = ParentGrid->GridData->ReturnProcessorNumber();
    }
  }

  // Assign tasks in a round-robin fashion if we're increasing the
  // processor count, but processors are grouped together
  // (0000111122223333).  Only used if LoadBalancing == 2.
  if (LoadBalancing == 2 && PreviousMaxTask < NumberOfProcessors-1)
    Task = Task * NumberOfProcessors / (PreviousMaxTask+1);

  // If provided load balancing of root grids based on subgrids, use
  // these instead.
  if (LoadBalancing > 1 && RootGridProcessors != NULL) {
    if (GridID <= NumberOfRootGrids) {
      Task = RootGridProcessors[GridID-1];
    } else if (StaticRefineRegionLevel[0] == INT_UNDEFINED) {
      // Load the child on the same processor as its parent only if
      // it's not a zoom-in calculation
      Task = ParentGrid->GridData->ReturnProcessorNumber();
    }
  }
  
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

  /* We pass in the global here as a parameter to the function */
  if(LoadGridDataAtStart){ 
    if (Grid->GridData->Group_ReadGrid(fptr, GridID, file_id, DataFilename, 
				       TRUE, TRUE, ReadParticlesOnly
#ifdef NEW_GRID_IO
                , CheckpointRestart
#endif
        ) == FAIL) {
      ENZO_VFAIL("Error in grid->Group_ReadGrid (grid %"ISYM").\n", GridID)
    }
  }else{
    if (Grid->GridData->Group_ReadGrid(fptr, GridID, file_id, DataFilename, TRUE, FALSE) == FAIL) {
      ENZO_VFAIL("Error in grid->Group_ReadGrid (grid %"ISYM").\n", GridID)
    }
    // Store grid and task id for later grid opening
    if (Grid->GridData->ReturnProcessorNumber() == MyProcessorNumber){
      OriginalGridID[Grid] = GridID;
#ifdef USE_HDF5_GROUPS
      OriginalTaskID[Grid] = Task;
#endif
    }
  }
 
  /* Initialize grid hydrodynamics parameters */
  Grid->GridData->SetHydroParameters(MetaData.CourantSafetyNumber,
                                     MetaData.PPMFlatteningParameter,
                                     MetaData.PPMDiffusionParameter,
                                     MetaData.PPMSteepeningParameter);

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
    if (Group_ReadDataHierarchy(fptr, Hfile_id, Grid->NextGridThisLevel,
                                MetaData, NextGridThisLevelID, ParentGrid,
                                file_id, NumberOfRootGrids, RootGridProcessors,
                                ReadParticlesOnly, log_fptr) == FAIL)
      ENZO_FAIL("Error in Group_ReadDataHierarchy(1).");
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
    if (Group_ReadDataHierarchy(fptr, Hfile_id, Grid->NextGridNextLevel,
                                MetaData, NextGridNextLevelID, Grid, file_id,
                                NumberOfRootGrids, RootGridProcessors,
                                ReadParticlesOnly, log_fptr) == FAIL)
      ENZO_FAIL("Error in Group_ReadDataHierarchy(2).");
  }
 
  return SUCCESS;
}
