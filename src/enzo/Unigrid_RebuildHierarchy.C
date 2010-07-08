/***********************************************************************
/
/  REBUILD HIERARCHY FUNCTION
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  August, 1995 by GB
/              Rewritten to rebuild an entire level (and below) at a time.
/  modified2:  February 2004, by Alexei Kritsuk; Added RandomForcing support.
/  modified3:  Robert Harkness
/  date:       March, 2008
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif
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
 
void AddLevel(LevelHierarchyEntry *LevelArray[], HierarchyEntry *Grid,
	      int level);
int FindSubgrids(HierarchyEntry *Grid, int level);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int  ReportMemoryUsage(char *header = NULL);
int CommunicationShareGrids(HierarchyEntry *GridHierarchyPointer[], int grids);
int CommunicationLoadBalanceGrids(HierarchyEntry *GridHierarchyPointer[],
				  int NumberOfGrids, int MoveParticles = TRUE);
int CommunicationTransferParticles(grid *GridPointer[], int NumberOfGrids,
				   int ShareParticles = TRUE);
 
 
/* RebuildHierarchy function */
 
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level)
{
 
  /* declarations */
 
  if (debug) printf("RebuildHierarchy: level = %"ISYM"\n", level);
  ReportMemoryUsage("Rebuild pos 1");
 
  int i, j, k, grids, grids2, subgrids;
  FLOAT ZeroVector[MAX_DIMENSION];
  LevelHierarchyEntry *Temp;
  for (i = 0; i < MAX_DIMENSION; i++)
    ZeroVector[i] = 0;
 
#ifdef MPI_INSTRUMENTATION
  double tmptime =  starttime;
  starttime = MPI_Wtime();
  timer[0] += starttime - tmptime;
  counter[0] ++;
#endif /* MPI_INSTRUMENTATION */
 
  /* --------------------------------------------------------------------- */
  /* For each grid on this level collect all the particles below it.
     Notice that this must be done even for static hierarchy's.  */
 
  HierarchyEntry *GridParent[MAX_NUMBER_OF_SUBGRIDS];
  grid           *GridPointer[MAX_NUMBER_OF_SUBGRIDS];
  grid           *ContigiousGridList[MAX_NUMBER_OF_SUBGRIDS];
  for (i = MAX_DEPTH_OF_HIERARCHY-1; i > level; i--) {
 
    Temp = LevelArray[i];
 
    /* Find the parents (at level=level) of all the grids (at level=i). */
 
    grids = 0;
    while (Temp != NULL) {
      GridPointer[grids] = Temp->GridData;
      GridParent[grids] = Temp->GridHierarchyEntry->ParentGrid;
      for (j = i-1; j > level; j--)
	GridParent[grids] = GridParent[grids]->ParentGrid;
      Temp = Temp->NextGridThisLevel;
      grids++;
    }
 
    /* Collect all the grids with the same parent and pass them all to
       MoveAllParticles (marking which ones have already been passed). */
 
    for (j = 0; j < grids; j++)
      if (GridPointer[j] != NULL) {
	grids2 = 0;
	for (k = j; k < grids; k++)
	  if (GridParent[k] == GridParent[j]) {
	    ContigiousGridList[grids2++] = GridPointer[k];
	    GridPointer[k] = NULL;
	  }
	if (GridParent[j]->GridData->MoveAllParticles(grids2,
						 ContigiousGridList) == FAIL) {
	  ENZO_FAIL("Error in grid->MoveAllParticles.\n");
	}
      }
 
  } // end: loop over levels
 
  /* --------------------------------------------------------------------- */
  /* if this is level 0 then transfer particles between grids. */
 
  ReportMemoryUsage("Rebuild pos 2");
  if (level == 0) {
    grids = 0;
    Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->DebugCheck("Before TransferParticles");
      GridPointer[grids++] = Temp->GridData;
      Temp = Temp->NextGridThisLevel;
    }
    if (CommunicationTransferParticles(GridPointer, grids) == FAIL) {
      ENZO_FAIL("Error in CommunicationTransferParticles.\n");
    }
  }
 
  /* --------------------------------------------------------------------- */
  /* For dynamic hierarchies, rebuild the grid structure. */
 
  ReportMemoryUsage("Rebuild pos 3");
  if (MetaData->StaticHierarchy == FALSE) {
 
//    if (debug) ReportMemoryUsage("Memory usage report: Rebuild 1");
 
    /* 1) Create a new TempLevelArray in which to keep the old grids. */
 
    LevelHierarchyEntry* TempLevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (i = level+1; i < MAX_DEPTH_OF_HIERARCHY; i++) {
      TempLevelArray[i] = LevelArray[i];
      LevelArray[i]     = NULL;
    }
    TempLevelArray[level] = LevelArray[level];
 
    /* 2) Clean up (delete excess baggage) all grids on this level and below.
          And delete the old hierarchy entries at the same time. */
 
    for (i = level; i < MAX_DEPTH_OF_HIERARCHY; i++) {
      Temp = TempLevelArray[i];
 
      while (Temp != NULL) {
	Temp->GridData->CleanUp();
	if (i > level)
	  delete Temp->GridHierarchyEntry;
	Temp = Temp->NextGridThisLevel;
      } // end: if (i > level)
 
    } // end: loop over levels
 
//    if (debug) ReportMemoryUsage("Memory usage report: Rebuild 3");
 
    /* 3) Rebuild all grids on this level and below.  Note: All the grids
          in LevelArray[level+] have been deleted. */
 
    for (i = level; i < MAX_DEPTH_OF_HIERARCHY-1; i++) {
 
      /* If there are no grids on this level, exit. */
 
      if (LevelArray[i] == NULL)
	break;
 
      /* 3a) Generate an array of grids on this level. */
 
      HierarchyEntry *GridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
      grids = 0;
      Temp = LevelArray[i];
      while (Temp != NULL) {
	GridHierarchyPointer[grids++] = Temp->GridHierarchyEntry;
	Temp                          = Temp->NextGridThisLevel;
      }
 
      /* 3b) Loop over grids creating new (but empty!) subgrids
	 (This also properly fills out the GridHierarchy tree). */
 
      for (j = 0; j < grids; j++)
	if (FindSubgrids(GridHierarchyPointer[j], i) == FAIL) {
	  ENZO_FAIL("Error in FindSubgrids.\n");
	}
 
      /* Create a temporary array of the new subgrids (which are on this
	 processor) for the next step. */
 
      HierarchyEntry *SubgridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS], *Temp2;
      subgrids = 0;
      for (j = 0; j < grids; j++) {
	Temp2 = GridHierarchyPointer[j]->NextGridNextLevel;
	while (Temp2 != NULL) {
	  SubgridHierarchyPointer[subgrids++] = Temp2;
	  Temp2 = Temp2->NextGridThisLevel;
	}
      }
 
      /* 3g) loop over parent, and copy particles to new grids
	     (all local to this processor) . */
 
      grid *ToGrids[MAX_NUMBER_OF_SUBGRIDS];
      for (j = 0; j < grids; j++)
 
	if (GridHierarchyPointer[j]->NextGridNextLevel != NULL) {
 
	  // if (debug) printf("grid %"ISYM":", j);
 
	  GridHierarchyPointer[j]->GridData->ZeroSolutionUnderSubgrid(
			   NULL, ZERO_UNDER_SUBGRID_FIELD);
 
	  for (k = 0; k < subgrids; k++) {
	    if (GridHierarchyPointer[j]->GridData->ZeroSolutionUnderSubgrid(
		              SubgridHierarchyPointer[k]->GridData,
		              ZERO_UNDER_SUBGRID_FIELD, float(k+1)) == FAIL) {
	      ENZO_FAIL("Error in grid->ZeroSolutionUnderSubgrid.\n");
	    }
	    ToGrids[k] = SubgridHierarchyPointer[k]->GridData;
	  }
 
	  if (GridHierarchyPointer[j]->GridData->MoveSubgridParticlesFast(
				 subgrids, ToGrids, TRUE) == FAIL) {
	    ENZO_FAIL("Error in grid->MoveSubgridParticlesFast.\n");
	  }
 
	}
 
      /* Share the new grids amoung processors. */
 
      CommunicationShareGrids(GridHierarchyPointer, grids);
 
      /* 3c) Combine the many linked-lists of subgrids into the LevelArray
	 linked list. */
 
      for (j = 0; j < grids; j++)
	if (GridHierarchyPointer[j]->NextGridNextLevel != NULL)
	  AddLevel(LevelArray, GridHierarchyPointer[j]->NextGridNextLevel,i+1);
 
      /* 3d) Create an array of the new subgrids. */
 
      subgrids = 0;
      Temp = LevelArray[i+1];
      while (Temp != NULL) {
	SubgridHierarchyPointer[subgrids++] = Temp->GridHierarchyEntry;
	Temp                                = Temp->NextGridThisLevel;
      }
 
      /* 3e) Loop over the new subgrids and record in the old subgrids how
	 many times they are needed (number of overlaps with new subgrids). */
 
      int Overlap, oldgrid = 0, NumberOfOverlaps[MAX_NUMBER_OF_SUBGRIDS];
      Temp = TempLevelArray[i+1];
      while (Temp != NULL) {
 
	NumberOfOverlaps[oldgrid] = 0;
	for (j = 0; j < subgrids; j++) {
	  if (SubgridHierarchyPointer[j]->GridData->CopyZonesFromGridCountOnly(
		                          Temp->GridData, Overlap) == FAIL) {
	    ENZO_FAIL("Error in grid->CopyZonesFromGridCountOnly.\n");
	  }
	  NumberOfOverlaps[oldgrid] += Overlap;
	}
 
	if (NumberOfOverlaps[oldgrid] == 0) {
	  delete Temp->GridData;
	  Temp->GridData = NULL;
	}
 
	Temp = Temp->NextGridThisLevel;
	oldgrid++;
      }
 
      /* 3f) For each new subgrid, interpolate from parent and then copy
	 from old subgrids.  Also, copy particles (if present).  For each
	 old subgrid, decrement the Overlap counter, deleting the grid
	 which it reaches zero. */
 
      for (j = 0; j < subgrids; j++) {
	SubgridHierarchyPointer[j]->ParentGrid->GridData->DebugCheck(
						        "Rebuild parent");
        if (RandomForcing) { //AK
          SubgridHierarchyPointer[j]->GridData->AppendForcingToBaryonFields();
          SubgridHierarchyPointer[j]->ParentGrid->GridData->AppendForcingToBaryonFields();
        }
	SubgridHierarchyPointer[j]->GridData->InterpolateFieldValues(
		       SubgridHierarchyPointer[j]->ParentGrid->GridData);
        if (RandomForcing) { //AK
          SubgridHierarchyPointer[j]->GridData->RemoveForcingFromBaryonFields();
          SubgridHierarchyPointer[j]->ParentGrid->GridData->RemoveForcingFromBaryonFields();
        }
	SubgridHierarchyPointer[j]->GridData->DebugCheck("Rebuild child");
      }
 
      /* Copy from old grids. */
 
      for (j = 0; j < subgrids; j++) {
 
	oldgrid = 0;
	Temp = TempLevelArray[i+1];
	while (Temp != NULL) {
 
	  if (Temp->GridData != NULL) {
 
	    /* Copy from old subgrid. */
 
	    if (SubgridHierarchyPointer[j]->GridData->CopyZonesFromGrid(
                                       Temp->GridData, ZeroVector) == FAIL) {
	      ENZO_FAIL("Error in grid->CopyZonesFromGrid.\n");
	    }
 
	    /* Check if we can delete the old subgrid yet. */
 
	    SubgridHierarchyPointer[j]->GridData->CopyZonesFromGridCountOnly(
		                                  Temp->GridData, Overlap);
 
	    if (Overlap == TRUE)
	      if (--NumberOfOverlaps[oldgrid] <= 0) {
		delete Temp->GridData;
		Temp->GridData = NULL;
	      }
 
	  } // end: if (Temp->GridData != NULL)
 
	  /* Next old subgrid. */
 
	  Temp = Temp->NextGridThisLevel;
	  oldgrid++;
 
	} // end: loop over old subgrids
 
//	if (debug) ReportMemoryUsage("Memory usage report: PostBuild");
 
      /* Extraneous time step calculation; not needed here;
	 see emails around Mar 24th, 2004; Contact Dave Collins. //AK
 
	 SubgridHierarchyPointer[j]->GridData->ComputeTimeStep();*/
 
      } // end: loop over new subgrids
 
      /* Redistribute grids over processors to Load balance. */

#ifdef ENABLE_LOAD_BALANCE 
      CommunicationLoadBalanceGrids(SubgridHierarchyPointer, subgrids);
#endif
 
      /* 3h) Clean up the LevelHierarchy entries for the old subgrids.
	     Also, we can check to see if any old subgrids were missed. */
 
      while (TempLevelArray[i+1] != NULL) {
	Temp = TempLevelArray[i+1]->NextGridThisLevel;
 
	if (TempLevelArray[i+1]->GridData != NULL) {
	  ENZO_FAIL("An old subgrid was not deleted.  Why?\n");
	}
 
	/* Remove the LevelHierarchy entry for that grid. */
 
	delete TempLevelArray[i+1];
	TempLevelArray[i+1] = Temp;
      }
 
    } // end: loop over levels
 
  } // end: if (StaticHierarchy == FALSE)
 
  /* --------------------------------------------------------------------- */
  /* Redistribute particles: for each grid, move the particles that belong in
     it's subgrids (this only has to be done if we didn't just do a rebuild
      since the rebuild does this as it goes). */
 
  if (MetaData->StaticHierarchy == TRUE) {
    for (i = level; i < MAX_DEPTH_OF_HIERARCHY-1; i++) {
 
      /* If there are no grids on this level, exit. */
 
      if (LevelArray[i] == NULL)
	break;
 
      /* 3a) Generate an array of grids on this level. */
 
      HierarchyEntry *GridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
      grids = 0;
      Temp = LevelArray[i];
      while (Temp != NULL) {
	GridHierarchyPointer[grids++] = Temp->GridHierarchyEntry;
	Temp                          = Temp->NextGridThisLevel;
      }
 
      /* 3d) Create an array of the subgrids. */
 
      HierarchyEntry *SubgridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
      subgrids = 0;
      Temp = LevelArray[i+1];
      while (Temp != NULL) {
	SubgridHierarchyPointer[subgrids++] = Temp->GridHierarchyEntry;
	Temp                                = Temp->NextGridThisLevel;
      }
 
      /* 3g) loop over parent, and copy particles to new grids. */
 
      grid *ToGrids[MAX_NUMBER_OF_SUBGRIDS];
      for (j = 0; j < grids; j++)
 
	if (GridHierarchyPointer[j]->NextGridNextLevel != NULL) {
 
	  GridHierarchyPointer[j]->GridData->ZeroSolutionUnderSubgrid(
			   NULL, ZERO_UNDER_SUBGRID_FIELD);
 
	  for (k = 0; k < subgrids; k++) {
	    if (GridHierarchyPointer[j]->GridData->ZeroSolutionUnderSubgrid(
		              SubgridHierarchyPointer[k]->GridData,
		              ZERO_UNDER_SUBGRID_FIELD, float(k+1)) == FAIL) {
	      ENZO_FAIL("Error in grid->ZeroSolutionUnderSubgrid.\n");
	    }
	    ToGrids[k] = SubgridHierarchyPointer[k]->GridData;
	  }
 
	  if (GridHierarchyPointer[j]->GridData->MoveSubgridParticlesFast(
				 subgrids, ToGrids, FALSE) == FAIL) {
	    ENZO_FAIL("Error in grid->MoveSubgridParticlesFast.\n");
	  }
 
	}
 
      /* Set boundary conditions. */
 
      LevelHierarchyEntry *Temp = LevelArray[i+1];
      while (Temp != NULL) {
 
	if (Temp->GridData->InterpolateBoundaryFromParent
	    (Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL) {
	  ENZO_FAIL("Error in grid->InterpolateBoundaryFromParent.\n");
	}
 
	Temp = Temp->NextGridThisLevel;
      }
 
    } // end: loop over levels
 
  } // end: if (StaticHierarchy == TRUE)

 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[1] += endtime - starttime;
  counter[1] ++;
#endif /* MPI_INSTRUMENTATION */
 
  /* Done for this level. */
 
  ReportMemoryUsage("Rebuild pos 4");
  return SUCCESS;
 
}
