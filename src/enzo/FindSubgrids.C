/***********************************************************************
/
/  FIND SUBGRIDS FUNCTION
/
/  written by: Greg Bryan
/  date:       April, 1996
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
 
int IdentifyNewSubgridsBySignature(ProtoSubgrid *SubgridList[],
				   int &NumberOfSubgrids);
 
static ProtoSubgrid *SubgridList[MAX_NUMBER_OF_SUBGRIDS];
 
 
int FindSubgrids(HierarchyEntry *Grid, int level, int &TotalFlaggedCells,
		 int &FlaggedGrids)
{
 
  /* declarations */
#ifdef MPI_INSTRUMENTATION
  int GridMemory,NumberOfCells,CellsTotal,Particles;
  float AxialRatio, GridVolume;
#endif /* MPI_INSTRUMENTATION */
  int NumberOfFlaggedCells = INT_UNDEFINED, i;
  grid *CurrentGrid = Grid->GridData;
 
  /* Clear pointer to lower grids. */
 
  Grid->NextGridNextLevel = NULL;
 
  /* If this is the lowest allowed level, then return. */
 
  if (level >= MaximumRefinementLevel)
    return SUCCESS;
 
  /* If this grid is not on this processor, then return. */
 
  if (MyProcessorNumber != CurrentGrid->ReturnProcessorNumber())
    return SUCCESS;
 
  /* Clear the flagging field. */
 
  CurrentGrid->ClearFlaggingField();
 
  /* Set the flagging field. */

  if (CurrentGrid->SetFlaggingField(NumberOfFlaggedCells, level) == FAIL) {
    ENZO_FAIL("Error in grid->SetFlaggingField.");
  }
 
  /* Add a buffer region around each flagged cell. */
 
  if (NumberOfFlaggedCells != 0) {

    /* check flagged cells are in allowed refined region */
    
    if (CurrentGrid->SetFlaggingFieldMultiRefineRegions(level) 
	== FAIL) {
      fprintf(stderr, "Error in grid->SetFlaggingFieldMultiRefineRegions.\n");
      return FAIL;
    }

    NumberOfFlaggedCells = CurrentGrid->FlagBufferZones();

  } 

  /* Set the static (permanent) regions. */
  if (MustRefineParticlesCreateParticles != 3 &&
      MustRefineParticlesCreateParticles != 4) {
     if (CurrentGrid->SetFlaggingFieldStaticRegions(level, NumberOfFlaggedCells)
	== FAIL) {
      ENZO_FAIL("Error in grid->SetFlaggingFieldStaticRegions.");
    }

  }
  TotalFlaggedCells += NumberOfFlaggedCells;
  if (NumberOfFlaggedCells > 0)
    FlaggedGrids++;
  if (debug1)
    printf("RebuildHierarchy[%"ISYM"]: NumberOfFlaggedCells = %"ISYM".\n",
	   level, NumberOfFlaggedCells);
 
#ifdef MPI_INSTRUMENTATION
  Grid->GridData->CollectGridInformation
    (GridMemory,GridVolume, NumberOfCells,AxialRatio,CellsTotal, Particles);
  flagging_count ++;
  flagging_pct += float(NumberOfFlaggedCells) / NumberOfCells;
#endif
 
  if (NumberOfFlaggedCells != 0) {
 
    /* Create the base ProtoSubgrid which contains the whole grid. */
 
    int NumberOfSubgrids = 1;
    SubgridList[0] = new ProtoSubgrid;
    
    SubgridList[0]->SetLevel(level+1);
 
    /* Copy the flagged zones into the ProtoSubgrid. */
 
    if (SubgridList[0]->CopyFlaggedZonesFromGrid(CurrentGrid) == FAIL) {
      ENZO_FAIL("Error in ProtoSubgrid->CopyFlaggedZonesFromGrid.");
    }
 
    /* Recursively break up this ProtoSubgrid and add new ones based on the
       flagged cells. */
 
    if (IdentifyNewSubgridsBySignature(SubgridList, NumberOfSubgrids) == FAIL){
      ENZO_FAIL("Error in IdentifyNewSubgridsBySignature.");
    }
 
    /* For each subgrid, create a new grid based on the current grid (i.e.
       same parameters, etc.) */
 
    HierarchyEntry *PreviousGrid = Grid, *ThisGrid;

    if ( NumberOfSubgrids > MAX_NUMBER_OF_SUBGRIDS ) {
      ENZO_VFAIL("PE %"ISYM" NumberOfSubgrids > MAX_NUMBER_OF_SUBGRIDS\n", MyProcessorNumber)
    }
 
    for (i = 0; i < NumberOfSubgrids; i++) {
 
      /* create hierarchy entry */
 
      ThisGrid = new HierarchyEntry;
 
      /* set hierarchy values */
 
      if (PreviousGrid == Grid)
	Grid->NextGridNextLevel = ThisGrid;
      else
	PreviousGrid->NextGridThisLevel = ThisGrid;
      ThisGrid->NextGridNextLevel = NULL;
      ThisGrid->NextGridThisLevel = NULL;
      ThisGrid->ParentGrid        = Grid;
 
      /* create new grid */
 
      ThisGrid->GridData = new grid;
 
      /* set some the new grid's properties (rank, field types, etc.)
	 based on the current grid */
 
      ThisGrid->GridData->InheritProperties(Grid->GridData);
 
      /* Set the new grid's positional parameters.
         (The zero indicates there are no particles (for now). */
 
      ThisGrid->GridData->PrepareGrid(SubgridList[i]->ReturnGridRank(),
				      SubgridList[i]->ReturnGridDimension(),
				      SubgridList[i]->ReturnGridLeftEdge(),
				      SubgridList[i]->ReturnGridRightEdge(),
				      0);
 
      ThisGrid->GridData->SetProcessorNumber(MyProcessorNumber);
 
#ifdef UNUSED
      /* Create a new LevelHierarchyEntry and fill it in. */
 
      LevelHierarchyEntry *Temp2 = new LevelHierarchyEntry;
      if (*ListOfNewGrids != NULL)

	Temp2->NextGridThisLevel = (*ListOfNewGrids)->NextGridThisLevel;
      else
	Temp2->NextGridThisLevel = NULL;
      Temp2->GridData           = ThisGrid->GridData;
      Temp2->GridHierarchyEntry = ThisGrid;
      *ListOfNewGrids = Temp2;
#endif /* UNUSED */
 
      /* Go on to the next subgrid */
 
      PreviousGrid = ThisGrid;
      delete SubgridList[i];
 
    } // next subgrid
 
  }
 
  /* De-allocate the flagging field. */
 
  CurrentGrid->DeleteFlaggingField();
 
  /* done for this grid */
 
  return SUCCESS;
 
}
