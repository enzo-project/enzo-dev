//
// SetUpSiblingList
// Uses the FastSiblingLocator to set up the list of siblings of grids.
//

 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "performance.h"
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

int FastSiblingLocatorInitializeStaticChainingMesh(ChainingMeshStructure *Mesh, int Rank,
						   int TopGridDims[]); 
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
static int StaticSiblingListInitialized = 0; 
#ifdef STATIC_SIBLING_LIST
static SiblingGridList StaticSiblingList[MAX_NUMBER_OF_SUBGRIDS];
#endif
//SetUpSiblingList(Grids, NumberOfGrids,
int CreateSiblingList(HierarchyEntry ** Grids, int NumberOfGrids, SiblingGridList *SiblingList, 
		     int StaticLevelZero,TopGridData * MetaData,int level){

  int grid1;

  if (NumberOfGrids == 0)
    return SUCCESS;

  if (NumberOfGrids == 1) {   // only "sibling" is itself.
    SiblingList[0].NumberOfSiblings = 1;
    SiblingList[0].GridList = new grid*;
    SiblingList[0].GridList[0] = Grids[0]->GridData;
    return SUCCESS;
  }

#ifdef STATIC_SIBLING_LIST
  if ( StaticLevelZero == 1 && level == 0 ) {
    
    if (!StaticSiblingListInitialized) {
      
      if (debug) fprintf(stderr, "INITIALIZE Level 0 StaticSiblingList\n");
      
      ChainingMeshStructure StaticChainingMesh;
      
      FastSiblingLocatorInitializeStaticChainingMesh
	(&StaticChainingMesh, MetaData->TopGridRank, MetaData->TopGridDims);
      
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
        Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&StaticChainingMesh);
      
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
        if (Grids[grid1]->GridData->FastSiblingLocatorFindSiblings(
                              &StaticChainingMesh, &StaticSiblingList[grid1],
                              MetaData->LeftFaceBoundaryCondition,
                              MetaData->RightFaceBoundaryCondition) == FAIL) {
          ENZO_FAIL("Error in grid->FastSiblingLocatorFindSiblings.\n");
        }

      /* Clean up the chaining mesh. */

      FastSiblingLocatorFinalize(&StaticChainingMesh);

      StaticSiblingListInitialized = 1;

    }

  } // if StaticLevelZero && level == 0
#endif

  ChainingMeshStructure ChainingMesh;

#ifdef STATIC_SIBLING_LIST
  if (StaticLevelZero == 1 && level == 0 ) {

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      SiblingList[grid1].NumberOfSiblings = StaticSiblingList[grid1].NumberOfSiblings;
      SiblingList[grid1].GridList = StaticSiblingList[grid1].GridList;
    }

  }
#endif

  if (( StaticLevelZero == 1 && level != 0 ) || StaticLevelZero == 0 ) {

  FastSiblingLocatorInitialize(&ChainingMesh, MetaData->TopGridRank,
			       MetaData->TopGridDims);
 
  /* Add all the grids to the chaining mesh. */

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);

  /* For each grid, get a list of possible siblings from the chaining mesh. */
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (Grids[grid1]->GridData->FastSiblingLocatorFindSiblings(
                              &ChainingMesh, &SiblingList[grid1],
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition) == FAIL) {
      ENZO_FAIL("Error in grid->FastSiblingLocatorFindSiblings.\n");

    }
 
  /* Clean up the chaining mesh. */
 
  FastSiblingLocatorFinalize(&ChainingMesh);

  }

  return SUCCESS;
}
