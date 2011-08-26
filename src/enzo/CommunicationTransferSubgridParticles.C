/***********************************************************************
/
/  COMMUNICATION ROUTINE: TRANSFER PARTICLES AND STARS TO SIBLING GRIDS
/
/  written by: John Wise
/  date:       September, 2009
/  modified:   
/
/  PURPOSE: Move particles to sibling grids if they have moved 
/           beyond the parent level-0 grid boundary in the last 
/           timestep.  CommunicationTransferParticles does not handle
/           this.
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
 
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
#include "CommunicationUtilities.h"
 
// function prototypes
 
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],
				       int NumberOfGrids);
int CommunicationShareParticles(int *NumberToMove, particle_data* &SendList,
				int &NumberOfReceives,
				particle_data* &SharedList);
int CommunicationShareStars(int *NumberToMove, star_data* &SendList,
			    int &NumberOfReceives, star_data* &SharedList);

int CommunicationTransferSubgridParticles(LevelHierarchyEntry *LevelArray[],
					  TopGridData *MetaData, int level)
{

  int proc, i, j, k, jstart, jend, TotalNumber;
  int particle_data_size, star_data_size;
  int Zero = 0;

  HierarchyEntry **Grids;
  grid** GridPointers;
  int grid1, grid2, NumberOfGrids, ID;
  ChainingMeshStructure ChainingMesh;
  SiblingGridList SiblingList;

  /* Star and particle lists for communication */

  particle_data *SendList = NULL;
  particle_data *SharedList = NULL;
  star_data *StarSendList = NULL;
  star_data *StarSharedList = NULL;

  int NumberOfReceives, StarNumberOfReceives;
  int *NumberToMove = new int[NumberOfProcessors];
  int *StarsToMove = new int[NumberOfProcessors];

  // Needed for TransferSubgridParticles.  This means that the
  // particles can be transferred to grids on other processors.
  bool KeepLocal = false;
  bool ParticlesAreLocal = true;

  /******************************************************************/

  for (i = 0; i < NumberOfProcessors; i++) {
    NumberToMove[i] = 0;
    StarsToMove[i] = 0;
  }

  /* Generate array of all grids on this level */

  NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);

  /* TransferSubgridParticles needs grid pointers instead of
     HierarchyEntry pointers. */

  GridPointers = new grid*[NumberOfGrids];
  for (i = 0; i < NumberOfGrids; i++)
    GridPointers[i] = Grids[i]->GridData;

  /* Initialize and fill out the fast sibling chaining mesh */

  FastSiblingLocatorInitialize(&ChainingMesh, MetaData->TopGridRank,
			       MetaData->TopGridDims);
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);

  /* Loop over grids and count particles that need transferring to
     sibling grids.  In the next step, we will allocate the memory and
     transfer them. */

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

    /* Get a list of possible siblings from the chaining mesh */

    Grids[grid1]->GridData->FastSiblingLocatorFindSiblings
      (&ChainingMesh, &SiblingList, MetaData->LeftFaceBoundaryCondition,
       MetaData->RightFaceBoundaryCondition);
    
    /* Fill out the UnderSubgrid field with the index of the grid that
       we want to transfer the particles to (excluding self and
       including ghost zones). */

    Grids[grid1]->GridData->ZeroSolutionUnderSubgrid
      (NULL, ZERO_UNDER_SUBGRID_FIELD);

    /* Last two arguments (AllProcessors = FALSE, IncludeGhostZones =
       TRUE) mean that each grid only marked on its host processor and
       mark cells in the ghost zones, respectively. */

    for (grid2 = 0; grid2 < SiblingList.NumberOfSiblings; grid2++) {
      ID = SiblingList.GridList[grid2]->GetGridID();
      if (Grids[grid1]->GridData != SiblingList.GridList[grid2])
	Grids[grid1]->GridData->ZeroSolutionUnderSubgrid
	  (SiblingList.GridList[grid2], ZERO_UNDER_SUBGRID_FIELD,
	   float(ID+1), FALSE, TRUE);
    }

    /* 
       Now collect any particles that have crossed the grid boundary.
       Re-use grid::TransferSubgridParticles(), but give it sibling
       grids instead of subgrids. The last argument (TRUE) makes it
       consider ghost zones.

       We use the whole grid list instead of the sibling list because
       we marked the subgrid field with grid IDs.  This is needed when
       we move the particles back into the grids because the sharing
       routine (below) assumes a complete grid list.  Note that this
       doesn't slow it down because we're only using the grid list for
       grid processor numbers.
    */

    Grids[grid1]->GridData->TransferSubgridStars
      (GridPointers, NumberOfGrids, StarsToMove, Zero, Zero, 
       StarSendList, KeepLocal, ParticlesAreLocal, COPY_OUT, TRUE);

    Grids[grid1]->GridData->TransferSubgridParticles
      (GridPointers, NumberOfGrids, NumberToMove, Zero, Zero, 
       SendList, KeepLocal, ParticlesAreLocal, COPY_OUT, TRUE, TRUE);

    delete [] SiblingList.GridList;

  } // ENDFOR grid1

  /* Allocate the memory for the move list and transfer the particles */

  TotalNumber = 0;
  for (j = 0; j < NumberOfProcessors; j++) {
    TotalNumber += NumberToMove[j];
    NumberToMove[j] = 0;  // Zero-out to use in the next step
  }
  SendList = new particle_data[TotalNumber];

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

    Grids[grid1]->GridData->TransferSubgridParticles
      (GridPointers, NumberOfGrids, NumberToMove, Zero, Zero, 
       SendList, KeepLocal, ParticlesAreLocal, COPY_OUT, TRUE, FALSE);

  } // ENDFOR grid1

  FastSiblingLocatorFinalize(&ChainingMesh);

  /* Now we have a list of particles to move to subgrids, so we
     communicate them with all processors. */

  CommunicationShareParticles(NumberToMove, SendList, NumberOfReceives,
			      SharedList);
  CommunicationShareStars(StarsToMove, StarSendList, StarNumberOfReceives,
			  StarSharedList);

  int TotalNumberOfReceives = NumberOfReceives;
  int TotalStarNumberOfReceives = StarNumberOfReceives;
  CommunicationSumValues(&TotalNumberOfReceives, 1);
  CommunicationSumValues(&TotalStarNumberOfReceives, 1);

  if (debug)
    printf("TransferSubgridParticles[%"ISYM"]: Moved %"ISYM" particles, "
	   "%"ISYM" stars.\n", level, TotalNumberOfReceives,
	   TotalStarNumberOfReceives);
	   

  /*******************************************************************/
  /****************** Copy particles back to grids. ******************/
  /*******************************************************************/

  jstart = 0;
  jend = 0;

  // Copy shared particles to grids, if any

  if (NumberOfReceives > 0)
    for (j = 0; j < NumberOfGrids && jend < NumberOfReceives; j++) {
      while (SharedList[jend].grid <= j) {
	jend++;
	if (jend == NumberOfReceives) break;
      }

      /*
      printf("j =%d, jstart =%d, jend =%d, NumberOfGrids =%d, " 
             "NumberToMove[] =%d/%d, NumberOfReceives =%d\n", 
	     j, jstart, jend, NumberOfGrids, 
	     NumberToMove[0], NumberToMove[1], NumberOfReceives); 
      */

      GridPointers[j]->TransferSubgridParticles
	(GridPointers, NumberOfGrids, NumberToMove, jstart, jend, 
	 SharedList, KeepLocal, ParticlesAreLocal, COPY_IN, TRUE);
      
      jstart = jend;
    } // ENDFOR grids

  /*******************************************************************/
  /******************** Copy stars back to grids. ********************/
  /*******************************************************************/

  jstart = 0;
  jend = 0;

  // Copy shared stars to grids, if any

  if (StarNumberOfReceives > 0)
    for (j = 0; j < NumberOfGrids && jend < StarNumberOfReceives; j++) {
      while (StarSharedList[jend].grid <= j) {
	jend++;
	if (jend == StarNumberOfReceives) break;
      }
      
      GridPointers[j]->TransferSubgridStars
	(GridPointers, NumberOfGrids, StarsToMove, jstart, jend, 
	 StarSharedList, KeepLocal, ParticlesAreLocal, COPY_IN, TRUE);
      
      jstart = jend;
    } // ENDFOR grids

  /************************************************************************
     Since the particles and stars are only on the grid's host
     processor, set number of particles so everybody agrees.
  ************************************************************************/

  CommunicationSyncNumberOfParticles(Grids, NumberOfGrids);

  /* Cleanup. */

  if (SendList != SharedList)
    delete [] SendList;
  delete [] SharedList;
  if (StarSendList != StarSharedList)
    delete [] StarSendList;
  delete [] StarSharedList;

  delete [] Grids;
  delete [] GridPointers;

  delete [] NumberToMove;
  delete [] StarsToMove;

  return SUCCESS;
}

