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
#include "ActiveParticle.h"

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
int CommunicationShareActiveParticles(int *NumberToMove,
        ActiveParticleList<ActiveParticleType> &SendList, int &NumberOfReceives,
        ActiveParticleList<ActiveParticleType> &SharedList);

int CommunicationTransferSubgridParticles(LevelHierarchyEntry *LevelArray[],
					  TopGridData *MetaData, int level)
{

  int proc, i, j, k, jstart, jend, TotalNumber, APTotalNumber;

  int particle_data_size, star_data_size;
  int Zero = 0;

  HierarchyEntry **Grids;
  grid** GridPointers;
  int grid1, grid2, NumberOfGrids, ID;
  ChainingMeshStructure ChainingMesh;
  SiblingGridList SiblingList;

  int *jstart = NULL;
  int *jend = NULL;
  int *kstart = NULL;
  int *kend = NULL;
  
  /* Star and particle lists for communication */

  particle_data *SendList = NULL;
  particle_data *SharedList = NULL;
  star_data *StarSendList = NULL;
  star_data *StarSharedList = NULL;
  ActiveParticleList<ActiveParticleType> APSendList;
  ActiveParticleList<ActiveParticleType> APSharedList;

  int NumberOfReceives=0, StarNumberOfReceives=0, APNumberOfReceives=0;
  int *NumberToMove = new int[NumberOfProcessors];
  int *StarsToMove = new int[NumberOfProcessors];
  int *APNumberToMove = new int[NumberOfProcessors];

  // Needed for TransferSubgridParticles.  This means that the
  // particles can be transferred to grids on other processors.
  bool KeepLocal = false;
  bool ParticlesAreLocal = true;

  /******************************************************************/

  for (i = 0; i < NumberOfProcessors; i++) {
    NumberToMove[i] = 0;
    StarsToMove[i] = 0;
    APNumberToMove[i] = 0;
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

  int ParticleCounter1 = 0;
  int ParticleCounter2 = 0;
#pragma omp parallel
{
#pragma omp for reduction(+:NumberToMove[:NumberOfProcessors], StarsToMove[:NumberOfProcessors]) private(ID, grid2, SiblingList)
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

    Grids[grid1]->GridData->TransferSubgridActiveParticles
      (GridPointers, NumberOfGrids, APNumberToMove, Zero, Zero,
       APSendList, KeepLocal, ParticlesAreLocal, COPY_OUT, TRUE, TRUE);

    Grids[grid1]->GridData->TransferSubgridParticles
      (GridPointers, NumberOfGrids, NumberToMove, ParticleCounter1, Zero, Zero, 
       SendList, KeepLocal, ParticlesAreLocal, COPY_OUT, TRUE, TRUE);

    delete [] SiblingList.GridList;

  } // ENDFOR grid1

  /* Allocate the memory for the move list and transfer the particles */

#pragma omp single
{
  TotalNumber = 0;
  APTotalNumber = 0;
  for (j = 0; j < NumberOfProcessors; j++) {
    TotalNumber += NumberToMove[j];
    APTotalNumber += APNumberToMove[j];
    NumberToMove[j] = 0;  // Zero-out to use in the next step
    APNumberToMove[j] = 0;
  }
  SendList = new particle_data[TotalNumber];
} // end omp single

#pragma omp for reduction(+:NumberToMove[:NumberOfProcessors])
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

    Grids[grid1]->GridData->TransferSubgridActiveParticles
      (GridPointers, NumberOfGrids, APNumberToMove, Zero, Zero,
       APSendList, KeepLocal, ParticlesAreLocal, COPY_OUT, TRUE, FALSE);

    Grids[grid1]->GridData->TransferSubgridParticles
      (GridPointers, NumberOfGrids, NumberToMove, ParticleCounter2, Zero, Zero, 
       SendList, KeepLocal, ParticlesAreLocal, COPY_OUT, TRUE, FALSE);

  } // ENDFOR grid1

} // end omp parallel section

  FastSiblingLocatorFinalize(&ChainingMesh);

  /* Now we have a list of particles to move to subgrids, so we
     communicate them with all processors. */

  CommunicationShareParticles(NumberToMove, SendList, NumberOfReceives,
			      SharedList);
  CommunicationShareStars(StarsToMove, StarSendList, StarNumberOfReceives,
			  StarSharedList);
  CommunicationShareActiveParticles(APNumberToMove, APSendList, APNumberOfReceives,
                                    APSharedList);

  int TotalNumberOfReceives = NumberOfReceives;
  int TotalStarNumberOfReceives = StarNumberOfReceives;
  int TotalAPNumberOfReceives = APNumberOfReceives;
  CommunicationSumValues(&TotalNumberOfReceives, 1);
  CommunicationSumValues(&TotalStarNumberOfReceives, 1);
  CommunicationSumValues(&TotalAPNumberOfReceives, 1);

  if (debug)
    printf("TransferSubgridParticles[%"ISYM"]: Moved %"ISYM" particles, "
	   "%"ISYM" stars, %"ISYM" active particles.\n", level, TotalNumberOfReceives,
           TotalStarNumberOfReceives, TotalAPNumberOfReceives);
	   

  /*******************************************************************/
  /****************** Copy particles back to grids. ******************/
  /*******************************************************************/

  jstart = new int[NumberOfGrids+1];
  jend = new int[NumberOfGrids+1];

  memset(jstart, 0, sizeof(int)*(NumberOfGrids+1));
  memset(jend, 0, sizeof(int)*(NumberOfGrids+1));
  jstart[0] = 0;
  jend[0] = 0;
  int ParticleIterations = 0;
  // Copy shared particles to grids, if any
  if (NumberOfReceives > 0) {
	for (j = 0; j < NumberOfGrids && jend[j] < NumberOfReceives; j++) {
		while (SharedList[jend[j]].grid <= j) {
			jend[j]++;
			if (jend[j] == NumberOfReceives) break;
      	}
      	ParticleIterations++;
	  	jstart[j+1] = jend[j];
	}
  }

  kstart = new int[NumberOfGrids+1];
  kend = new int[NumberOfGrids+1];
  memset(kstart, 0, sizeof(int)*(NumberOfGrids+1));
  memset(kend, 0, sizeof(int)*(NumberOfGrids+1));
  kstart[0] = 0;
  kend[0] = 0;
  int StarParticleIterations = 0;
  // Copy shared particles to grids, if any
  if (StarNumberOfReceives > 0) {
    for (k = 0; k < NumberOfGrids && kend[k] < StarNumberOfReceives; k++) {
      while (StarSharedList[kend[k]].grid <= k) {
		kend[k]++;
		if (kend[k] == StarNumberOfReceives) break;
      }
      StarParticleIterations++;
	  kstart[k+1] = kend[k];
    }
 }

  ParticleCounter1 = 0;
  ParticleCounter2 = 0;
#pragma omp parallel
{
#pragma omp for reduction(+:NumberToMove[:NumberOfProcessors])
    for (j = 0; j < ParticleIterations; j++) {

      
/*      printf("--> j =%d, ParticleIterations= %d, jstart =%d, jend =%d, NumberOfGrids =%d, " 
             "NumberToMove[] =%d/%d, NumberOfReceives =%d\n", 
	     j, ParticleIterations, jstart[j], jend[j], NumberOfGrids, 
	     NumberToMove[0], NumberToMove[1], NumberOfReceives); */
      

      GridPointers[j]->TransferSubgridParticles
	(GridPointers, NumberOfGrids, NumberToMove, ParticleCounter1, jstart[j], jend[j], 
	 SharedList, KeepLocal, ParticlesAreLocal, COPY_IN, TRUE);
      
      //jstart = jend;
    } // ENDFOR grids
  
  /*******************************************************************/
  /******************** Copy stars back to grids. ********************/
  /*******************************************************************/

  // Copy shared stars to grids, if any
#pragma omp for reduction(+:StarsToMove[:NumberOfProcessors])
    for (k = 0; k < StarParticleIterations; k++) {

/*      printf("--> k =%d, StarParticleIterations= %d, kstart =%d, kend =%d, NumberOfGrids =%d, " 
             "NumberToMove[] =%d/%d, NumberOfReceives =%d\n", 
	     k, StarParticleIterations, kstart[k], kend[k], NumberOfGrids, 
	     StarsToMove[0], StarsToMove[1], StarNumberOfReceives);  */
      
      
      GridPointers[k]->TransferSubgridStars
	(GridPointers, NumberOfGrids, StarsToMove, kstart[k], kend[k], 
	 StarSharedList, KeepLocal, ParticlesAreLocal, COPY_IN, TRUE);
      
      //jstart = jend;
    } // ENDFOR grids

} // end omp parallel

  /*******************************************************************/
  /************** Copy Active Particles back to grids. ***************/
  /*******************************************************************/

  jstart = 0;
  jend = 0;

  // Copy shared stars to grids, if any

  if (APNumberOfReceives > 0)
    for (j = 0; j < NumberOfGrids && jend < APNumberOfReceives; j++) {
      while (APSharedList[jend]->ReturnGridID() <= j) {
        jend++;
        if (jend == APNumberOfReceives) break;
      }
      
      GridPointers[j]->TransferSubgridActiveParticles
                    (GridPointers, NumberOfGrids, APNumberToMove, jstart, jend,
                     APSharedList, KeepLocal, ParticlesAreLocal, COPY_IN, TRUE);
      
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
  delete [] APNumberToMove;
  return SUCCESS;
}

