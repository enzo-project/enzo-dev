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
#include <time.h>

#include "EnzoTiming.h" 
#include "ErrorExceptions.h"
#include "performance.h"
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
 
/* function prototypes */
 
void AddLevel(LevelHierarchyEntry *LevelArray[], HierarchyEntry *Grid,
	      int level);
int FindSubgrids(HierarchyEntry *Grid, int level, int &TotalFlaggedCells,
		 int &FlaggedGrids);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int ReportMemoryUsage(char *header = NULL);
int DepositParticleMassFlaggingField(LevelHierarchyEntry* LevelArray[],
				     int level, bool AllLocal);
int CommunicationShareGrids(HierarchyEntry *GridHierarchyPointer[], int grids,
			    int ShareParticles = TRUE); 
int CommunicationLoadBalanceGrids(HierarchyEntry *GridHierarchyPointer[],
				  int NumberOfGrids, int MoveParticles = TRUE);
int LoadBalanceHilbertCurve(HierarchyEntry *GridHierarchyPointer[],
			    int NumberOfGrids, int MoveParticles = TRUE);
int CommunicationTransferSubgridParticles(LevelHierarchyEntry *LevelArray[],
					  TopGridData *MetaData, int level);
int DetermineSubgridSizeExtrema(long_int NumberOfCells, int level, int MaximumStaticSubgridLevel);
int CommunicationTransferParticles(grid *GridPointer[], int NumberOfGrids,
				   int TopGridDims[]);
int CommunicationTransferStars(grid *GridPointer[], int NumberOfGrids,
			       int TopGridDims[]);
int CommunicationCollectParticles(LevelHierarchyEntry *LevelArray[], int level,
				  bool ParticlesAreLocal,
				  bool SyncNumberOfParticles, 
				  bool MoveStars, int CollectMode);
int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],
				       int NumberOfGrids);
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
int CopyZonesFromOldGrids(LevelHierarchyEntry *OldGrids, 
			  TopGridData *MetaData,
			  ChainingMeshStructure ChainingMesh);
#ifdef TRANSFER
int SetSubgridMarker(TopGridData &MetaData, 
		     LevelHierarchyEntry *LevelArray[], int level,
		     int UpdateReplicatedGridsOnly);
#endif
double ReturnWallTime(void);

void fpcol(Eflt64 *x, int n, int m, FILE *log_fptr);
bool _first = true;
static double RHperf[16];

#define NO_RH_PERF


/* RebuildHierarchy function */
 
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level)
{

  if (LevelSubCycleCount[level] % RebuildHierarchyCycleSkip[level]) {
    return SUCCESS;
  }

  if (ConductionDynamicRebuildHierarchy) {
    if (TimeSinceRebuildHierarchy[level] < dtRebuildHierarchy[level]) {
      return SUCCESS;
    }
    else {
      for (int i = level;i <= MaximumRefinementLevel;i++) {
        dtRebuildHierarchy[i] = -1.0;
        TimeSinceRebuildHierarchy[i] = 0.0;
      }
    }
  }

  double tt0, tt1, tt2, tt3;
 
  /* declarations */

  int dbx = 0;
 
  LCAPERF_START("RebuildHierarchy");
  TIMER_START("RebuildHierarchy");

  if (debug) printf("RebuildHierarchy: level = %"ISYM"\n", level);
  ReportMemoryUsage("Rebuild pos 1");
 
  bool ParticlesAreLocal, SyncNumberOfParticles = true;
  bool MoveStars = true;
  int i, j, k, grids, grids2, subgrids, MoveParticles, ncells;
  int TotalFlaggedCells, FlaggedGrids;
  FLOAT ZeroVector[MAX_DIMENSION];
  LevelHierarchyEntry *Temp;
  HierarchyEntry *GridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
 
  for (i = 0; i < MAX_DIMENSION; i++)
    ZeroVector[i] = 0;

  if (_first) {
    for (i = 0; i < 16; i++)
      RHperf[i] = 0;
    _first = false;
  }
 
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


  /* Because we're storing particles in "empty" grids that are local
     to the subgrid, keep track of the number of particles stored
     locally.  Zero out NumberOfParticles on other processors, then
     collect the subgrid particles. */

  for (i = level; i < MAX_DEPTH_OF_HIERARCHY; i++)
    for (Temp = LevelArray[i]; Temp; Temp = Temp->NextGridThisLevel)
      if (MyProcessorNumber != Temp->GridData->ReturnProcessorNumber()) {
	Temp->GridData->SetNumberOfParticles(0);
	Temp->GridData->SetNumberOfStars(0);
      }

  /* The dynamic grids should be distributed enough to store the
     particles on each grid, so we'll collect the particles at the
     finest static subgrid level, which we find here. */

  int MaximumStaticSubgridLevel = -1;
  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    MaximumStaticSubgridLevel = max(MaximumStaticSubgridLevel,
				    StaticRefineRegionLevel[i]);

  /* Calculate number of cells on each level */

  long_int NumberOfCells[MAX_DEPTH_OF_HIERARCHY];
  if (SubgridSizeAutoAdjust == TRUE) {
    for (i = level; i < MAX_DEPTH_OF_HIERARCHY; i++) {
      NumberOfCells[i] = 0;
      for (Temp = LevelArray[i]; Temp; Temp = Temp->NextGridThisLevel)
	if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber())
	  NumberOfCells[i] += Temp->GridData->GetActiveSize();
    }
    CommunicationAllSumValues(NumberOfCells, MAX_DEPTH_OF_HIERARCHY);
  }

  tt0 = ReturnWallTime();
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

	GridParent[j]->GridData->MoveAllStars(grids2, ContigiousGridList, 
					      MetaData->TopGridDims[0]);
	GridParent[j]->GridData->MoveAllParticles(grids2, ContigiousGridList);

#ifdef TRANSFER   
	/* Rescue all PhotonPackages before the subgrids are deleted. */
	GridParent[j]->GridData->MoveAllPhotonPackages(grids2, ContigiousGridList);
#endif // TRANSFER
	
      } // end: if grid pointer valid
 
  } // end: loop over levels
  tt1 = ReturnWallTime();
  RHperf[0] += tt1-tt0;


  /* If the initial level is finer than the finest level with static
     subgrids, we must collect all of the particles on the grids' host
     processor before rebuilding.  Before MoveAllParticles did
     this. */

  tt0 = ReturnWallTime();
  if (level > MaximumStaticSubgridLevel) {
    ParticlesAreLocal = false;
    SyncNumberOfParticles = false;
    CommunicationCollectParticles(LevelArray, level, ParticlesAreLocal, 
				  SyncNumberOfParticles, MoveStars,
				  SIBLINGS_ONLY);
    ParticlesAreLocal = true;
    SyncNumberOfParticles = true;
  }
  tt1 = ReturnWallTime();
  RHperf[2] += tt1-tt0;
 
  /* --------------------------------------------------------------------- */
  /* if this is level 0 then transfer particles between grids. */

  tt0 = ReturnWallTime();
  if (dbx) fprintf(stderr, "Rebuild pos 2\n");
  ReportMemoryUsage("Rebuild pos 2");
  if (level == 0) {
    grids = 0;
    Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->DebugCheck("Before TransferParticles");
      GridPointer[grids++] = Temp->GridData;
      Temp = Temp->NextGridThisLevel;
    }

    CommunicationTransferParticles(GridPointer, grids, MetaData->TopGridDims);
    CommunicationTransferStars(GridPointer, grids, MetaData->TopGridDims);

    /* We need to collect particles again */

    if (level > MaximumStaticSubgridLevel) {
      ParticlesAreLocal = false;
      SyncNumberOfParticles = true;
      CommunicationCollectParticles(LevelArray, level, ParticlesAreLocal, 
				    SyncNumberOfParticles, MoveStars,
				    SIBLINGS_ONLY);
      ParticlesAreLocal = true;
      SyncNumberOfParticles = true;
    }


  } // ENDIF level 0
  tt1 = ReturnWallTime();
  RHperf[1] += tt1-tt0;


  /* --------------------------------------------------------------------- */
  /* Transfer particle between grids on this level to make sure that
     each grid contains all of the particles that it should (in case
     some particles have moved outside of this grid's boundaries).
     (This does the same as the previous routine but is not optimized
     for level 0 and does not make sure that all particles have been
     transfered).  This must be done after CommunicationCollectParticles.
  */

  if (MoveParticlesBetweenSiblings && 
      level > max(MaximumStaticSubgridLevel,0))
    CommunicationTransferSubgridParticles(LevelArray, MetaData, level);



  /* --------------------------------------------------------------------- */
  /* For dynamic hierarchies, rebuild the grid structure. */
 
  if (dbx) fprintf(stderr, "Rebuild pos 3\n");
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


      /* Determine the subgrid minimum and maximum sizes, if
         requested.

         If we are initializing and on our first trip through
         the hierachy, use the number of cells on the parent
         level to estimate the grid efficiency parameters
      */

      if (NumberOfCells[i+1] == 0)
        ncells = NumberOfCells[i];
      else
        ncells = NumberOfCells[i+1];

      DetermineSubgridSizeExtrema(ncells, i+1, MaximumStaticSubgridLevel+1);

      DetermineSubgridSizeExtrema(NumberOfCells[i+1], i+1, 
				  MaximumStaticSubgridLevel+1);

      /* 3a) Generate an array of grids on this level. */
 
//??      HierarchyEntry *GridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
 
      grids = 0;
      Temp = LevelArray[i];
      while (Temp != NULL) {
	GridHierarchyPointer[grids++] = Temp->GridHierarchyEntry;
	Temp                          = Temp->NextGridThisLevel;
      }

      /* 3b.1) Loop over grids, creating the particle mass flagging
	 field by considering particles on all processors.  They
	 aren't on the local processor anymore to distribute memory
	 usage when running large nested grid runs.  We want to do
	 this in a separate loop from FindSubgrids to improve parallel
	 efficiency. */

      ParticlesAreLocal = (i > MaximumStaticSubgridLevel);
      MoveParticles = (ParticlesAreLocal) ? TRUE : FALSE;

      tt0 = ReturnWallTime();
      DepositParticleMassFlaggingField(LevelArray, i, ParticlesAreLocal);
      tt1 = ReturnWallTime();
      RHperf[3] += tt1-tt0;

      /* 3b.2) Loop over grids creating new (but empty!) subgrids
	 (This also properly fills out the GridHierarchy tree). */

      tt0 = ReturnWallTime();
      TotalFlaggedCells = FlaggedGrids = 0;
      for (j = 0; j < grids; j++)
	FindSubgrids(GridHierarchyPointer[j], i, TotalFlaggedCells, FlaggedGrids);
      CommunicationSumValues(&TotalFlaggedCells, 1);
      CommunicationSumValues(&FlaggedGrids, 1);
      if (debug)
	printf("RebuildHierarchy[%"ISYM"]: "
	       "Flagged %"ISYM"/%"ISYM" grids. %"ISYM" flagged cells\n", 
	       i, FlaggedGrids, grids, TotalFlaggedCells);
      tt1 = ReturnWallTime();
      RHperf[4] += tt1-tt0;

      /* Create a temporary array of the new subgrids (which are on this 
	 processor) for the next step. */

      HierarchyEntry *Temp2;
      HierarchyEntry *SubgridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
      for (j = 0; j < MAX_NUMBER_OF_SUBGRIDS; j++){
        SubgridHierarchyPointer[j] = NULL;
      }
      subgrids = 0;
      for (j = 0; j < grids; j++) {
	Temp2 = GridHierarchyPointer[j]->NextGridNextLevel;
	while (Temp2 != NULL) {
	  SubgridHierarchyPointer[subgrids++] = Temp2;
	  Temp2 = Temp2->NextGridThisLevel;
	}
      }

      /* Share the new grids among processors. */

      tt0 = ReturnWallTime();
      CommunicationShareGrids(GridHierarchyPointer, grids, MoveParticles); 
      tt1 = ReturnWallTime();
      RHperf[5] += tt1-tt0;

      /* 3c) Combine the many linked-lists of subgrids into the LevelArray
	 linked list. */
 
      tt0 = ReturnWallTime();
      for (j = 0; j < grids; j++)
	if (GridHierarchyPointer[j]->NextGridNextLevel != NULL)
	  AddLevel(LevelArray, GridHierarchyPointer[j]->NextGridNextLevel,i+1);
      tt1 = ReturnWallTime();
      RHperf[6] += tt1-tt0;

      /* 3g) loop over parent, and copy particles to new grids
	     (all local to this processor) . */

      /* JHW (May 2009) For levels with static subgrids, the particles
	 are still on the same processor as they were before we
	 entered RebuildHierarchy.  We only collect them on the
	 correct processor after everything is rebuilt when we reach
	 the finest level with static subgrids and after load
	 balancing to distribute memory usage.. */

      tt0 = ReturnWallTime();
      CommunicationCollectParticles(LevelArray, i, ParticlesAreLocal,
				    SyncNumberOfParticles, MoveStars,
				    SUBGRIDS_LOCAL);
      tt1 = ReturnWallTime();
      RHperf[7] += tt1-tt0;

      /* 3d) Create an array of the new subgrids. */
 
      subgrids = 0;
      Temp = LevelArray[i+1];
      while (Temp != NULL) {
	SubgridHierarchyPointer[subgrids++] = Temp->GridHierarchyEntry;
	Temp                                = Temp->NextGridThisLevel;
      }
 
      //Old fine grids are necessary during the interpolation for ensuring DivB = 0 with MHDCT
      //Note that this is a loop the size of N_{new sub grids} * N_{old sub grids}.  Fast Sib locator 
      //needs to be employed here.
      if( UseMHDCT ){
        for (j = 0; j < subgrids; j++) {
           if(SubgridHierarchyPointer[j]->GridData->MHD_SendOldFineGrids(
                  TempLevelArray[i+1],SubgridHierarchyPointer[j]->ParentGrid->GridData, MetaData) == FALSE ){
             ENZO_FAIL("Error in SendOldFineGrids");
              }
        }
      }

      /* 3e) For each new subgrid, interpolate from parent and then
	 copy from old subgrids.  For each old subgrid, decrement the
	 Overlap counter, deleting the grid which it reaches zero. */
      
      tt0 = ReturnWallTime();
      for (j = 0; j < subgrids; j++) {
	SubgridHierarchyPointer[j]->ParentGrid->GridData->
	  DebugCheck("Rebuild parent");

        if (RandomForcing) { //AK
          SubgridHierarchyPointer[j]->GridData->AppendForcingToBaryonFields();
          SubgridHierarchyPointer[j]->ParentGrid->GridData->
	    AppendForcingToBaryonFields();
        }

	SubgridHierarchyPointer[j]->GridData->InterpolateFieldValues
	  (SubgridHierarchyPointer[j]->ParentGrid->GridData
		,TempLevelArray[i+1],
	   MetaData);

        if (RandomForcing) { //AK
          SubgridHierarchyPointer[j]->GridData->RemoveForcingFromBaryonFields();
          SubgridHierarchyPointer[j]->ParentGrid->GridData->
	    RemoveForcingFromBaryonFields();
        }

	SubgridHierarchyPointer[j]->GridData->DebugCheck("Rebuild child");
      }
      tt1 = ReturnWallTime();
      RHperf[8] += tt1-tt0;
 
      /* 3f) Loop over the old grids and copy the data into the new grids.
             This is done in two steps in order to speed up the search,
             first we generate a chaining mesh with a linked list that lists
	     all new subgrids by location.  Then, we use that to generate a
	     sibling list for each of the old grids, and copy the data. */

      tt0 = ReturnWallTime();
      if (dbx) fprintf(stderr, "RH: Initialize FSL \n");
      ChainingMeshStructure ChainingMesh;
      FastSiblingLocatorInitialize(&ChainingMesh, MetaData->TopGridRank,
				   MetaData->TopGridDims);
      tt1 = ReturnWallTime();
      RHperf[9] += tt1-tt0;
 
      /*  Add all the new subgrids to the chaining mesh. */

      if (dbx) fprintf(stderr, "RH: FSL AddGrid entry \n");

      tt0 = ReturnWallTime();
      for (j = 0; j < subgrids; j++)
	SubgridHierarchyPointer[j]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);
      tt1 = ReturnWallTime();
      RHperf[10] += tt1-tt0;

      if (dbx) fprintf(stderr, "RH: FSL AddGrid exit \n");

      /* Copy data from old to new grids */
 
      tt0 = ReturnWallTime();
      CopyZonesFromOldGrids(TempLevelArray[i+1], MetaData, ChainingMesh);
      tt1 = ReturnWallTime();
      RHperf[12] += tt1-tt0;
 
      /* Clean up chaining mesh. */
 
      FastSiblingLocatorFinalize(&ChainingMesh);

      tt0 = ReturnWallTime();
      /* Redistribute grids over processors to Load balance. */
      switch( LoadBalancing ){
      case 1:
      case 2:
      case 3:
	if (i >= LoadBalancingMinLevel && i <= LoadBalancingMaxLevel)
	  CommunicationLoadBalanceGrids(SubgridHierarchyPointer, subgrids, 
					MoveParticles);
	break;
      case 4:
	if (i >= LoadBalancingMinLevel && i <= LoadBalancingMaxLevel)
	  LoadBalanceHilbertCurve(SubgridHierarchyPointer, subgrids, 
				  MoveParticles);
	break;
      default:
	break;
      }
      tt1 = ReturnWallTime();
      RHperf[13] += tt1-tt0;

      /* If this is the finest level with static subgrids, the grids
	 should be distributed enough to collect the particles on each
	 host processor. */

      tt0 = ReturnWallTime();
      if (i == MaximumStaticSubgridLevel)
	for (j = level; j <= MaximumStaticSubgridLevel+1; j++)
	  if (LevelArray[j] != NULL)
	    CommunicationCollectParticles(LevelArray, j, ParticlesAreLocal,
					  SyncNumberOfParticles, MoveStars,
					  SIBLINGS_ONLY);
      tt1 = ReturnWallTime();
      RHperf[14] += tt1-tt0;

      /* 3h) Clean up the LevelHierarchy entries for the old subgrids.
	     Also, we can check to see if any old subgrids were missed. */
 
      while (TempLevelArray[i+1] != NULL) {
	Temp = TempLevelArray[i+1]->NextGridThisLevel;
 
	if (TempLevelArray[i+1]->GridData != NULL)
	  ENZO_FAIL("An old subgrid was not deleted.  Why?");
 
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
		              ZERO_UNDER_SUBGRID_FIELD, float(k+1)) == FAIL)
	      ENZO_FAIL("Error in grid->ZeroSolutionUnderSubgrid.");

	    ToGrids[k] = SubgridHierarchyPointer[k]->GridData;
	  }
 
	  if (GridHierarchyPointer[j]->GridData->MoveSubgridStars(
				 subgrids, ToGrids, FALSE) == FAIL)
	    ENZO_FAIL("Error in grid->MoveSubgridStars.");

	  if (GridHierarchyPointer[j]->GridData->MoveSubgridParticlesFast(
				 subgrids, ToGrids, FALSE) == FAIL)
	    ENZO_FAIL("Error in grid->MoveSubgridParticlesFast.");
 
	}
 
      /* Set boundary conditions. */
 
      LevelHierarchyEntry *Temp = LevelArray[i+1];
      while (Temp != NULL) {
 
	if (Temp->GridData->InterpolateBoundaryFromParent
	    (Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL)
	  ENZO_FAIL("Error in grid->InterpolateBoundaryFromParent.");
 
	Temp = Temp->NextGridThisLevel;
      }
 
    } // end: loop over levels
 
  } // end: if (StaticHierarchy == TRUE)

  /* set grid IDs */

  for (i = level; i < MAX_DEPTH_OF_HIERARCHY-1; i++)
    for (Temp = LevelArray[i], j = 0; Temp; Temp = Temp->NextGridThisLevel, j++)
      Temp->GridData->SetGridID(j);

  /* update all SubgridMarkers */

#ifdef TRANSFER
  SetSubgridMarker(*MetaData, LevelArray, level, FALSE);
#endif /* TRANSFER  */
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[1] += endtime - starttime;
  counter[1] ++;
#endif /* MPI_INSTRUMENTATION */
 
  /* Done for this level. */

#ifdef RH_PERF
#ifdef USE_MPI
  CommunicationReduceValues(RHperf, 16, MPI_MAX);
#endif
  //CommunicationSumValues(RHperf, 16);
  if (debug) fpcol(RHperf, 16, 16, stdout);
#endif /* RH_PERF */
  ReportMemoryUsage("Rebuild pos 4");
  TIMER_STOP("RebuildHierarchy");
  LCAPERF_STOP("RebuildHierarchy");
  return SUCCESS;
 
}
