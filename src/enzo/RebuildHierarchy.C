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
int DepositParticleMassFlaggingField(LevelHierarchyEntry* LevelArray[],
				     int level, bool AllLocal);
int CommunicationShareGrids(HierarchyEntry *GridHierarchyPointer[], int grids,
			    int ShareParticles = TRUE);
int CommunicationLoadBalanceGrids(HierarchyEntry *GridHierarchyPointer[],
				  int NumberOfGrids, int MoveParticles = TRUE);
int CommunicationTransferParticles(grid *GridPointer[], int NumberOfGrids);
int CommunicationCollectParticles(LevelHierarchyEntry *LevelArray[], int level,
				  bool ParticlesAreLocal, int CollectMode);
int CommunicationSyncNumberOfParticles(HierarchyEntry *GridHierarchyPointer[],
				       int NumberOfGrids);
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
#ifdef TRANSFER
int SetSubgridMarker(TopGridData &MetaData, 
		     LevelHierarchyEntry *LevelArray[], int level);
#endif
double ReturnWallTime(void);

#define NO_TIME_MESSAGING
 
/* RebuildHierarchy function */
 
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level)
{
 
  /* declarations */

  time_t rawtime;
  struct tm* timeinfo;
  int garbage;

  int dbx = 0;
 
  if (debug) printf("RebuildHierarchy: level = %"ISYM"\n", level);
  ReportMemoryUsage("Rebuild pos 1");
 
  bool ParticlesAreLocal;
  int i, j, k, grids, grids2, subgrids, MoveParticles;
  FLOAT ZeroVector[MAX_DIMENSION];
  LevelHierarchyEntry *Temp;
  HierarchyEntry *GridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
 
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

  /* Because we're storing particles in "empty" grids that are local
     to the subgrid, keep track of the number of particles stored
     locally.  Zero out NumberOfParticles on other processors, then
     collect the subgrid particles. */

  for (i = level; i < MAX_DEPTH_OF_HIERARCHY; i++)
    for (Temp = LevelArray[i]; Temp; Temp = Temp->NextGridThisLevel)
      if (MyProcessorNumber != Temp->GridData->ReturnProcessorNumber())
	Temp->GridData->SetNumberOfParticles(0);

  /* The dynamic grids should be distributed enough to store the
     particles on each grid, so we'll collect the particles at the
     finest static subgrid level, which we find here. */

  int MaximumStaticSubgridLevel = -1;
  for (i = 0; i < MAX_STATIC_REGIONS; i++)
    MaximumStaticSubgridLevel = max(MaximumStaticSubgridLevel,
				    StaticRefineRegionLevel[i]);

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

	if (GridParent[j]->GridData->
	    MoveAllStars(grids2, ContigiousGridList, 
			 MetaData->TopGridDims[0]) == FAIL) {
	  fprintf(stderr, "Error in grid->MoveAllStars.\n");
	  return FAIL;
	}

	if (GridParent[j]->GridData->MoveAllParticles(grids2,
					   ContigiousGridList) == FAIL) {
	  fprintf(stderr, "Error in grid->MoveAllParticles.\n");
	  return FAIL;
	}

#ifdef TRANSFER   
	/* Rescue all PhotonPackages before the subgrids are deleted. */
	if (GridParent[j]->GridData->
	    MoveAllPhotonPackages(grids2, ContigiousGridList) == FAIL) {
	  fprintf(stderr, "Error in grid->MoveAllPhotonPackages(%"ISYM").\n", level);
	  return FAIL;
	}
#endif // TRANSFER
	
      } // end: if grid pointer valid
 
  } // end: loop over levels
 
  /* --------------------------------------------------------------------- */
  /* if this is level 0 then transfer particles between grids. */

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

    TIME_MSG("Entering CommunicationTransferParticles.");

    if (CommunicationTransferParticles(GridPointer, grids) == FAIL) {
      fprintf(stderr, "Error in CommunicationTransferParticles.\n");
      return FAIL;
    }

  } // ENDIF level 0

  /* If the initial level is finer than the finest level with static
     subgrids, we must collect all of the particles on the grids' host
     processor before rebuilding.  Before MoveAllParticles did
     this. */

  if (level > MaximumStaticSubgridLevel) {
    TIME_MSG("Entering CommunicationCollectParticles.");
    ParticlesAreLocal = false;
    if (CommunicationCollectParticles(LevelArray, level, ParticlesAreLocal, 
				      SIBLINGS_ONLY) == FAIL) {
      fprintf(stderr, "Error in CommunicationCollectParticles(root).\n");
      return FAIL;
    }
    ParticlesAreLocal = true;
    TIME_MSG("After CommunicationCollectParticles.");
  }

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

    double tt0, tt1, tt2;
 
    for (i = level; i < MAX_DEPTH_OF_HIERARCHY-1; i++) {
 
      /* If there are no grids on this level, exit. */
 
      if (LevelArray[i] == NULL)
	break;

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

      TIME_MSG("Depositing particle mass flagging field");

      if (DepositParticleMassFlaggingField(LevelArray, i, 
					   ParticlesAreLocal) == FAIL) {
	fprintf(stderr, "Error in DepositParticleMassFlaggingField.\n");
	return FAIL;
      }

      /* 3b.2) Loop over grids creating new (but empty!) subgrids
	 (This also properly fills out the GridHierarchy tree). */

      TIME_MSG("Finding new subgrids");

      for (j = 0; j < grids; j++)
	if (FindSubgrids(GridHierarchyPointer[j], i) == FAIL) {
	  fprintf(stderr, "Error in FindSubgrids.\n");
	  return FAIL;
	}

      TIME_MSG("Found new subgrids");

      /* Create a temporary array of the new subgrids (which are on this
	 processor) for the next step. */
 
      HierarchyEntry *Temp2;
      HierarchyEntry *SubgridHierarchyPointer[MAX_NUMBER_OF_SUBGRIDS];
      subgrids = 0;
      for (j = 0; j < grids; j++) {
	Temp2 = GridHierarchyPointer[j]->NextGridNextLevel;
	while (Temp2 != NULL) {
	  SubgridHierarchyPointer[subgrids++] = Temp2;
	  Temp2 = Temp2->NextGridThisLevel;
	}
      }
 
      /* Share the new grids amoung processors. */

      TIME_MSG("Sharing grids");
      CommunicationShareGrids(GridHierarchyPointer, grids, MoveParticles);
      TIME_MSG("Finished sharing grids");

      /* 3c) Combine the many linked-lists of subgrids into the LevelArray
	 linked list. */
 
      for (j = 0; j < grids; j++)
	if (GridHierarchyPointer[j]->NextGridNextLevel != NULL)
	  AddLevel(LevelArray, GridHierarchyPointer[j]->NextGridNextLevel,i+1);
 
 
      /* 3g) loop over parent, and copy particles to new grids
	     (all local to this processor) . */

      /* JHW (May 2009) For levels with static subgrids, the particles
	 are still on the same processor as they were before we
	 entered RebuildHierarchy.  We only collect them on the
	 correct processor after everything is rebuilt when we reach
	 the finest level with static subgrids and after load
	 balancing to distribute memory usage.. */

      tt0 = ReturnWallTime();
      if (CommunicationCollectParticles(LevelArray, i, ParticlesAreLocal,
					SUBGRIDS_LOCAL) == FAIL) {
	fprintf(stderr, "Error in CommunicationCollectParticles(subgrids).\n");
	return FAIL;
      }
      TIME_MSG("Moved subgrid particles");
      if (MyProcessorNumber == ROOT_PROCESSOR) {
	tt1 = ReturnWallTime();
	printf("RebuildHierarchy[AA]: Took %lg seconds to move particles to"
	       " subgrids.\n", tt1-tt0);
      }

      /* 3d) Create an array of the new subgrids. */
 
      subgrids = 0;
      Temp = LevelArray[i+1];
      while (Temp != NULL) {
	SubgridHierarchyPointer[subgrids++] = Temp->GridHierarchyEntry;
	Temp                                = Temp->NextGridThisLevel;
      }
 
      /* 3e) For each new subgrid, interpolate from parent and then
	 copy from old subgrids.  For each old subgrid, decrement the
	 Overlap counter, deleting the grid which it reaches zero. */

      TIME_MSG("Copying zones");
 
      for (j = 0; j < subgrids; j++) {
	SubgridHierarchyPointer[j]->ParentGrid->GridData->
	  DebugCheck("Rebuild parent");

        if (RandomForcing) { //AK
          SubgridHierarchyPointer[j]->GridData->AppendForcingToBaryonFields();
          SubgridHierarchyPointer[j]->ParentGrid->GridData->
	    AppendForcingToBaryonFields();
        }

	SubgridHierarchyPointer[j]->GridData->InterpolateFieldValues
	  (SubgridHierarchyPointer[j]->ParentGrid->GridData);

        if (RandomForcing) { //AK
          SubgridHierarchyPointer[j]->GridData->RemoveForcingFromBaryonFields();
          SubgridHierarchyPointer[j]->ParentGrid->GridData->
	    RemoveForcingFromBaryonFields();
        }

	SubgridHierarchyPointer[j]->GridData->DebugCheck("Rebuild child");
      }
 
      /* 3f) Loop over the old grids and copy the data into the new grids.
             This is done in two steps in order to speed up the search,
             first we generate a chaining mesh with a linked list that lists
	     all new subgrids by location.  Then, we use that to generate a
	     sibling list for each of the old grids, and copy the data. */

      if (dbx) fprintf(stderr, "RH: Initialize FSL \n");
      ChainingMeshStructure ChainingMesh;
      FastSiblingLocatorInitialize(&ChainingMesh, MetaData->TopGridRank,
				   MetaData->TopGridDims);
      SiblingGridList SiblingList;
 
      /*  Add all the new subgrids to the chaining mesh. */

      if (dbx) fprintf(stderr, "RH: FSL AddGrid entry \n");

      for (j = 0; j < subgrids; j++)
	SubgridHierarchyPointer[j]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);

      if (dbx) fprintf(stderr, "RH: FSL AddGrid exit \n");
 
      /* Loop over the old grids. */
 
      Temp = TempLevelArray[i+1];
      while (Temp != NULL) {
 
	/* Find sibling grids. */
 
	if (Temp->GridData->FastSiblingLocatorFindSiblings(
                              &ChainingMesh, &SiblingList,
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition) == FAIL) {
	  fprintf(stderr, "Error in grid->FastSiblingLocatorFindSiblings.\n");
	  return FAIL;
	}
 
	/* For each of the sibling grids, copy data. */
 
	for (j = 0; j < SiblingList.NumberOfSiblings; j++) {
	  if (SiblingList.GridList[j]->CopyZonesFromGrid(
		                        Temp->GridData, ZeroVector) == FAIL) {
	    fprintf(stderr, "Error in grid->CopyZonesFromGridCountOnly.\n");
	    return FAIL;
	  }
	}
 
	/* delete old grid and sibling data. */
 
	delete Temp->GridData;
	Temp->GridData = NULL;
	delete [] SiblingList.GridList;
 
	/* Next old grid. */
 
	Temp = Temp->NextGridThisLevel;
      }
 
      /* Clean up chaining mesh. */
 
      FastSiblingLocatorFinalize(&ChainingMesh);
 
      /* Redistribute grids over processors to Load balance. */
      switch( LoadBalancing ){
      case 1:
	TIME_MSG("Load balancing");
	CommunicationLoadBalanceGrids(SubgridHierarchyPointer, subgrids, 
				      MoveParticles);
	TIME_MSG("Finished load balancing");
	break;
      default:
	
	break;
      }


      /* If this is the finest level with static subgrids, the grids
	 should be distributed enough to collect the particles on each
	 host processor. */

      if (i == MaximumStaticSubgridLevel) {
	tt0 = ReturnWallTime();
	TIME_MSG("Collecting particles");
	for (j = level; j <= MaximumStaticSubgridLevel+1; j++)
	  if (LevelArray[j] != NULL)
	    if (CommunicationCollectParticles(LevelArray, j, ParticlesAreLocal,
					      SIBLINGS_ONLY) == FAIL) {
	      fprintf(stderr, "Error in CommunicationCollectParticles(siblings).\n");
	      return FAIL;
	    }
	//CommunicationSyncNumberOfParticles(SubgridHierarchyPointer, subgrids);
	TIME_MSG("Finished collecting particles");
	if (MyProcessorNumber == ROOT_PROCESSOR) {
	  tt1 = ReturnWallTime();
	  printf("RebuildHierarchy[AA]: Took %lg seconds to move particles "
		 "to correct processor.\n", tt1-tt0);
	}
      }

      /* 3h) Clean up the LevelHierarchy entries for the old subgrids.
	     Also, we can check to see if any old subgrids were missed. */
 
      while (TempLevelArray[i+1] != NULL) {
	Temp = TempLevelArray[i+1]->NextGridThisLevel;
 
	if (TempLevelArray[i+1]->GridData != NULL) {
	  fprintf(stderr, "An old subgrid was not deleted.  Why?\n");
	  return FAIL;
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
	      fprintf(stderr, "Error in grid->ZeroSolutionUnderSubgrid.\n");
	      return FAIL;
	    }
	    ToGrids[k] = SubgridHierarchyPointer[k]->GridData;
	  }
 
	  if (GridHierarchyPointer[j]->GridData->MoveSubgridStars(
				 subgrids, ToGrids, FALSE) == FAIL) {
	    fprintf(stderr, "Error in grid->MoveSubgridStars.\n");
	    return FAIL;
	  }

	  if (GridHierarchyPointer[j]->GridData->MoveSubgridParticlesFast(
				 subgrids, ToGrids, FALSE) == FAIL) {
	    fprintf(stderr, "Error in grid->MoveSubgridParticlesFast.\n");
	    return FAIL;
	  }
 
	}
 
      /* Set boundary conditions. */
 
      LevelHierarchyEntry *Temp = LevelArray[i+1];
      while (Temp != NULL) {
 
	if (Temp->GridData->InterpolateBoundaryFromParent
	    (Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->InterpolateBoundaryFromParent.\n");
	  return FAIL;
	}
 
	Temp = Temp->NextGridThisLevel;
      }
 
    } // end: loop over levels
 
  } // end: if (StaticHierarchy == TRUE)

  /* update all SubgridMarkers */

#ifdef TRANSFER
  if (SetSubgridMarker(*MetaData, LevelArray, 0) == FAIL) {
    fprintf(stderr, "Error in SetSubgridMarker from RebuildHierarchy.\n");
    return FAIL;
  }
#endif /* TRANSFER  */
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[1] += endtime - starttime;
  counter[1] ++;
#endif /* MPI_INSTRUMENTATION */
 
  /* Done for this level. */
 
  ReportMemoryUsage("Rebuild pos 4");
  return SUCCESS;
 
}
