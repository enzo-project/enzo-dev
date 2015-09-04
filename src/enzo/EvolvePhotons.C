/***********************************************************************
/
/  EVOLVE PHOTONS FUNCTION
/
/  written by: Tom Abel
/  date:       May 2004
/  modified1:  November 2005 by John Wise (parallelized it)
/                
/
/  PURPOSE:
/    This routine is the main photon evolution function. 
/    From here we call first the routines that control the emission:
/    grid::Shine
/  while (stil_photons_to_update)
/    transport all photon packages local to a grid
/    communicate all surviving photon packages to their new parent grids. 
/  endwhile  
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
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
#include "GroupPhotonList.h"
#include "PhotonCommunication.h"
#include "CommunicationUtilities.h"

/* function prototypes */
void my_exit(int status);
int RadiationFieldCalculateRates(FLOAT Time);
int CommunicationReceiverPhotons(LevelHierarchyEntry *LevelArray[],
				 bool local_transport,
				 int &keep_transporting);
int CommunicationTransferPhotons(LevelHierarchyEntry *LevelArray[], 
				 ListOfPhotonsToMove **AllPhotons, 
				 char *kt_global,
				 int &keep_transporting);
int RadiativeTransferLoadBalanceRevert(HierarchyEntry **Grids[], int *NumberOfGrids);
int CommunicationLoadBalancePhotonGrids(HierarchyEntry **Grids[], int *NumberOfGrids,
					int FirstTimeAfterRestart);
int RadiativeTransferMoveLocalPhotons(ListOfPhotonsToMove **AllPhotons,
				      int &keep_transporting);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int InitializePhotonMessages(void);
int InitializePhotonCommunication(void);
int FinalizePhotonCommunication(void);
int KeepTransportingInitialize(char* &kt_global, bool initial_call);
int KeepTransportingFinalize(char* &kt_global, int keep_transporting);
int KeepTransportingCheck(char* &kt_global, int &keep_transporting);
int KeepTransportingSend(int keep_transporting);
RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS);
PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
int CreateSourceClusteringTree(int nShine, SuperSourceData *SourceList,
			       LevelHierarchyEntry *LevelArray[]);
void PrintSourceClusteringTree(SuperSourceEntry *leaf);
int CommunicationSyncNumberOfPhotons(LevelHierarchyEntry *LevelArray[]);
int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove,
				     int level);
int StarParticleFindAll(LevelHierarchyEntry *LevelArray[], Star *&AllStars);
int SetSubgridMarker(TopGridData &MetaData, 
		     LevelHierarchyEntry *LevelArray[], int level,
		     int UpdateReplicatedGridsOnly);
void PrintMemoryUsage(char *str);
void fpcol(Eflt64 *x, int n, int m, FILE *log_fptr);
double ReturnWallTime();

#ifdef USE_MPI
int InitializePhotonReceive(int max_size, bool local_transport,
			    MPI_Datatype MPI_PhotonType);
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_PhotonList;
#endif

//#define NONBLOCKING_RT_OFF  // moved to a compile-time define
#define REPORT_PERF
#define MAX_ITERATIONS 5

#ifdef REPORT_PERF
#define START_PERF() tt0 = ReturnWallTime();
#else
#define START_PERF() ;
#endif
#ifdef REPORT_PERF
#define END_PERF(A) \
  tt1 = ReturnWallTime(); \
  PerfCounter[A] += tt1-tt0;
#else
#define END_PERF(A) ;
#endif

/* EvolvePhotons function */
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		  Star *&AllStars, FLOAT GridTime, int level, int LoopTime)
{

  struct ThinGridList
  {
    grid *ThisGrid;
    ThinGridList *NextGrid;
  };

  bool FirstTime = true;

#ifdef REPORT_PERF
  double ep0, tt0, tt1, PerfCounter[14];
  ep0 = ReturnWallTime();
  for (int i = 0; i < 14; i++)
    PerfCounter[i] = 0;
#endif

  if (!RadiativeTransfer)
    return SUCCESS;

  /* Only call on the finest level */

  if (LevelArray[level+1] != NULL && LoopTime)
    return SUCCESS;

//  printf("GridTime = %f, PhotonTime = %f, dtPhoton = %g (Loop = %d)\n",
//	 GridTime, PhotonTime, dtPhoton, (GridTime >= PhotonTime));

  if (dtPhoton < 0)
    return SUCCESS;  

  if (GridTime <= PhotonTime)
    return SUCCESS;

  LCAPERF_START("EvolvePhotons");

  /* Declarations */

#ifdef USE_MPI
  if (FirstTimeCalled) {
    MPI_Type_contiguous(sizeof(GroupPhotonList), MPI_BYTE, &MPI_PhotonList);
    MPI_Type_commit(&MPI_PhotonList);
    FirstTimeCalled = FALSE;
  }
#endif

  /* For early termination with a background, calculate background
     intensities */

  RadiationFieldCalculateRates(PhotonTime+0.5*dtPhoton);

  int i, lvl, GridNum;
  grid *Helper;
  LevelHierarchyEntry *Temp;
  RadiationSourceEntry *RS;

  /* Find the finest grid that hosts each radiation source, and store
     them in a GridList */

  ThinGridList *RS_GridList, *NewNode, *TempGridList, *TailNode, *LastNode, 
    *Destroyer;
  HierarchyEntry **Grids[MAX_DEPTH_OF_HIERARCHY];
  int nGrids[MAX_DEPTH_OF_HIERARCHY];
  for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    if (LevelArray[lvl] != NULL)
      nGrids[lvl] = GenerateGridArray(LevelArray, lvl, &Grids[lvl]);
    else
      nGrids[lvl] = 0;

  RS_GridList = new ThinGridList;
  RS_GridList->NextGrid = NULL;
  TailNode = RS_GridList;
  for (RS = GlobalRadiationSources->NextSource; RS; RS = RS->NextSource) {
    // Search for grid, if not defined. Skip sources that aren't born
    // for PhotonTest
    if (RS->GridID == INT_UNDEFINED && 
	!(RS->CreationTime > PhotonTime && ProblemType == 50)) {
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0; lvl--) {
	for (GridNum = 0; GridNum < nGrids[lvl]; GridNum++) {
	  if (MyProcessorNumber == Grids[lvl][GridNum]->GridData->
	      ReturnProcessorNumber())
	    if (Grids[lvl][GridNum]->GridData->PointInGrid(RS->Position)) {
	      RS->GridID = GridNum;
	      RS->GridLevel = lvl;
	    } // ENDIF PointInGrid
	} // ENDFOR grids
      } // ENDFOR level
    } // ENDIF undefined grid

    /* Now we know which grid, we can add to the GridList */

    NewNode = new ThinGridList;
    if (RS->GridLevel != INT_UNDEFINED)
      NewNode->ThisGrid = Grids[RS->GridLevel][RS->GridID]->GridData;
    else
      NewNode->ThisGrid = NULL;
    NewNode->NextGrid = NULL;
    TailNode->NextGrid = NewNode;
    TailNode = NewNode;

  } // ENDFOR RS

  /**********************************************************************
                       MAIN RADIATION TRANSPORT LOOP
   **********************************************************************/
    
  while (GridTime > PhotonTime) {

  /* Temporarily load balance grids according to the number of ray
     segments.  We'll move the grids back at the end of this
     routine */

    if (RadiativeTransferLoadBalance) {
      CommunicationLoadBalancePhotonGrids(Grids, nGrids, 
					  MetaData->FirstTimestepAfterRestart);
    }

    /* Recalculate timestep if this isn't the first loop.  We already
       did this in RadiativeTransferPrepare */

    FLOAT dtLevelAbove;
    if (!FirstTime) {
      dtLevelAbove = LevelArray[level]->GridData->ReturnTimeStep();
      RadiativeTransferComputeTimestep(LevelArray, MetaData, dtLevelAbove, level);
    }

    if (debug && LoopTime == TRUE)
      printf("EvolvePhotons[%"ISYM"]: dt = %"GSYM", Time = %"FSYM", ", 
	     level, dtPhoton, PhotonTime);
      
    /* delete source if we are passed (or before) their lifetime (only
       if not restarting).  We must remove the host grid from the grid
       list as well. */

    RS = GlobalRadiationSources->NextSource;
    LastNode = RS_GridList;
    TempGridList = RS_GridList->NextGrid;
    int NumberOfSources = 0;
    while (RS != NULL) {
      if ( ((RS->CreationTime + RS->LifeTime) < PhotonTime) && LoopTime == TRUE) {  
	if (debug) {
	  fprintf(stdout, "\nEvolvePhotons: Deleted Source on lifetime limit \n");
	  fprintf(stdout, "EvolvePhotons:  %"GSYM" %"GSYM" %"GSYM" \n",
		  RS->CreationTime, RS->LifeTime, PhotonTime);
	}
	RS = DeleteRadiationSource(RS);
	
	// Remove the grid from the list
	LastNode->NextGrid = TempGridList->NextGrid;
	Destroyer = TempGridList;
	TempGridList = TempGridList->NextGrid;
	delete Destroyer;

      } else {
	if (RS->CreationTime < PhotonTime)
	  NumberOfSources++;                 // count sources
	RS = RS->NextSource;
	LastNode = TempGridList;
	TempGridList = TempGridList->NextGrid;
      }
    }

    if (debug) fprintf(stdout, "%"ISYM" SRC(s)\n", NumberOfSources);

    /* Initialize radiation fields */

    START_PERF();
    for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) 
	if (Temp->GridData->InitializeRadiativeTransferFields() == FAIL) {
	  ENZO_FAIL("Error in InitializeRadiativeTransferFields.\n");
	}
    END_PERF(0);

    /* create temperature fields for Compton heating */  

    if (RadiationXRayComptonHeating)  
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) 
	  if (Temp->GridData->InitializeTemperatureFieldForComptonHeating() == FAIL) {  
	    ENZO_FAIL("Error in InitializeTemperatureFieldForComptonHeating.\n");
	  }	

    for (i = 0; i < 4; i++)
      EscapedPhotonCount[i] = 0.0;

    for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++)
      FieldsToInterpolate[i] = FALSE;

    if (NumberOfSources == 0) {
      PhotonTime += dtPhoton;
      continue;
    }    

    /* Create tree that clusters the sources if requested.  While
       creating tree (type SuperSource), compute position of the super
       source in each leaf. */

    START_PERF();
    if (RadiativeTransferSourceClustering == TRUE) {
      CreateSourceClusteringTree(NULL, NULL, LevelArray);
      //PrintSourceClusteringTree(SourceClusteringTree);
    }
    END_PERF(1);

    // first identify sources and let them radiate 
    RS = GlobalRadiationSources->NextSource;
    TempGridList = RS_GridList->NextGrid;

    START_PERF(); 
    while (RS != NULL) {
      if (TempGridList->ThisGrid != NULL)
	TempGridList->ThisGrid->Shine(RS);
      TempGridList = TempGridList->NextGrid;
      RS = RS->NextSource;
    }    // while still sources 
    END_PERF(2);

#ifdef USE_MPI
    if (RadiativeTransferInterpolateField)
      CommunicationAllReduceValues(FieldsToInterpolate,
				   MAX_NUMBER_OF_BARYON_FIELDS, MPI_MAX);
#endif /* USE_MPI */  

    /* Initialize interpolated radiation fields */  

    //  if (RadiativeTransferInterpolateField)
    //    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    //      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
    //	if (Temp->GridData->AllocateInterpolatedRadiation() == FAIL) {
    //	  ENZO_FAIL("Error in grid->AllocateInterpolatedRadiation.\n");
    //	}


    /* Evolve all photons by fixed timestep. */
  
    ListOfPhotonsToMove *PhotonsToMove = new ListOfPhotonsToMove;
    PhotonsToMove->NextPackageToMove = NULL;

    int keep_transporting = 1;
    int local_keep_transporting = 1, last_keep_transporting;
    int secondary_kt_check = TRUE, iteration = 0;
    bool initial_call = true;
    char *kt_global = NULL;

    HierarchyEntry **Temp0;
    int nGrids0 = GenerateGridArray(LevelArray, 0, &Temp0);
    grid **Grids0 = new grid*[nGrids0];
    for (i = 0; i < nGrids0; i++)
      Grids0[i] = Temp0[i]->GridData;

    /* Initialize nonblocking communication */

    START_PERF();
    InitializePhotonCommunication();
    END_PERF(3);

    /* Transport the rays! */

    PrintMemoryUsage("EvolvePhotons -- before loop");

    while (secondary_kt_check == TRUE && iteration++ < MAX_ITERATIONS) {

#ifdef NONBLOCKING_RT
    KeepTransportingInitialize(kt_global, initial_call);
    initial_call = false;
#endif

    while (keep_transporting != NO_TRANSPORT && 
	   keep_transporting != HALT_TRANSPORT) {
#ifndef NONBLOCKING_RT
      InitializePhotonMessages();
#endif
      last_keep_transporting = local_keep_transporting;
      keep_transporting = 0;
      PhotonsToMove->NextPackageToMove = NULL;
      START_PERF();
#ifndef NONBLOCKING_RT
      keep_transporting = 1;
#endif /* !NONBLOCKING_RT */

      if (local_keep_transporting)
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--) {

	//NumberOfGrids = GenerateGridArray(LevelArray, lvl, &Grids);
	for (Temp = LevelArray[lvl], GridNum = 0;
	     Temp; Temp = Temp->NextGridThisLevel, GridNum++) {
	  //for (GridNum = 0; GridNum < NumberOfGrids; GridNum++) {

	  if (Temp->GridHierarchyEntry->ParentGrid != NULL) 
	    Helper = Temp->GridHierarchyEntry->ParentGrid->GridData;
	  else
	    Helper = NULL;

#ifdef BITWISE_IDENTICALITY
	  Temp->GridData->PhotonSortLinkedLists();
#endif
	  Temp->GridData->TransportPhotonPackages
	    (lvl, &PhotonsToMove, GridNum, Grids0, nGrids0, Helper, 
	     Temp->GridData);

	} // ENDFOR grids

	//delete [] Grids;

      }                          // loop over levels
      END_PERF(4);

      if (PhotonsToMove->NextPackageToMove != NULL)
	keep_transporting = 1;
      else
	keep_transporting = 0;

//    printf("EvoPH[P%"ISYM"]: keep_transporting = %"ISYM", PhotonsToMove = %x\n",
//	     MyProcessorNumber, keep_transporting, PhotonsToMove->NextPackageToMove);

      /* Check if there are any photons leaving this grid.  If so, move them. */
      
      START_PERF();
      CommunicationTransferPhotons(LevelArray, &PhotonsToMove, kt_global,
				   keep_transporting);
      END_PERF(5);

      /* When all photons have been traced, all of the paused (to be
	 merged) photons are in their correct grid, merge them */

      int nmerges = 0;
      if (RadiativeTransferSourceClustering && keep_transporting == 0) {
	for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0; lvl--)
	  for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) {
	    nmerges += Temp->GridData->MergePausedPhotonPackages();
	  } // ENDFOR grids
	if (nmerges > 0) keep_transporting = TRUE;
      }

      /* Receive keep_transporting messages and take the MAX */

      START_PERF();
      local_keep_transporting = keep_transporting;
#ifdef NONBLOCKING_RT
      if (keep_transporting != last_keep_transporting)
	KeepTransportingSend(keep_transporting);
      KeepTransportingCheck(kt_global, keep_transporting);
#else /* NONBLOCKING_RT */
      keep_transporting = CommunicationMaxValue(keep_transporting);
#endif
      END_PERF(6);

    }                           //  end while keep_transporting

#ifdef NONBLOCKING_RT    
    KeepTransportingFinalize(kt_global, keep_transporting);
#ifdef USE_MPI
    InitializePhotonReceive(PHOTON_BUFFER_SIZE, true, MPI_PhotonList);
#endif
    CommunicationReceiverPhotons(LevelArray, false, local_keep_transporting);
    secondary_kt_check = CommunicationMaxValue(local_keep_transporting);
#else /* NONBLOCKING_RT */
    secondary_kt_check = FALSE;
#endif

    }  // ENDWHILE secondary keep_transporting check

    FinalizePhotonCommunication();

    /* Move all finished photon packages back to their original place,
       PhotonPackages.  For the adaptive timestep, we don't carryover
       any photons to the next timestep. */

    START_PERF();
    if (RadiativeTransferAdaptiveTimestep)  
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  Temp->GridData->DeletePhotonPackages();  
    else
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  Temp->GridData->MoveFinishedPhotonsBack();
    END_PERF(7);
    PrintMemoryUsage("EvolvePhotons -- deleted photons");

    /* If we're keeping track of photon escape fractions on multiple
       processors, collect photon counts from all processors */

    FILE *fptr;

    if (RadiativeTransferPhotonEscapeRadius > 0) {
      CommunicationSumValues(EscapedPhotonCount, 4);
      if (MyProcessorNumber == ROOT_PROCESSOR) {

	/* Open f_esc file for writing */

	if (TotalEscapedPhotonCount[0] <= 0) {
	  if ((fptr = fopen(PhotonEscapeFilename, "w")) == NULL) {
	    ENZO_VFAIL("Error opening file %s\n", PhotonEscapeFilename)
	  }
	  fprintf(fptr, 
		  "# Time TotalPhotons fesc(0.5rvir) fesc(rvir) fesc(2rvir)\n");
	} else {
	  if ((fptr = fopen(PhotonEscapeFilename, "a")) == NULL) {
	    ENZO_VFAIL("Error opening file %s\n", PhotonEscapeFilename)
	  }
	}

	for (i = 0; i < 4; i++)
	  TotalEscapedPhotonCount[i] += EscapedPhotonCount[i];

	fprintf(fptr, "%"GOUTSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", 
		PhotonTime, TotalEscapedPhotonCount[0],
		TotalEscapedPhotonCount[1] / TotalEscapedPhotonCount[0], 
		TotalEscapedPhotonCount[2] / TotalEscapedPhotonCount[0], 
		TotalEscapedPhotonCount[3] / TotalEscapedPhotonCount[0]);

	fclose(fptr);

      } // ENDIF ROOT_PROCESSOR

    } // ENDIF RTPhotonEscapeRadius

    PhotonTime += dtPhoton;

    delete PhotonsToMove;
    delete [] Grids0;
    delete [] Temp0;

    /* Sync photon counts after ray tracing */

    START_PERF();
    CommunicationSyncNumberOfPhotons(LevelArray);
    END_PERF(11);

    /* Delete baryon fields on temporary "fake" grid on
       ProcessorNumber and revert it back to OriginalProcessorNumber,
       which was saved in CommunicationLoadBalancePhotonGrids.  Photon
       packages must be moved back, too. */

    if (RadiativeTransferLoadBalance)
      RadiativeTransferLoadBalanceRevert(Grids, nGrids);

    /************************************************************************/
    /********************* Coupled rate & energy solver *********************/
    /************************************************************************/

    int debug_store = debug;
    debug = FALSE;

    // Divide the photo-ionization and photo-heating rates by the
    // number of particles (rho * dx^3)
    START_PERF();
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	if (Temp->GridData->RadiationPresent() == TRUE)
	  Temp->GridData->FinalizeRadiationFields();
    END_PERF(8);

    /* Set the optically-thin H2 dissociation rates */

    START_PERF();
    if (RadiativeTransferOpticallyThinH2)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  Temp->GridData->AddH2Dissociation(AllStars);
    END_PERF(10);

    START_PERF();
    if (RadiativeTransferCoupledRateSolver)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->RadiationPresent() == TRUE) {

	    int RTCoupledSolverIntermediateStep = TRUE;
	    Temp->GridData->SolveRateAndCoolEquations(RTCoupledSolverIntermediateStep);

	  } /* ENDIF radiation */
    END_PERF(9);

    /* Clean up temperature field */

    if (RadiationXRayComptonHeating)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->FinalizeTemperatureFieldForComptonHeating() == FAIL) {  
	    ENZO_FAIL("Error in FinalizeTemperatureFieldForComptonHeating.\n");
	  }	
    
    debug = debug_store;

    /* If we're using the HII restricted timestep, get the global
       maximum kph in I-fronts. */

    START_PERF();
    if (RadiativeTransferHIIRestrictedTimestep) {
      float LocalMaximumkph = -1e20;
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  LocalMaximumkph = max(LocalMaximumkph,
				Temp->GridData->ReturnMaximumkphIfront());
      LocalMaximumkph = CommunicationMaxValue(LocalMaximumkph);
      MetaData->GlobalMaximumkphIfront = LocalMaximumkph;
    }
    END_PERF(12);

#ifdef DEBUG
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->ErrorCheckPhotonNumber(lvl);
#endif

    if (!LoopTime)
      break;
    
    FirstTime = false;

  } // ENDWHILE GridTime >= PhotonTime

  /* Cleanup photon memory pool if we're deleting all photons between
     timesteps, i.e. no need to save photons */

  START_PERF();
#ifdef MEMORY_POOL
  const int PhotonMemorySize = MEMORY_POOL_SIZE;
  int PhotonSize = sizeof(PhotonPackageEntry);
  if (RadiativeTransferAdaptiveTimestep) {
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->DeletePhotonPackages(TRUE);
    delete PhotonMemoryPool;
    PhotonMemoryPool = new MPool::MemoryPool(PhotonMemorySize*PhotonSize,
					     PhotonSize,
					     PhotonMemorySize*PhotonSize/4);
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->InitializePhotonPackages();
  }
#endif
  END_PERF(13);

#ifdef REPORT_PERF
  if (!FirstTime) {
    if (debug) printf("EvolvePhotons: total time = %g\n", ReturnWallTime()-ep0);
    for (int i = 0; i < 1; i++) {  // Only report on ROOT, not all
      CommunicationBarrier();
      if (MyProcessorNumber == i) {

	printf("P%d:", MyProcessorNumber);
	fpcol(PerfCounter, 14, 14, stdout);
	fflush(stdout);
      }
    }
  }
#endif

  /* Delete GridList */

  ThinGridList *Orphan;
  TempGridList = RS_GridList;
  while (TempGridList != NULL) {
    Orphan = TempGridList;
    TempGridList = TempGridList->NextGrid;
    if (Orphan != NULL) delete Orphan;
  }

  /* Delete grid lists */

  for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    if (nGrids[lvl] > 0) delete [] Grids[lvl];

  LCAPERF_STOP("EvolvePhotons");
  return SUCCESS;

}
