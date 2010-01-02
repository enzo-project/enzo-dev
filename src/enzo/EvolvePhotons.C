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

/* function prototypes */
void my_exit(int status);
int CommunicationTransferPhotons(LevelHierarchyEntry *LevelArray[], 
				 ListOfPhotonsToMove **AllPhotons,
				 int &keep_transporting);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
int InitiateKeepTransportingCheck(int keep_transporting);
int StopKeepTransportingCheck();
int InitializePhotonCommunication();
int KeepTransportingCheck(int &keep_transporting);
RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS);
PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
int CreateSourceClusteringTree(int nShine, SuperSourceData *SourceList,
			       LevelHierarchyEntry *LevelArray[]);
void PrintSourceClusteringTree(SuperSourceEntry *leaf);
int CommunicationSyncNumberOfPhotons(LevelHierarchyEntry *LevelArray[]);
int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove,
				     int level);

/* EvolvePhotons function */
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		  Star *AllStars, FLOAT GridTime, int level, int LoopTime)
{

  bool FirstTime = true;

  if (!RadiativeTransfer)
    return SUCCESS;

  /* Only call on the finest level */

  if (LevelArray[level+1] != NULL)
    return SUCCESS;

//  printf("GridTime = %f, PhotonTime = %f, dtPhoton = %g (Loop = %d)\n",
//	 GridTime, PhotonTime, dtPhoton, (GridTime >= PhotonTime));

  //if (dtPhoton < 0)
  //  return SUCCESS;

  //while (GridTime >= PhotonTime) {
  while (GridTime > PhotonTime) {

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
      
    /* Declarations */

    grid *Helper;
    int RefinementFactors[MAX_DIMENSION];

    /* Create an array (Grids) of all the grids. */

    typedef HierarchyEntry* HierarchyEntryPointer;
    HierarchyEntry **Grids;
    HierarchyEntry **Parents;
    LevelHierarchyEntry *Temp;
    int GridNum = 0, value, i, proc, lvl;
    int NumberOfGrids = 0;  

    // delete source if we are passed (or before) their lifetime (only
    // if not restarting)
    RadiationSourceEntry *RS;
    RS = GlobalRadiationSources->NextSource;
    int NumberOfSources = 0;
    while (RS != NULL) {
      if ( ((RS->CreationTime + RS->LifeTime) < PhotonTime ||
	    (RS->CreationTime > PhotonTime + dtPhoton)) &&
	   LoopTime == TRUE) {
	if (debug) {
	  fprintf(stdout, "\nEvolvePhotons: Deleted Source on lifetime limit \n");
	  fprintf(stdout, "EvolvePhotons:  %"GSYM" %"GSYM" %"GSYM" \n",
		  RS->CreationTime, RS->LifeTime, PhotonTime);
	}
	RS = DeleteRadiationSource(RS);
      } else {
	NumberOfSources++;                 // count sources
	RS = RS->NextSource;
      }
    }

    if (debug) fprintf(stdout, "%"ISYM" SRC(s)\n", NumberOfSources);

    int Rank, Dims[MAX_DIMENSION];
    FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];
  
    /* Initialize radiation fields */

    for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) 
	if (Temp->GridData->InitializeRadiativeTransferFields() == FAIL) {
	  fprintf(stderr, "Error in InitializeRadiativeTransferFields.\n");
	  ENZO_FAIL("");
	}

    /* create temperature fields for Compton heating */  

    if (RadiationXRayComptonHeating)  
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) 
	  if (Temp->GridData->InitializeTemperatureFieldForComptonHeating() == FAIL) {  
	    fprintf(stderr, "Error in InitializeTemperatureFieldForComptonHeating.\n");
	    ENZO_FAIL("");
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

    if (RadiativeTransferSourceClustering == TRUE) {
      if (CreateSourceClusteringTree(NULL, NULL, LevelArray) == FAIL) {
	fprintf(stderr, "Error in CreateSourceClusteringTree.\n");
	ENZO_FAIL("");
      }
      //PrintSourceClusteringTree(SourceClusteringTree);
    }

    // first identify sources and let them radiate 
    RS = GlobalRadiationSources->NextSource;
 
    while (RS != NULL) {
      int Continue = 1;
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; (lvl >= 0 && Continue); lvl--) {
	for (Temp = LevelArray[lvl]; (Temp && Continue); 
	     Temp = Temp->NextGridThisLevel)
	  if (MyProcessorNumber == Temp->GridData->ReturnProcessorNumber())
	    if (Temp->GridData->PointInGrid(RS->Position)) {
	      if (Temp->GridData->Shine(RS) == FAIL) {
		fprintf(stderr, "Error in Shine.\n");
		ENZO_FAIL("");
	      }
	      Continue = FALSE; // do not continue with this source
	    } // If source in grid

	/* For MPI, communicate to minimum value of Continue to ensure
	   that the source's host grid was found. */

	Continue = CommunicationMinValue(Continue);

	if (lvl == 0 && Continue) {  // this should never happen ... 
	  fprintf(stderr, "Could not find grid for source %x: "
		  "Pos: %"FSYM" %"FSYM" %"FSYM"\n",
		  RS, RS->Position[0], RS->Position[1], RS->Position[2]);
	  ENZO_FAIL("");
	}
      }    // Loop through levels 
      RS = RS->NextSource;
    }    // while still sources 

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
    //	  fprintf(stderr, "Error in grid->AllocateInterpolatedRadiation.\n");
    //	  ENZO_FAIL("");
    //	}

    /* Evolve all photons by fixed timestep. */
  
    ListOfPhotonsToMove *PhotonsToMove = new ListOfPhotonsToMove;
    PhotonsToMove->NextPackageToMove = NULL;

    int keep_transporting = 1;
    int ThisProcessor;

    HierarchyEntry **Temp0;
    int nGrids0 = GenerateGridArray(LevelArray, 0, &Temp0);
    grid **Grids0 = new grid*[nGrids0];
    for (i = 0; i < nGrids0; i++)
      Grids0[i] = Temp0[i]->GridData;

    /* Initialize nonblocking MPI routine */

    InitializePhotonCommunication();

    /* Transport the rays! */

    while (keep_transporting) {
      keep_transporting = 0;
      PhotonsToMove->NextPackageToMove = NULL;
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--) {

	NumberOfGrids = GenerateGridArray(LevelArray, lvl, &Grids);
	for (GridNum = 0; GridNum < NumberOfGrids; GridNum++) {

	  if (Grids[GridNum]->ParentGrid != NULL) 
	    Helper = Grids[GridNum]->ParentGrid->GridData;
	  else
	    Helper = NULL;
	  if (Grids[GridNum]->GridData->TransportPhotonPackages
	      (lvl, &PhotonsToMove, GridNum, Grids0, 
	       nGrids0, Helper, Grids[GridNum]->GridData) == FAIL) {
	    fprintf(stderr, "Error in %"ISYM" th grid. "
		    "grid->TransportPhotonPackages.\n",GridNum);
	    ENZO_FAIL("");
	  }

	} // ENDFOR grids

	delete [] Grids;

      }                          // loop over levels

      if (PhotonsToMove->NextPackageToMove != NULL)
	keep_transporting = 1;

//    printf("EvoPH[P%"ISYM"]: keep_transporting = %"ISYM", PhotonsToMove = %x\n",
//	   MyProcessorNumber, keep_transporting, PhotonsToMove->NextPackageToMove);

      /* Check if there are any photons leaving this grid.  If so, move them. */
      
      if (CommunicationTransferPhotons(LevelArray, &PhotonsToMove, 
				       keep_transporting) == FAIL) {
	fprintf(stderr, "Error in CommunicationTransferPhotons.\n");
	ENZO_FAIL("");
      }

      /* Receive keep_transporting messages and take the MAX */

#ifdef NONBLOCKING
      InitiateKeepTransportingCheck(keep_transporting);
      KeepTransportingCheck(keep_transporting);
#else /* NON_BLOCKING */
      keep_transporting = CommunicationMaxValue(keep_transporting);
#endif

    }                           //  end while keep_transporting

    //  StopKeepTransportingCheck();

    /* Move all finished photon packages back to their original place,
       PhotonPackages.  For the adaptive timestep, we don't carryover
       any photons to the next timestep. */

    if (RadiativeTransferAdaptiveTimestep)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  Temp->GridData->DeletePhotonPackages();
    else
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  Temp->GridData->MoveFinishedPhotonsBack();

    /* If we're keeping track of photon escape fractions on multiple
       processors, collect photon counts from all processors */

    FILE *fptr;

    if (RadiativeTransferPhotonEscapeRadius > 0) {
      CommunicationSumValues(EscapedPhotonCount, 4);
      if (MyProcessorNumber == ROOT_PROCESSOR) {

	/* Open f_esc file for writing */

	if (TotalEscapedPhotonCount[0] <= 0) {
	  if ((fptr = fopen(PhotonEscapeFilename, "w")) == NULL) {
	    fprintf(stderr, "Error opening file %s\n", PhotonEscapeFilename);
	    ENZO_FAIL("");
	  }
	  fprintf(fptr, 
		  "# Time TotalPhotons fesc(0.5rvir) fesc(rvir) fesc(2rvir)\n");
	} else {
	  if ((fptr = fopen(PhotonEscapeFilename, "a")) == NULL) {
	    fprintf(stderr, "Error opening file %s\n", PhotonEscapeFilename);
	    ENZO_FAIL("");
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

    /* Coupled rate & energy solver */

    float dtMin, dtGrid;
    float dtLastLevel = 1e20;

    int debug_store = debug;
    debug = FALSE;

    // Divide the photo-ionization and photo-heating rates by the
    // number of particles (rho * dx^3)
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	if (Temp->GridData->RadiationPresent() == TRUE)
	  if (Temp->GridData->FinalizeRadiationFields() == FAIL) {
	    fprintf(stderr, "Error in FinalizeRadiationFields.\n");
	    ENZO_FAIL("");
	  }

    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	if (Temp->GridData->RadiationPresent() == TRUE) {

	  if (RadiativeTransferCoupledRateSolver && 
	      RadiativeTransferOpticallyThinH2)
	    if (Temp->GridData->AddH2Dissociation(AllStars) == FAIL) {
	      fprintf(stderr, "Error in AddH2Dissociation.\n");
	      ENZO_FAIL("");
	    }

	  if (RadiativeTransferCoupledRateSolver)
	    if (Temp->GridData->SolveCoupledRateEquations() == FAIL) {
	      fprintf(stderr, "Error in grid->SolveCoupledRateEquations.\n");
	      ENZO_FAIL("");
	    }

	  if (RadiativeTransferCoupledRateSolver &&
	      RadiativeTransferInterpolateField)
	    Temp->GridData->DeleteInterpolatedFields();

	} /* ENDIF radiation */

    /* For the non-coupled (i.e. cells without radiation) rate & energy
       solver, we have to set the H2 dissociation rates */

    if (RadiativeTransferOpticallyThinH2)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->RadiationPresent() == FALSE)
	    if (Temp->GridData->AddH2Dissociation(AllStars) == FAIL) {
	      fprintf(stderr, "Error in AddH2Dissociation.\n");
	      ENZO_FAIL("");
	    }

    /* Clean up temperature field */

    if (RadiationXRayComptonHeating)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->FinalizeTemperatureFieldForComptonHeating() == FAIL) {  
	    fprintf(stderr, "Error in FinalizeTemperatureFieldForComptonHeating.\n");
	    ENZO_FAIL("");
	  }	
    
    debug = debug_store;

    /* We don't rely on the count NumberOfPhotonPackages here, so they
       aren't synchronized across processors.  But in RebuildHierarchy, 
       this number is needed.  Synchronize them now. */

    CommunicationSyncNumberOfPhotons(LevelArray);

    /* If we're using the HII restricted timestep, get the global
       maximum kph in I-fronts. */

    if (RadiativeTransferHIIRestrictedTimestep) {
      float LocalMaximumkph = -1e20;
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  LocalMaximumkph = max(LocalMaximumkph,
				Temp->GridData->ReturnMaximumkphIfront());
      LocalMaximumkph = CommunicationMaxValue(LocalMaximumkph);
      MetaData->GlobalMaximumkphIfront = LocalMaximumkph;
    }

#ifdef DEBUG
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->ErrorCheckPhotonNumber(lvl);
#endif

    if (!LoopTime)
      break;
    
    FirstTime = false;

  } // ENDWHILE GridTime >= PhotonTime

  return SUCCESS;

}
