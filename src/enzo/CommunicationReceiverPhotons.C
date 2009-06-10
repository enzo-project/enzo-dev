#define DEBUG 0
/***********************************************************************
/
/  COMMUNICATION ROUTINE: TRANSFER PHOTONS
/
/  written by: John H. Wise
/  date:       November, 2005
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "GroupPhotonList.h"
#include "PhotonCommunication.h"

#ifdef USE_MPI
static Eint32 PH_ListOfIndices[MAX_PH_RECEIVE_BUFFERS];
static MPI_Status PH_ListOfStatuses[MAX_PH_RECEIVE_BUFFERS];
#endif /* USE_MPI */

void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int FindSuperSource(PhotonPackageEntry **PP, int &LeafID, 
		    int SearchNewTree = TRUE);

int CommunicationReceiverPhotons(LevelHierarchyEntry *LevelArray[],
				 int &keep_transporting)
{

#ifdef USE_MPI

  /* Get any buffered receives */

  int NoErrorSoFar = TRUE;
  int FirstTime = TRUE;
  Eint32 ReceivesCompletedToDate = 0, NumberOfCompletedRequests, index, errcode;
  Eint32 TotalReceives = PH_CommunicationReceiveIndex;
  int TotalReceivedPhotons = 0;
  int *CompletedRequests = new int[TotalReceives];
  PhotonPackageEntry *NewPack = NULL;
  PhotonPackageEntry *ToPP = NULL;
  int lvl, gi, dim, i, count, NumberOfActiveRequests;
  grid *ToGrid;
  int ret;

  HierarchyEntry **Grids[MAX_DEPTH_OF_HIERARCHY];
  int nGrids[MAX_DEPTH_OF_HIERARCHY];
  int level;

  if (DEBUG)
    printf("P(%"ISYM") in PH_CRH with %"ISYM" requests\n", MyProcessorNumber,
	   PH_CommunicationReceiveIndex);

  for (i = 0; i < TotalReceives; i++)
    CompletedRequests[i] = FALSE;

  while (ReceivesCompletedToDate < TotalReceives) {

    /* Call the MPI wait handler */

    NumberOfCompletedRequests = 0;
    if (DEBUG) {
      printf("PH_CRH[%"ISYM"][a] : %"ISYM" %"ISYM" %"ISYM"\n", MyProcessorNumber, TotalReceives, 
	     ReceivesCompletedToDate, NumberOfCompletedRequests);
      printf("PH_CRH[%"ISYM"][a1]: %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n", MyProcessorNumber,
	     PH_CommunicationReceiveMPI_Request[0],
	     PH_CommunicationReceiveMPI_Request[1],
	     PH_CommunicationReceiveMPI_Request[2],
	     PH_CommunicationReceiveMPI_Request[3],
	     PH_CommunicationReceiveMPI_Request[4]);
      fflush(stdout);
    }
    MPI_Waitsome(TotalReceives, PH_CommunicationReceiveMPI_Request,
		 &NumberOfCompletedRequests, 
		 PH_ListOfIndices, PH_ListOfStatuses);

    if (DEBUG) {
      printf("PH_CRH[%"ISYM"][b]: %"ISYM" %"ISYM" %"ISYM" (%"ISYM" %"ISYM" %"ISYM")\n", MyProcessorNumber,
	     TotalReceives, ReceivesCompletedToDate, NumberOfCompletedRequests, 
	     PH_ListOfIndices[0], PH_ListOfIndices[1], PH_ListOfIndices[2]);
      fflush(stdout);
    }

    /* Loop over receive handles, looking for completed (i.e. null)
       requests. */

    GroupPhotonList *RecvBuffer = NULL;
    int irecv, NumberReceives, ToCount;

    for (irecv = 0; irecv < NumberOfCompletedRequests; irecv++) {

      index = PH_ListOfIndices[irecv];

      if (PH_ListOfStatuses[index].MPI_ERROR != 0) {
	fprintf(stderr, "MPI Error on processor %"ISYM". "
		"Error number %"ISYM" on request %"ISYM"\n",
		MyProcessorNumber, PH_ListOfStatuses[index].MPI_ERROR, index);
	fprintf(stdout, "P(%"ISYM") index %"ISYM" -- mpi error %"ISYM"\n", MyProcessorNumber,
		index, PH_ListOfStatuses[index].MPI_ERROR);
      }

      if (CompletedRequests[index] == TRUE)
	continue;

      if (DEBUG) {
	printf("PH_CRH[P%"ISYM"][%"ISYM"]: processing request %"ISYM" (addr %"ISYM")\n",
	       MyProcessorNumber, irecv, index, 
	       PH_CommunicationReceiveMPI_Request[index]);
	fflush(stdout);
      }

      /* Get grid lists */
  
      if (FirstTime == TRUE) {
	for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
	  if (LevelArray[level] != NULL)
	    nGrids[level] = 
	      GenerateGridArray(LevelArray, level, &Grids[level]);
	FirstTime = FALSE;
      }

      /* Locate received data */

      RecvBuffer = (GroupPhotonList *) PH_CommunicationReceiveBuffer[index];

      /* Count photons first ... in case, I figure out a way to
	 allocate the memory in a contiguous block (i.e. new
	 photon_t[NumberReceives]).  I believe this will be much
	 faster. */

      NumberReceives = 0;
      while (RecvBuffer[NumberReceives].ToLevel != BUFFER_END)
	NumberReceives++;
      TotalReceivedPhotons += NumberReceives;

      if (DEBUG)
	printf("CTPhR(P%"ISYM"): Received %"ISYM" photons\n", MyProcessorNumber, 
	       NumberReceives);

      /* Insert received photons in the photon list of the receiving
	 grid */

      for (i = 0; i < NumberReceives; i++) {
	  
	lvl	 = RecvBuffer[i].ToLevel;
	gi	 = RecvBuffer[i].ToGrid;
	ToGrid = Grids[lvl][gi]->GridData;
	ToPP	 = ToGrid->ReturnPhotonPackagePointer();

	NewPack = new PhotonPackageEntry;
	NewPack->Photons		= RecvBuffer[i].buffer.Photons;
	NewPack->Type			= RecvBuffer[i].buffer.Type;
	NewPack->Energy		= RecvBuffer[i].buffer.Energy;
	NewPack->EmissionTimeInterval = 
	  RecvBuffer[i].buffer.EmissionTimeInterval;
	NewPack->EmissionTime		= RecvBuffer[i].buffer.EmissionTime;
	NewPack->CurrentTime		= RecvBuffer[i].buffer.CurrentTime;
	NewPack->ColumnDensity	= RecvBuffer[i].buffer.ColumnDensity;
	NewPack->CrossSection		= RecvBuffer[i].buffer.CrossSection;
	NewPack->Radius		= RecvBuffer[i].buffer.Radius;
	NewPack->ipix			= RecvBuffer[i].buffer.ipix;
	NewPack->level		= RecvBuffer[i].buffer.level;

	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  NewPack->SourcePosition[dim] = 
	    RecvBuffer[i].buffer.SourcePosition[dim];

	NewPack->SourcePositionDiff = RecvBuffer[i].buffer.SourcePositionDiff;

	/* Search for the corresponding SuperSource, given a source ID
	   on the tree */

//	printf("CTPhR(P%"ISYM"): Photon %"ISYM" :: lvl %"ISYM", grid %"ISYM", srcid=%"ISYM", L = %"GSYM"\n",
//	       MyProcessorNumber, i, lvl, gi, 
//	       RecvBuffer[i].buffer.SuperSourceID, NewPack->Photons);

	if (RadiativeTransferSourceClustering) {
	  if (FindSuperSource(&NewPack, RecvBuffer[i].buffer.SuperSourceID) 
	      == FAIL) {
	    fprintf(stderr, "Error in FindSuperSource.\n");
	    ENZO_FAIL("");
	  }
	} else
	  NewPack->CurrentSource = NULL;

	InsertPhotonAfter(ToPP, NewPack);

	/* Update photon count */

	ToCount = ToGrid->ReturnNumberOfPhotonPackages();
	ToGrid->SetNumberOfPhotonPackages(ToCount+1);

	if (DEBUG > 1)
	  printf("CTPhR(P%"ISYM"): Photon %"ISYM" :: lvl %"ISYM", grid %"ISYM", srcid=%"ISYM", L = %"GSYM" (%"GSYM")\n",
		 MyProcessorNumber, i, lvl, gi, 
		 RecvBuffer[i].buffer.SuperSourceID, NewPack->Photons,
		 ToPP->NextPackage->Photons);

      } // ENDFOR transferred photons (i)

      delete [] RecvBuffer;
      CompletedRequests[index] = TRUE;
      //      if (PH_CommunicationReceiveMPI_Request[index] != MPI_REQUEST_NULL)
      //        MPI_Request_free(PH_CommunicationReceiveMPI_Request+index);
      //      PH_CommunicationReceiveMPI_Request[index] = NULL;
      ReceivesCompletedToDate++;

    } // ENDFOR completed requests (index)

  } // ENDWHILE receiving  

  PH_CommunicationReceiveIndex = 0;  

  if (FirstTime == FALSE)
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      if (LevelArray[level] != NULL)
	delete [] Grids[level];

  if (TotalReceives > 0)
    delete [] CompletedRequests;

  if (TotalReceivedPhotons > 0)
    keep_transporting = 1;
  
//  if (DEBUG)
//    printf("P(%"ISYM") out of PH_CRH with %"ISYM" requests\n", MyProcessorNumber,
//	   PH_CommunicationReceiveIndex);

#endif /* USE_MPI */

  return SUCCESS;

}
