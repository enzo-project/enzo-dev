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
void CommunicationCheckForErrors(int NumberOfStatuses, MPI_Status *statuses,
				 char *msg=NULL);
int CommunicationFindOpenRequest(MPI_Request *requests, Eint32 last_free,
				 Eint32 nrequests, Eint32 index, 
				 Eint32 &max_index);
#endif /* USE_MPI */

void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int FindSuperSource(PhotonPackageEntry **PP, int &LeafID, 
		    int SearchNewTree = TRUE);

#define NO_DEBUG_CRP
#define NO_DEBUG_CRP2

int CommunicationReceiverPhotons(LevelHierarchyEntry *LevelArray[],
				 bool local_transport,
				 int &keep_transporting)
{

#ifdef USE_MPI

  /* Get any buffered receives */

  int NoErrorSoFar = TRUE;
  int FirstTime = TRUE;
  MPI_Arg ReceivesCompletedToDate = 0, NumberOfCompletedRequests, index, errcode;
  MPI_Arg TotalReceives = PH_CommunicationReceiveMaxIndex;
  int TotalReceivedPhotons = 0;
  bool *CompletedRequests = NULL;
  PhotonPackageEntry *NewPack, *ToPP;
  int lvl, gi, dim, i, count, NumberOfActiveRequests;
  grid *ToGrid;
  int ret, level;

  HierarchyEntry **Grids[MAX_DEPTH_OF_HIERARCHY];
  int nGrids[MAX_DEPTH_OF_HIERARCHY];

  if (TotalReceives == 0)
    return SUCCESS;

#ifdef DEBUG_CRP
  printf("P(%"ISYM") in PH_CRH with %"ISYM" requests (local=%d)\n", 
	 MyProcessorNumber, TotalReceives, local_transport);
#endif

  CompletedRequests = new bool[TotalReceives];
  for (i = 0; i < TotalReceives; i++)
    CompletedRequests[i] = false;

#ifndef NONBLOCKING_RT
  if (TotalReceives > 0)
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      if (LevelArray[level] != NULL)
	nGrids[level] = GenerateGridArray(LevelArray, level, &Grids[level]);
  while (ReceivesCompletedToDate < TotalReceives) {
#endif

  NumberOfCompletedRequests = 0;
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  /* Wait for >1 receives */

  if (local_transport)
    MPI_Testsome(TotalReceives, PH_CommunicationReceiveMPI_Request,
		 &NumberOfCompletedRequests, 
		 PH_ListOfIndices, PH_ListOfStatuses);
  else
    MPI_Waitsome(TotalReceives, PH_CommunicationReceiveMPI_Request,
		 &NumberOfCompletedRequests, 
		 PH_ListOfIndices, PH_ListOfStatuses);
  
#ifdef DEBUG_CRP
  printf("PH_CRH[%"ISYM"][b]: %"ISYM" %"ISYM" %"ISYM" (%"ISYM" %"ISYM" %"ISYM")\n", 
	 MyProcessorNumber, TotalReceives, ReceivesCompletedToDate, 
	 NumberOfCompletedRequests, 
	 PH_ListOfIndices[0], PH_ListOfIndices[1], PH_ListOfIndices[2]);
  fflush(stdout);
#endif

  if (NumberOfCompletedRequests > 0)
    CommunicationCheckForErrors(TotalReceives, PH_ListOfStatuses,
				"CommunicationReceiverPhotons");
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);

  /* Get grid lists */

#ifdef NONBLOCKING_RT
  if (NumberOfCompletedRequests > 0) {
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      if (LevelArray[level] != NULL)
	nGrids[level] = GenerateGridArray(LevelArray, level, &Grids[level]);
  }
#endif

  /* Loop over receive handles, looking for completed (i.e. null)
     requests. */

  GroupPhotonList *RecvBuffer = NULL;
  int irecv, NumberReceives, ToCount;

  for (irecv = 0; irecv < NumberOfCompletedRequests; irecv++) {

    index = PH_ListOfIndices[irecv];
    if (CompletedRequests[index])
      continue;

#ifdef DEBUG_CRP
    printf("PH_CRH[P%"ISYM"][%"ISYM"]: processing request %"ISYM"\n",
	   MyProcessorNumber, irecv, index);
    fflush(stdout);
#endif

    /* Locate received data */

    RecvBuffer = (GroupPhotonList *) PH_CommunicationReceiveBuffer[index];

    /* Count photons first ... in case, I figure out a way to allocate
       the memory in a contiguous block (i.e. new
       photon_t[NumberReceives]).  I believe this will be much
       faster. */

    NumberReceives = 0;
    while (RecvBuffer[NumberReceives].ToLevel != BUFFER_END &&
	   NumberReceives < PHOTON_BUFFER_SIZE)
      NumberReceives++;
    TotalReceivedPhotons += NumberReceives;

#ifdef DEBUG_CRP
    printf("CTPhR(P%"ISYM"): Received %"ISYM" photons\n", MyProcessorNumber, 
	   NumberReceives);
#endif

    /* Insert received photons in the photon list of the receiving
       grid */

    for (i = 0; i < NumberReceives; i++) {
	  
      lvl	 = RecvBuffer[i].ToLevel;
      gi	 = RecvBuffer[i].ToGrid;

      // Double check if the grids exists on this processor and if the
      // grid number is valid.  If not, skip and warn the user.
      if (gi >= nGrids[lvl]) {
	printf("P%d: WARNING: CommunicationReceiverPhotons: Bad grid number = %d\n"
	       "\t Receive %d, level %d, NumberOfGrids = %d.  SKIPPING!\n",
	       MyProcessorNumber, gi, i, lvl, nGrids[lvl]);
	continue;
      }
      else if (Grids[lvl][gi]->GridData->ReturnProcessorNumber() != 
	       MyProcessorNumber) {
	printf("P%d: WARNING: CommunicationReceiverPhotons: This grid isn't on this processor!\n"
	       "\t Grid %d (P%d), Receive %d, level %d, NumberOfGrids = %d. SKIPPING!\n",
	       MyProcessorNumber, gi, Grids[lvl][gi]->GridData->ReturnProcessorNumber(), 
	       i, lvl, nGrids[lvl]);
	continue;
      }
      ToGrid = Grids[lvl][gi]->GridData;
      if (RecvBuffer[i].PausedPhoton == FALSE)
	ToPP	 = ToGrid->ReturnPhotonPackagePointer();
      else
	ToPP	 = ToGrid->ReturnPausedPackagePointer();

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

#ifdef DEBUG_CRP2
      printf("CTPhR(P%"ISYM"): Photon %"ISYM" :: lvl %"ISYM", grid %"ISYM
	     ", srcid=%"ISYM", L = %"GSYM"\n",
	     MyProcessorNumber, i, lvl, gi, 
	     RecvBuffer[i].buffer.SuperSourceID, NewPack->Photons);
#endif

      if (RadiativeTransferSourceClustering)
	FindSuperSource(&NewPack, RecvBuffer[i].buffer.SuperSourceID);
      else
	NewPack->CurrentSource = NULL;

      InsertPhotonAfter(ToPP, NewPack);

      /* Update photon count */

      ToCount = ToGrid->ReturnNumberOfPhotonPackages();
      ToGrid->SetNumberOfPhotonPackages(ToCount+1);

#ifdef DEBUG_CRP2
      printf("CTPhR(P%"ISYM"): Photon %"ISYM" :: lvl %"ISYM", grid %"ISYM
	     ", srcid=%"ISYM", L = %"GSYM" (%"GSYM")\n",
	     MyProcessorNumber, i, lvl, gi, 
	     RecvBuffer[i].buffer.SuperSourceID, NewPack->Photons,
	     ToPP->NextPackage->Photons);
#endif

    } // ENDFOR transferred photons (i)

    delete [] RecvBuffer;
    CompletedRequests[index] = true;
    //      if (PH_CommunicationReceiveMPI_Request[index] != MPI_REQUEST_NULL)
    //        MPI_Request_free(PH_CommunicationReceiveMPI_Request+index);
    //      PH_CommunicationReceiveMPI_Request[index] = NULL;
    ReceivesCompletedToDate++;

  } // ENDFOR completed requests (index)

#ifndef NONBLOCKING_RT
  } // ENDWHILE receiving
  PH_CommunicationReceiveIndex = 0;
  PH_CommunicationReceiveMaxIndex = 0;
#else
  PH_CommunicationReceiveIndex = 
    CommunicationFindOpenRequest(PH_CommunicationReceiveMPI_Request, NO_HINT,
				 MAX_PH_RECEIVE_BUFFERS,
				 PH_CommunicationReceiveIndex,
				 PH_CommunicationReceiveMaxIndex);
#endif

#ifdef NONBLOCKING_RT
  if (NumberOfCompletedRequests > 0)
#else
  if (TotalReceives > 0)
#endif
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      if (LevelArray[level] != NULL)
	delete [] Grids[level];

  delete [] CompletedRequests;
  if (TotalReceivedPhotons > 0)
    keep_transporting = 1;
  
#ifdef DEBUG_CRP
  printf("P(%"ISYM") out of PH_CRH with %"ISYM" requests. nphotons=%d, kt=%d\n",
	 MyProcessorNumber, PH_CommunicationReceiveMaxIndex, 
	 TotalReceivedPhotons, keep_transporting);
#endif

#endif /* USE_MPI */

  return SUCCESS;

}
