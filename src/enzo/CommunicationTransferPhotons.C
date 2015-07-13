#define DEBUG 0
/***********************************************************************
/
/  COMMUNICATION ROUTINE: TRANSFER PHOTONS
/
/  written by: John H. Wise
/  date:       November, 2005
/  modified1:  JHW (December 2007) -- Triple-phase communication
/
/  PURPOSE:
/
/  TODO: Stop using MPI_BYTE for communication, and create an MPI
/        structure
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

PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int CommunicationReceiverPhotons(LevelHierarchyEntry *LevelArray[],
				 bool local_transport,
				 int &keep_transporting);
int CommunicationNumberOfPhotonSends(int *nPhoton, int size);
//int InitiatePhotonNumberSend(int *nPhoton);
//int InitializePhotonReceive(int group_size);
void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);
#ifdef USE_MPI
int InitializePhotonReceive(int max_size, bool local_transport, 
			    MPI_Datatype MPI_PhotonType);
int CommunicationBufferPurge(void);
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
                              int Target, int Tag, MPI_Comm CommWorld, 
			      int BufferSize);

static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_PhotonList;

#endif /* USE_MPI */

int CommunicationTransferPhotons(LevelHierarchyEntry *LevelArray[], 
				 ListOfPhotonsToMove **AllPhotons,
				 char *kt_global,
				 int &keep_transporting)
{

  ListOfPhotonsToMove *Mover, *Destroyer;
  PhotonPackageEntry *ToGridPackages = NULL;
  int ToGridNumber, FromGridNumber;

  /* Serial case */

  if (NumberOfProcessors == 1) {
    Mover = (*AllPhotons)->NextPackageToMove;
    while (Mover != NULL) {  // insert PhotonPackage in the correct grid

      ToGridNumber = Mover->ToGrid->ReturnNumberOfPhotonPackages();
      FromGridNumber = Mover->FromGrid->ReturnNumberOfPhotonPackages();
      Mover->ToGrid->SetNumberOfPhotonPackages(ToGridNumber+1);
      Mover->FromGrid->SetNumberOfPhotonPackages(FromGridNumber-1);

      if (Mover->PausedPhoton == FALSE)
	ToGridPackages = Mover->ToGrid->ReturnPhotonPackagePointer();
      else
	ToGridPackages = Mover->ToGrid->ReturnPausedPackagePointer();
      Mover->PhotonPackage->NextPackage = ToGridPackages->NextPackage;
      if (ToGridPackages->NextPackage != NULL) 
	ToGridPackages->NextPackage->PreviousPackage = Mover->PhotonPackage;
      ToGridPackages->NextPackage = Mover->PhotonPackage;
      Mover->PhotonPackage->PreviousPackage = ToGridPackages;
      Mover = Mover->NextPackageToMove;                // next one

    } // end      while Mover != Null 

    Mover = (*AllPhotons)->NextPackageToMove;
    while (Mover != NULL) {  // free memory
      Destroyer = Mover;
      Mover = Mover->NextPackageToMove;                // next one	
      delete Destroyer;
    }

    return SUCCESS;

  }  /* ENDIF serial */

#ifdef USE_MPI

  MPI_Status status;
  
  /* Generate a new MPI type corresponding to the PhotonList struct. */
  
  if (FirstTimeCalled) {
    MPI_Type_contiguous(sizeof(GroupPhotonList), MPI_BYTE, &MPI_PhotonList);
    MPI_Type_commit(&MPI_PhotonList);
    FirstTimeCalled = FALSE;
  }

  /* If parallel, Partition photons into linked lists that are
     transferred to the same grid */

  float value;
  int ivalue, GridNum, level, i, dim, proc;
  int NumberOfGrids = 0, count = 0;

  GroupPhotonList **SendList = new GroupPhotonList*[NumberOfProcessors];
  int *PhotonCounter = new int[NumberOfProcessors];
  int *nPhoton = new int[NumberOfProcessors];

  /* Initialize counters */

  for (proc = 0; proc < NumberOfProcessors; proc++) {
    nPhoton[proc] = 0;
    PhotonCounter[proc] = 0;
  }

  /* Count photons to move */

  int NumberToMove = 0;
  ListOfPhotonsToMove *TempList = (*AllPhotons)->NextPackageToMove;
  while (TempList != NULL) {
    nPhoton[TempList->ToProcessor]++;
    NumberToMove++;
    TempList = TempList->NextPackageToMove;
  }
  nPhoton[MyProcessorNumber] = 0;

  /* Last entry in comm. lists will be a marker to indicate the end.
     Adjust nPhoton for this. */

  for (proc = 0; proc < NumberOfProcessors; proc++)
    if (nPhoton[proc] > 0)
      nPhoton[proc]++;

  /* Communicate number of photons to move */

  //InitiatePhotonNumberSend(nPhoton);
  CommunicationNumberOfPhotonSends(nPhoton, PHOTON_BUFFER_SIZE);


  /* Allocate memory */

  for (proc = 0; proc < NumberOfProcessors; proc++)
    if (proc != MyProcessorNumber) 
      SendList[proc] = new GroupPhotonList[nPhoton[proc]];
    else
      SendList[proc] = NULL;

  /* Collect photons into lists groups by ToProcessor */

  ListOfPhotonsToMove *LastMover;
  PhotonPackageEntry *dummy;
  int ToProc, ToCount, TempLevel, TempGridNum, FromNumber;
  int localCounter = 0;

  LastMover = *AllPhotons;
  Mover = (*AllPhotons)->NextPackageToMove;
  while (Mover != NULL) {

    ToProc = Mover->ToProcessor;
    ToCount = PhotonCounter[ToProc];

    /* If the photons are being transported to a grid on the same
       processor, transfer them as in the serial case */

    if (ToProc == MyProcessorNumber) {

      ToGridNumber = Mover->ToGrid->ReturnNumberOfPhotonPackages();
      FromGridNumber = Mover->FromGrid->ReturnNumberOfPhotonPackages();
      Mover->ToGrid->SetNumberOfPhotonPackages(ToGridNumber+1);
      Mover->FromGrid->SetNumberOfPhotonPackages(FromGridNumber-1);

      if (Mover->PausedPhoton == FALSE)
	ToGridPackages = Mover->ToGrid->ReturnPhotonPackagePointer();
      else
	ToGridPackages = Mover->ToGrid->ReturnPausedPackagePointer();
      Mover->PhotonPackage->NextPackage = ToGridPackages->NextPackage;
      if (ToGridPackages->NextPackage != NULL) 
	ToGridPackages->NextPackage->PreviousPackage = Mover->PhotonPackage;
      ToGridPackages->NextPackage = Mover->PhotonPackage;
      Mover->PhotonPackage->PreviousPackage = ToGridPackages;

      localCounter++;

    } /* ENDIF same processor */

    else {

      TempLevel = Mover->ToLevel;
      TempGridNum = Mover->ToGridNum;

      /* Add the photon information to the send list */

      SendList[ToProc][ToCount].ToLevel	       = TempLevel;
      SendList[ToProc][ToCount].ToGrid	       = TempGridNum;
      SendList[ToProc][ToCount].PausedPhoton   = Mover->PausedPhoton;
      SendList[ToProc][ToCount].buffer.Photons = Mover->PhotonPackage->Photons;
      SendList[ToProc][ToCount].buffer.Type    = Mover->PhotonPackage->Type;
      SendList[ToProc][ToCount].buffer.Energy  = Mover->PhotonPackage->Energy;

      SendList[ToProc][ToCount].buffer.EmissionTimeInterval = 
	Mover->PhotonPackage->EmissionTimeInterval;
      SendList[ToProc][ToCount].buffer.EmissionTime	    = 
	Mover->PhotonPackage->EmissionTime;
      SendList[ToProc][ToCount].buffer.CurrentTime	    = 
	Mover->PhotonPackage->CurrentTime;
      SendList[ToProc][ToCount].buffer.ColumnDensity        = 
	Mover->PhotonPackage->ColumnDensity;
      SendList[ToProc][ToCount].buffer.CrossSection         = 
	Mover->PhotonPackage->CrossSection;

      SendList[ToProc][ToCount].buffer.Radius = Mover->PhotonPackage->Radius;
      SendList[ToProc][ToCount].buffer.ipix   = Mover->PhotonPackage->ipix;
      SendList[ToProc][ToCount].buffer.level  = Mover->PhotonPackage->level;

      for (dim = 0; dim < MAX_DIMENSION; dim++)
	SendList[ToProc][ToCount].buffer.SourcePosition[dim] = 
	  Mover->PhotonPackage->SourcePosition[dim];
      SendList[ToProc][ToCount].buffer.SourcePositionDiff = 
	Mover->PhotonPackage->SourcePositionDiff;
      if (Mover->PhotonPackage->CurrentSource != NULL)
	SendList[ToProc][ToCount].buffer.SuperSourceID =
	  Mover->PhotonPackage->CurrentSource->LeafID;
      else
	SendList[ToProc][ToCount].buffer.SuperSourceID = -1;

      PhotonCounter[ToProc]++;

      /* Update photon count */

      FromNumber = Mover->FromGrid->ReturnNumberOfPhotonPackages();
      Mover->FromGrid->SetNumberOfPhotonPackages(FromNumber-1);

//      printf("CTPh(P%"ISYM"): Photon %"ISYM" (=>P%"ISYM") :: lvl %"ISYM", grid %"ISYM", srcid=%"ISYM", L = %"GSYM"\n",
//	     MyProcessorNumber, ToCount, ToProc, TempLevel, TempGridNum,
//	     SendList[ToProc][ToCount].buffer.SuperSourceID,
//	     SendList[ToProc][ToCount].buffer.Photons);
      if (DEBUG > 1)
	printf("CTPh(P%"ISYM"): Photon %"ISYM" (=>P%"ISYM") :: lvl %"ISYM", grid %"ISYM", L = %"GSYM"\n",
	       MyProcessorNumber, ToCount, ToProc, TempLevel, TempGridNum,
	       SendList[ToProc][ToCount].buffer.Photons);

    }  /* ENDELSE different processor */

    LastMover = Mover;
    Mover = Mover->NextPackageToMove;

  } /* ENDWHILE Mover */

  /* Mark the last buffer entry so we can count photons when we
     receive. */

  for (proc = 0; proc < NumberOfProcessors; proc++)
    if (nPhoton[proc] > 0)
      SendList[proc][PhotonCounter[proc]++].ToLevel = BUFFER_END;

  if (DEBUG) {
    for (proc = 0; proc < NumberOfProcessors; proc++)
      if (PhotonCounter[proc])
	printf("CTPh(P%"ISYM"): Packed %"ISYM" photons into COMM buffer #%"ISYM".\n", 
	       MyProcessorNumber, PhotonCounter[proc], proc);
    if (localCounter)
      printf("CTPh(P%"ISYM"): Transferred %"ISYM" photons locally.\n", 
	     MyProcessorNumber, localCounter);  
  }

  /* Now that we're done reading the photon move list, delete it */

  Mover = (*AllPhotons)->NextPackageToMove;
  while (Mover != NULL) {
    Destroyer = Mover;
    dummy = Mover->PhotonPackage;

    // Only delete the photon if it's being transferred to another processor
    if (Mover->ToProcessor != MyProcessorNumber)
      delete dummy;

    Mover = Mover->NextPackageToMove;                // next one
    delete Destroyer;
  }

  /***************************************************************/
  /*                  TRIPLE-PHASE COMMUNICATION                 */
  /* Modeled after the scheme in communication.h and EvolveLevel */
  /***************************************************************/

  bool local_transport;
  int NumberOfMessages, Offset;
  Eint32 tag, SizeOfGroupPhotonList, Size;
  
  local_transport = (localCounter > 0);
  MPI_Type_size(MPI_PhotonList, &SizeOfGroupPhotonList);
  //SizeOfGroupPhotonList = sizeof(GroupPhotonList);

  /* First stage: check for any received nPhoton messages.  For the
     completed messages, post receive calls from all processors with
     photons to receive */

#ifdef NONBLOCKING_RT
  InitializePhotonReceive(PHOTON_BUFFER_SIZE, true, MPI_PhotonList);
#else
  InitializePhotonReceive(PHOTON_BUFFER_SIZE, false, MPI_PhotonList);
#endif

  /* Second stage: post sends to processors */
     
  for (proc = 0; proc < NumberOfProcessors; proc++)
    if (proc != MyProcessorNumber) {
      NumberOfMessages = nPhoton[proc] / PHOTON_BUFFER_SIZE;
      if (nPhoton[proc] % PHOTON_BUFFER_SIZE > 0) NumberOfMessages++;
      //tag = MPI_PHOTONGROUP_TAG*10 + nPhoton[proc];
      tag = MPI_PHOTONGROUP_TAG;
      if (DEBUG && NumberOfMessages > 0)
	printf("CTPh(P%"ISYM"): Sending %"ISYM" photons to P%"ISYM" (%d messages, TAG=%d)\n", 
	       MyProcessorNumber, nPhoton[proc], proc, NumberOfMessages, tag);
      for (i = 0; i < NumberOfMessages; i++) {
	Offset = i*PHOTON_BUFFER_SIZE;
	Size = (i < NumberOfMessages-1) ? PHOTON_BUFFER_SIZE : (nPhoton[proc]-Offset);
	CommunicationBufferedSend(SendList[proc]+Offset, 
				  Size, MPI_PhotonList, proc, tag, 
				  MPI_COMM_WORLD,
				  Size*sizeof(GroupPhotonList));
      } // ENDFOR messages
      delete [] SendList[proc];
    } // ENDIF other processor

  /* Third stage: finally receive the data and transfer them to their
     respective grids  */

#ifndef NONBLOCKING_RT
  local_transport = false;
#endif
  CommunicationReceiverPhotons(LevelArray, local_transport, 
			       keep_transporting);

  if (kt_global != NULL)
    for (proc = 0; proc < NumberOfProcessors; proc++)
      if (proc != MyProcessorNumber && nPhoton[proc] > 0 &&
	  kt_global[proc] != HALT_TRANSPORT) {
	kt_global[proc] = SENT_DATA;
	keep_transporting = 1;
      }

  /* Clean up */

  CommunicationBufferPurge();
  delete [] SendList;
  delete [] nPhoton;
  delete [] PhotonCounter;

#endif /* USE_MPI */

  return SUCCESS;

}
