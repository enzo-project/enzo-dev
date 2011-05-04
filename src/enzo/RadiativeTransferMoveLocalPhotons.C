#define DEBUG 0
/***********************************************************************
/
/  RAY TRACING ROUTINE: TRANSFER PHOTONS TO LOCAL GRIDS
/
/  written by: John H. Wise
/  date:       September, 2010
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

int RadiativeTransferMoveLocalPhotons(ListOfPhotonsToMove **AllPhotons,
				      int &keep_transporting)
{

  ListOfPhotonsToMove *Mover, *Destroyer, *Last;
  PhotonPackageEntry *ToGridPackages = NULL;
  int ToGridNumber, FromGridNumber;

  /* insert PhotonPackage in the correct grid if it's on the same
     processor */

  keep_transporting = FALSE;
  Last = *AllPhotons;
  Mover = (*AllPhotons)->NextPackageToMove;
  while (Mover != NULL) {  
    if (MyProcessorNumber == Mover->ToGrid->ReturnProcessorNumber()) {
      ToGridNumber = Mover->ToGrid->ReturnNumberOfPhotonPackages();
      FromGridNumber = Mover->FromGrid->ReturnNumberOfPhotonPackages();
      Mover->ToGrid->SetNumberOfPhotonPackages(ToGridNumber+1);
      Mover->FromGrid->SetNumberOfPhotonPackages(FromGridNumber-1);

      if (Mover->PausedPhoton)
	ToGridPackages = Mover->ToGrid->ReturnPausedPackagePointer();
      else
	ToGridPackages = Mover->ToGrid->ReturnPhotonPackagePointer();
      Mover->PhotonPackage->NextPackage = ToGridPackages->NextPackage;
      if (ToGridPackages->NextPackage != NULL) 
	ToGridPackages->NextPackage->PreviousPackage = Mover->PhotonPackage;
      ToGridPackages->NextPackage = Mover->PhotonPackage;
      Mover->PhotonPackage->PreviousPackage = ToGridPackages;

      if (Mover->PausedPhoton == FALSE)
	keep_transporting = TRUE;

      // Remove the photon from the move list and free memory
      Last->NextPackageToMove = Mover->NextPackageToMove;
      Destroyer = Mover;
      Mover = Mover->NextPackageToMove;
      delete Destroyer;

    } // ENDIF same processor
    else {
      Last = Mover;
      Mover = Mover->NextPackageToMove;                // next one
    }

  } // ENDWHILE Mover != NULL 

  return SUCCESS;

}
