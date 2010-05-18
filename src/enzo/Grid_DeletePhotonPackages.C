/***********************************************************************
/
/  GRID CLASS (DELETES ALL PHOTONS)
/
/  written by: John H. Wise
/  date:       November, 2005
/  modified1:
/
/  PURPOSE:  
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);

int grid::DeletePhotonPackages(int DeleteHeadPointer) {

  PhotonPackageEntry *PP;

  for (PP = PhotonPackages->NextPackage; PP; PP = PP->NextPackage)
    PP = DeletePhotonPackage(PP);
  for (PP = FinishedPhotonPackages->NextPackage; PP; PP = PP->NextPackage)
    PP = DeletePhotonPackage(PP);

  if (DeleteHeadPointer) {
    delete PhotonPackages;
    delete FinishedPhotonPackages;
    delete PausedPhotonPackages;
    PhotonPackages = NULL;
    FinishedPhotonPackages = NULL;
    PausedPhotonPackages = NULL;
  }
  else {
    PhotonPackages->NextPackage = NULL;
    PhotonPackages->PreviousPackage = NULL;
    FinishedPhotonPackages->NextPackage = NULL;
    FinishedPhotonPackages->PreviousPackage = NULL;
  }

  this->NumberOfPhotonPackages = 0;

  return SUCCESS;
  
}
