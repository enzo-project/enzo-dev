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
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);

int grid::DeletePhotonPackages() {

  PhotonPackageEntry *PP = PhotonPackages;

  while (PP != NULL) PP = DeletePhotonPackage(PP);

  PhotonPackages = new PhotonPackageEntry;
  PhotonPackages->NextPackage = NULL;
  PhotonPackages->PreviousPackage = NULL;

  return SUCCESS;
  
}
