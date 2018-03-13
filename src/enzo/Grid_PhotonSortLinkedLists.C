/***********************************************************************
/
/  GRID CLASS (SORT LINKED LISTS OF PHOTONS)
/
/  written by: Stephen Skory
/  date:       June, 2011
/  modified1:
/
/  PURPOSE:  Sorts the linked lists of arrays.
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

Eint32 compare_ss(const void *a, const void *b);
PhotonPackageEntry *LinkedListToArray(PhotonPackageEntry *Node, int n);
PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);

int grid::PhotonSortLinkedLists(void)
{
  // Photons of all three types are sorted following the method used
  // in Grid_MergePausedPhotonPackages.

  if (MyProcessorNumber != ProcessorNumber) return SUCCESS;
  
  int nphotons, dim, count;
  PhotonPackageEntry *PP, *TempPP, *NewPack;
  
  // PhotonPackages
  nphotons = 0;
  PP = PhotonPackages->NextPackage;
  while (PP != NULL) {
    PP = PP->NextPackage;
    nphotons++;
  }
  
  PP = PhotonPackages->NextPackage;
  TempPP = LinkedListToArray(PP, nphotons);
  // Sort the list.
  qsort(TempPP, nphotons, sizeof(PhotonPackageEntry), compare_ss);
  
  // Now we need to clear out PhotonPackages before we fill it back up again.
  for (PP = PhotonPackages->NextPackage; PP; PP = PP->NextPackage)
    PP = DeletePhotonPackage(PP);
  PhotonPackages->NextPackage = NULL;
  PhotonPackages->PreviousPackage = NULL;

  // Fill it back in from the sorted list.
  PP = PhotonPackages;
  for (count = 0; count < nphotons; count++) {
    NewPack = new PhotonPackageEntry;
    NewPack->CurrentSource = TempPP[count].CurrentSource;
    NewPack->Photons = TempPP[count].Photons;
    NewPack->Type = TempPP[count].Type;
    NewPack->Energy = TempPP[count].Energy;
    NewPack->CrossSection = TempPP[count].CrossSection; //
    NewPack->EmissionTimeInterval = TempPP[count].EmissionTimeInterval;
    NewPack->EmissionTime = TempPP[count].EmissionTime;
    NewPack->CurrentTime = TempPP[count].CurrentTime;
    NewPack->Radius = TempPP[count].Radius;
    NewPack->ColumnDensity = TempPP[count].ColumnDensity;
    NewPack->ipix = TempPP[count].ipix;
    NewPack->level = TempPP[count].level;
    NewPack->SourcePositionDiff = TempPP[count].SourcePositionDiff;
    for (dim = 0; dim < 3; dim++) {
      NewPack->SourcePosition[dim] = TempPP[count].SourcePosition[dim];
    }
    InsertPhotonAfter(PP, NewPack);
    PP = PP->NextPackage;
  }
  
  // clean up
  delete [] TempPP;
  

  // FinishedPhotonPackages
  nphotons = 0;
  PP = FinishedPhotonPackages->NextPackage;
  while (PP != NULL) {
    PP = PP->NextPackage;
    nphotons++;
  }
  
  PP = FinishedPhotonPackages->NextPackage;
  TempPP = LinkedListToArray(PP, nphotons);
  // Sort the list.
  qsort(TempPP, nphotons, sizeof(PhotonPackageEntry), compare_ss);
  
  // Now we need to clear out FinishedPhotonPackages
  // before we fill it back up again.
  for (PP = FinishedPhotonPackages->NextPackage; PP; PP = PP->NextPackage)
    PP = DeletePhotonPackage(PP);
  FinishedPhotonPackages->NextPackage = NULL;
  FinishedPhotonPackages->PreviousPackage = NULL;

  // Fill it back in from the sorted list.
  PP = FinishedPhotonPackages;
  for (count = 0; count < nphotons; count++) {
    NewPack = new PhotonPackageEntry;
    NewPack->CurrentSource = TempPP[count].CurrentSource;
    NewPack->Photons = TempPP[count].Photons;
    NewPack->Type = TempPP[count].Type;
    NewPack->Energy = TempPP[count].Energy;
    NewPack->CrossSection = TempPP[count].CrossSection;
    NewPack->EmissionTimeInterval = TempPP[count].EmissionTimeInterval;
    NewPack->EmissionTime = TempPP[count].EmissionTime;
    NewPack->CurrentTime = TempPP[count].CurrentTime;
    NewPack->Radius = TempPP[count].Radius;
    NewPack->ColumnDensity = TempPP[count].ColumnDensity;
    NewPack->ipix = TempPP[count].ipix;
    NewPack->level = TempPP[count].level;
    NewPack->SourcePositionDiff = TempPP[count].SourcePositionDiff;
    for (dim = 0; dim < 3; dim++) {
      NewPack->SourcePosition[dim] = TempPP[count].SourcePosition[dim];
    }
    InsertPhotonAfter(PP, NewPack);
    PP = PP->NextPackage;
  }
  
  // clean up
  delete [] TempPP;

  // PausedPhotonPackages
  nphotons = 0;
  PP = PausedPhotonPackages->NextPackage;
  while (PP != NULL) {
    PP = PP->NextPackage;
    nphotons++;
  }
  
  PP = PausedPhotonPackages->NextPackage;
  TempPP = LinkedListToArray(PP, nphotons);
  // Sort the list.
  qsort(TempPP, nphotons, sizeof(PhotonPackageEntry), compare_ss);
  
  // Now we need to clear out PausedPhotonPackages
  // before we fill it back up again.
  for (PP = PausedPhotonPackages->NextPackage; PP; PP = PP->NextPackage)
    PP = DeletePhotonPackage(PP);
  PausedPhotonPackages->NextPackage = NULL;
  PausedPhotonPackages->PreviousPackage = NULL;

  // Fill it back in from the sorted list.
  PP = PausedPhotonPackages;
  for (count = 0; count < nphotons; count++) {
    NewPack = new PhotonPackageEntry;
    NewPack->CurrentSource = TempPP[count].CurrentSource;
    NewPack->Photons = TempPP[count].Photons;
    NewPack->Type = TempPP[count].Type;
    NewPack->Energy = TempPP[count].Energy;
    NewPack->CrossSection = TempPP[count].CrossSection;
    NewPack->EmissionTimeInterval = TempPP[count].EmissionTimeInterval;
    NewPack->EmissionTime = TempPP[count].EmissionTime;
    NewPack->CurrentTime = TempPP[count].CurrentTime;
    NewPack->Radius = TempPP[count].Radius;
    NewPack->ColumnDensity = TempPP[count].ColumnDensity;
    NewPack->ipix = TempPP[count].ipix;
    NewPack->level = TempPP[count].level;
    NewPack->SourcePositionDiff = TempPP[count].SourcePositionDiff;
    for (dim = 0; dim < 3; dim++) {
      NewPack->SourcePosition[dim] = TempPP[count].SourcePosition[dim];
    }
    InsertPhotonAfter(PP, NewPack);
    PP = PP->NextPackage;
  }
  
  // clean up
  delete [] TempPP;

  return SUCCESS;
}
