#define DEBUG 0
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "PhotonPackage.h"

int SplitPhotonPackage(PhotonPackageEntry *PP)
{

  int childrays, dim;
  long long nipix;
  PhotonPackageEntry *NewPack;

  if (PP->PreviousPackage == NULL)
    fprintf(stderr, "SplitPhotonPackage: Warning Previous Photon Package == NULL");
  
  if (DEBUG) 
    fprintf(stdout, "split package ipix:%ld level:%ld\n", (long)PP->ipix, (long)PP->level);

  nipix = (PP->ipix)*4;
  for (childrays=0; childrays < 4; childrays++) {
    NewPack = new PhotonPackageEntry;

    // insert pointer in list just after PP
    NewPack->PreviousPackage = PP;
    NewPack->NextPackage     = PP->NextPackage;
    if (PP->NextPackage != NULL) PP->NextPackage->PreviousPackage = NewPack;
    PP->NextPackage      = NewPack;
    if (DEBUG) 
      fprintf(stdout, "split package %ld %p %p %p %p\n",
	      (long)PP->ipix, (void*)PP->PreviousPackage,
	      (void*)NewPack->PreviousPackage, (void*)NewPack, (void*)NewPack->NextPackage);

    NewPack->Photons         = 0.25*PP->Photons;
    NewPack->Type            = PP->Type;
    NewPack->EmissionTimeInterval = PP->EmissionTimeInterval;
    NewPack->EmissionTime    = PP->EmissionTime;
    NewPack->CurrentTime     = PP->CurrentTime;
    NewPack->ColumnDensity   = PP->ColumnDensity;
    NewPack->Radius          = PP->Radius;
    NewPack->ipix            = nipix++;
    NewPack->level           = PP->level+1;
    NewPack->Energy          = PP->Energy;
    NewPack->CrossSection    = PP->CrossSection;
    for (dim = 0; dim < 3; dim++)
      NewPack->SourcePosition[dim] = PP->SourcePosition[dim];
    NewPack->SourcePositionDiff  = PP->SourcePositionDiff;
    NewPack->CurrentSource   = PP->CurrentSource;

    if ((NewPack->PreviousPackage->NextPackage != NewPack)) {
      ENZO_VFAIL("SplitPhotonPackage: Problem splitting %p %p %p %p %p \n",
	      (void*)PP, (void*)NewPack, (void*)NewPack->PreviousPackage,
	      (void*)NewPack->PreviousPackage->NextPackage, (void*)NewPack->NextPackage)

    }
  } // for childrays=0,3
    
  
  return SUCCESS;
}
