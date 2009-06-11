/***********************************************************************
/
/  DELETES PHOTON PACKAGE
/
/  written by: Tom Abel
/  date:       April, 2004
/  modified1:
/
/  PURPOSE:  Deletes photon package from linked list and frees memory
/
/  INPUTS:
/    PhotonPackageEntry - to delete
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "PhotonPackage.h"


PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP)
{
//  if (PP->PreviousPackage == NULL) {
//    fprintf(stderr, "DeletePhotonPackage: Warning Previous Photon Package == NULL");
//  }

  PhotonPackageEntry* dummy;
  //end of list ? 
  if (PP->NextPackage == NULL) {
    PP->PreviousPackage->NextPackage = NULL;
    dummy = PP->PreviousPackage;
    delete PP;
  } else  // beginning
    if (PP->PreviousPackage == NULL) { 
      //      fprintf(stderr,"delete head of list???\n");
      return PP;
    } else {
      {
	(PP->PreviousPackage)->NextPackage = PP->NextPackage;
	(PP->NextPackage)->PreviousPackage = PP->PreviousPackage;
	dummy = PP->PreviousPackage;
	delete PP;
      }
    }
  return dummy;
}
