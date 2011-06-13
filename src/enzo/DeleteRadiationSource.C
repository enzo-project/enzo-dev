/***********************************************************************
/
/  DELETES RADIATION SOURCES
/
/  written by: Tom Abel
/  date:       April, 2004
/  modified1:
/
/  PURPOSE:  Deletes radiation source from linked list and frees memory
/            Identical to DeletePhotonPackage
/  INPUTS:
/    RadiationSourceEntry - to delete
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "RadiationSource.h"

RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS)
{
  if (RS->PreviousSource == NULL) {
    fprintf(stderr, "DeleteRadiationSource: Warning Previous Radiation Source == NULL");
  }

  RadiationSourceEntry* dummy;
  //end of list ? 
  if (RS->NextSource == NULL) {
    RS->PreviousSource->NextSource = NULL;
    dummy = NULL;
    delete [] RS->Position;
    delete [] RS->SED;
    delete [] RS->Energy;
    delete [] RS->Orientation;
    delete RS;
  } 
  // beginning
  else if (RS->PreviousSource == NULL) { 
    fprintf(stderr,"delete head of list???\n");
    return RS;
  } else {
    (RS->PreviousSource)->NextSource = RS->NextSource;
    (RS->NextSource)->PreviousSource = RS->PreviousSource;
    dummy = RS->NextSource;
    delete [] RS->Position;
    delete [] RS->SED;
    delete [] RS->Energy;
    delete [] RS->Orientation;
    delete RS;
  }
  return dummy;
}
