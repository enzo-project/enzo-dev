/***********************************************************************
/
/  PHOTON PACKAGE ROUTINES
/
/  written by: John Wise
/  date:       February, 2010
/  modified1:  
/
/  PURPOSE: Constructs a Linked List of Photon Packages including data
/
************************************************************************/
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "PhotonPackage.h"

PhotonPackageEntry::PhotonPackageEntry(void)
{
  NextPackage = NULL;
  PreviousPackage = NULL;
  CurrentSource = NULL;
  Photons = 0.0;
  Type = 0;
  Energy = 0.0;
  CrossSection = 0.0;
  EmissionTimeInterval = 0.0;
  EmissionTime = 0.0;
  CurrentTime = 0.0;
  Radius = 0.0;
  ColumnDensity = 0.0;
  ipix = 0;
  level = 0;
  SourcePosition[0] = 0.0;
  SourcePosition[1] = 0.0;
  SourcePosition[2] = 0.0;
  SourcePositionDiff = 0.0;
}

/**********************************************************************/

#ifdef MEMORY_POOL
void* PhotonPackageEntry::operator new(size_t object_size)
{
  return PhotonMemoryPool->GetMemory(object_size);
}

void PhotonPackageEntry::operator delete(void* object)
{
  PhotonMemoryPool->FreeMemory(object);
}
#endif
