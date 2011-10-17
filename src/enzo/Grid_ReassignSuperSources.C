#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

int FindSuperSource(PhotonPackageEntry **PP, int &LeafID, int SearchNewTree);
int FindSuperSourceByPosition(PhotonPackageEntry **PP);

int grid::ReassignSuperSources(void)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (PhotonPackages->NextPackage == NULL)
    return SUCCESS;

  PhotonPackageEntry *PP;
  int dim, LeafID;
  float radius2, dx;
  FLOAT OldPosition[MAX_DIMENSION];

  for (PP = PhotonPackages->NextPackage; PP; PP = PP->NextPackage) {

    if (PP->CurrentSource == NULL)
      continue;

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      OldPosition[dim] = PP->CurrentSource->Position[dim];

    // Reassign super source by leaf ID
    LeafID = PP->CurrentSource->LeafID;
    if (FindSuperSource(&PP, LeafID, FALSE) == FAIL) {
      ENZO_FAIL("Error in FindSuperSource.\n");
    }

    radius2 = 0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      dx = PP->CurrentSource->Position[dim] - OldPosition[dim];
      radius2 += dx*dx;
    }

    /* In the case where the leaf ID has changed (check by change in
       position), find super source by position.  If LeafID is
       undefined, then FindSuperSource couldn't locate the leaf by its
       ID, so we search by position. */

    if (radius2 > PP->CurrentSource->ClusteringRadius *
	PP->CurrentSource->ClusteringRadius ||
	LeafID == INT_UNDEFINED)
      if (FindSuperSourceByPosition(&PP) == FAIL) {
	ENZO_FAIL("Error in FindSuperSourceByPosition.\n");

      }

  } // ENDFOR photons

  return SUCCESS;

}
