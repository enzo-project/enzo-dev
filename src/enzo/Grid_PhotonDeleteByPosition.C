#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (CHECKS IF PHOTON NEEDS TO BE DELETED USING POSITION CHECK)
/
/  written by: Andrew Emerick
/  date:       Feb, 2017
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: TRUE  = photon should stay alive
/           FALSE = photon should be deleted
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"


int grid::PhotonDeleteByPosition(int &cindex, FLOAT *r,
                                 PhotonPackageEntry* &PP, grid* &MoveToGrid,
                                 int &DeleteMe){

  float radius;
  int dim;

  // be careful and make sure we need to do this check
  if( !RadiativeTransferDeletePhotonByPosition) return TRUE;

  /* Check if we've left the photon sphere */
  for (dim = 0, radius = 0.0; dim <MAX_DIMENSION; dim++)
    radius += (r[dim] - 0.5 * (DomainRightEdge[dim] - DomainLeftEdge[dim]))*
              (r[dim] - 0.5 * (DomainRightEdge[dim] - DomainLeftEdge[dim]));
  radius = sqrt(radius);

  if (radius > RadiativeTransferDeletePhotonRadius){
    PP->Photons = -1;
    MoveToGrid = NULL;
    DeleteMe   = TRUE;
    return FALSE;
  }

  return TRUE;
}
