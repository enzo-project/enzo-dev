/***********************************************************************
/
/  GRID CLASS (SET BOUNDARY FLUXES TO ZERO)
/
/  written by: Peng Wang
/  date:       September, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"

void grid::ZeroFluxes(fluxes *SubgridFluxes[], int NumberOfSubgrids)
{
  if (ProcessorNumber != MyProcessorNumber)
    return;

  int fluxsize;
  for (int subgrid = 0; subgrid < NumberOfSubgrids; subgrid++) {
    for (int flux = 0; flux < GridRank; flux++) {
      fluxsize = 1;
      for (int j = 0; j < GridRank; j++) {
        fluxsize *= SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[flux][j] -
          SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[flux][j] + 1;
      }
      for (int field = 0; field < NumberOfBaryonFields; field++) {
        for (int n = 0; n < fluxsize; n++) {
          SubgridFluxes[subgrid]->LeftFluxes[field][flux][n] = 0.0;
          SubgridFluxes[subgrid]->RightFluxes[field][flux][n] = 0.0;
        }
      }
    }
  }

  return;
}
