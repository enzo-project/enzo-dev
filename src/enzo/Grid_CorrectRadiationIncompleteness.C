#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (CORRECT THE RADIATION FIELD FOR INCOMPLETE SAMPLING)
/
/  written by: John Wise
/  date:       June, 2009
/  modified1:
/
/  PURPOSE: 
/
************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"

int grid::CorrectRadiationIncompleteness(void)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

#ifdef TRANSFER

  int i, dim, size = 1;

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find radiative transfer fields. */

  int kphHINum, gammaHINum, kphHeINum, gammaHeINum, kphHeIINum, gammaHeIINum,
    kdissH2INum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				      gammaHeINum, kphHeIINum, gammaHeIINum, 
				      kdissH2INum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifyRadiativeTransferFields.\n");
    return FAIL;
  }

  int nsrc;

  for (i = 0; i < size; i++)
    if (BaryonField[kphHeIINum][i] > 0) {
      nsrc = max(int(BaryonField[kphHeIINum][i]-0.1), 1);
      BaryonField[kphHINum][i] /= BaryonField[kphHeIINum][i] / nsrc;
      BaryonField[gammaHINum][i] /= BaryonField[kphHeIINum][i] / nsrc;
    }

#endif /* TRANSFER */  
  
  return SUCCESS;
}
