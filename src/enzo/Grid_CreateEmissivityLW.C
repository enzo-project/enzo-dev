#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (COMPUTE EMISSIVITY FROM STARS)
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
/  PURPOSE: 
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

int FindField(int f, int farray[], int n);

int grid::CreateEmissivityLW(Star *AllStars)
{

  if (RadiativeTransfer == FALSE || RadiativeTransferFLD == FALSE) 
    return SUCCESS;

  Star *cstar;
  int i, j, k, index, dim, size, EtaNum;
  float energies[4], E_LW;
  double Luminosity[4], L_LW;
  const float ev_erg = 1.602176e-12;

  /* Allocate and compute emissivity field */

  EtaNum = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Emissivity0;

  if (MyProcessorNumber == ProcessorNumber) {

    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

    BaryonField[EtaNum] = new float[size];
    for (i = 0; i < size; i++)
      BaryonField[EtaNum][i] = 0.0;

    /* Loop over all stars (because these may be contained in child
       subgrids) and fill in emissivity field */

    for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

      if (this->PointInGrid(cstar->pos) == FALSE)
	continue;

      cstar->ComputePhotonRates(energies, Luminosity);
      E_LW = energies[3];
      L_LW = Luminosity[3];

      i = int((cstar->pos[0] - CellLeftEdge[0][0]) / CellWidth[0][0]);
      j = int((cstar->pos[1] - CellLeftEdge[1][0]) / CellWidth[1][0]);
      k = int((cstar->pos[2] - CellLeftEdge[2][0]) / CellWidth[2][0]);
      index = GRIDINDEX(i,j,k);

      BaryonField[EtaNum][index] += L_LW * E_LW * ev_erg;

    } // ENDFOR stars

  } // ENDIF MyProcessor == ProcessorNumber
  
	 
  return SUCCESS;
}
