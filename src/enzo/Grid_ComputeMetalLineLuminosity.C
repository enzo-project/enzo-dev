/***********************************************************************
/
/  GRID CLASS (COMPUTE THE METAL LINE EMISSIVITY FIELD)
/
/  written by: John Wise
/  date:       May, 2008
/  modified1:
/
/  PURPOSE:  Compute metal line emissivity for the n-th line in the 
/            lookup table.
/
/  RETURNS:
/
************************************************************************/

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
#include "fortran.def"
#include "Grid.h"
#include "CosmologyParameters.h"

#define NUMBER_OF_COOLANTS 11

/* function prototypes */

int FindField(int f, int farray[], int n);

int grid::ComputeMetalLineLuminosity(float *total_luminosity, float *all_emis, 
				     float *temperature)
{

  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber || NumberOfBaryonFields == 0)
    return SUCCESS;

  /* Compute the size of the fields. */

  int dim, i, j, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find electron and total density to calculate electron fraction */

  int DensNum, DeNum;
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0)
    ENZO_FAIL("Cannot find density.");

  if ((DeNum = FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) < 0)
    ENZO_FAIL("Cannot find electron density.");

  float dlogtem, dlogxe, x_e;
  int xebin, tbin, bin;

  dlogtem = log(CoolData.MR_TemperatureEnd / CoolData.MR_TemperatureStart) 
    / float(CoolData.MR_NumberOfTemperatureBins-1);
  dlogxe = log(CoolData.MR_ElectronFracEnd / CoolData.MR_ElectronFracStart)
    / float(CoolData.MR_NumberOfElectronFracBins-1);

  for (i = 0; i < size; i++) {

    x_e = BaryonField[DeNum][i] / BaryonField[DensNum][i];

    /* Find table cell that corresponds to this x_e and temperature.
       No interpolation for now. */

    tbin = log(temperature[i] / CoolData.MR_TemperatureStart) / dlogtem;
    xebin = log(x_e / CoolData.MR_ElectronFracStart) / dlogxe;
    bin = (xebin * CoolData.MR_NumberOfTemperatureBins + tbin) * 
      NUMBER_OF_COOLANTS;

    for (j = 0; j < NUMBER_OF_COOLANTS; j++)
      all_emis[j*size+i] = total_luminosity[i] * CoolData.metal_ratios[bin+j];

  } // ENDFOR i

  return SUCCESS;
}
