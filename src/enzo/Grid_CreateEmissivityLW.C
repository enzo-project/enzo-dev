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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::CreateEmissivityLW(Star *AllStars, FLOAT TimeFLD, float dtFLD)
{

  if (RadiativeTransfer == FALSE || RadiativeTransferFLD == FALSE) 
    return SUCCESS;

  Star *cstar;
  int i, j, k, index, dim, size, EtaNum, nbins;
  float energies[MAX_ENERGY_BINS], E_LW, TimeFraction;
  double Luminosity[MAX_ENERGY_BINS], L_LW, CellVolume;
  const float ev_erg = 1.602176e-12;

  /* Allocate and compute emissivity field */

  EtaNum = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Emissivity0;

  if (MyProcessorNumber == ProcessorNumber) {

    size = 1;
    CellVolume = 1.0;
    for (dim = 0; dim < GridRank; dim++) {
      size *= GridDimension[dim];
      CellVolume *= CellWidth[dim][0];
    }

    /* Get units */

    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
      VelocityUnits;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
	     &VelocityUnits, Time);

    CellVolume *= POW(LengthUnits, 3.0);

    BaryonField[EtaNum] = new float[size];
    for (i = 0; i < size; i++)
      BaryonField[EtaNum][i] = 0.0;

    /* Loop over all stars (because these may be contained in child
       subgrids) and fill in emissivity field */

    if (ProblemType != 50) { // Photon test

      for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

	if (this->PointInGrid(cstar->pos) == FALSE || 
	    cstar->IsActive() == FALSE)
	  continue;

	/* Get if the star will die in the next timestep and calculate
	   factor to reduce the emissivity */

	if (cstar->BirthTime + cstar->LifeTime < TimeFLD + dtFLD)
	  TimeFraction = (cstar->BirthTime + cstar->LifeTime - TimeFLD) / 
	    dtFLD;
	else
	  TimeFraction = 1.0;

	cstar->ComputePhotonRates(nbins, energies, Luminosity);
	E_LW = energies[3];
	L_LW = Luminosity[3];

	i = int((cstar->pos[0] - GridLeftEdge[0]) / CellWidth[0][0]);
	j = int((cstar->pos[1] - GridLeftEdge[1]) / CellWidth[1][0]);
	k = int((cstar->pos[2] - GridLeftEdge[2]) / CellWidth[2][0]);
	index = GRIDINDEX(i,j,k);

	BaryonField[EtaNum][index] += L_LW * E_LW * ev_erg / CellVolume * 
	  TimeFraction;

      } // ENDFOR stars

    } // ENDIF ProblemType != 50

    // Photon test problem
    else {

      RadiationSourceEntry *RS;

      // Convert from #/s to RT units
      double LConv = (double) TimeUnits / pow(LengthUnits,3);

      for (RS = GlobalRadiationSources->NextSource; RS; RS = RS->NextSource) {

	if (this->PointInGrid(RS->Position) == FALSE)
	  continue;

	/* Get if the star will die in the next timestep and calculate
	   factor to reduce the emissivity */

	if (RS->CreationTime + RS->LifeTime < TimeFLD + dtFLD)
	  TimeFraction = (RS->CreationTime + RS->LifeTime - TimeFLD) / dtFLD;
	else
	  TimeFraction = 1.0;

	E_LW = RS->Energy[3];
	L_LW = RS->Luminosity * RS->SED[3] / LConv;

	i = int((RS->Position[0] - GridLeftEdge[0]) / CellWidth[0][0]);
	j = int((RS->Position[1] - GridLeftEdge[1]) / CellWidth[1][0]);
	k = int((RS->Position[2] - GridLeftEdge[2]) / CellWidth[2][0]);
	index = GRIDINDEX(i,j,k);

	BaryonField[EtaNum][index] += L_LW * E_LW * ev_erg / CellVolume * 
	  TimeFraction;

      } // ENDFOR stars

    }

  } // ENDIF MyProcessor == ProcessorNumber
  
	 
  return SUCCESS;
}
