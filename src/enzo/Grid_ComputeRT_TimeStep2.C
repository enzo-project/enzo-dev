/***********************************************************************
/
/  GRID CLASS (COMPUTE TIME STEP FOR RADIATIVE TRANSFER)
/
/  written by: John Wise
/  date:       September, 2009
/  modified1:
/
/  PURPOSE: Calculates the photon timestep by restricting the change in
/           ionizing hydrogen to 10% in cells with radiation.
/
/  RETURNS:
/    dt   - photon timestep
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
#include "Grid.h"
#include "StarParticleData.h"
#include "RadiativeTransferParameters.h"

#define MAX_CHANGE 0.5

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

float grid::ComputeRT_TimeStep2(void)
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return huge_number;

//  if (this->HasRadiation == FALSE)
//    return huge_number;

  int i, j, k, dim, index, size;
  float dt, *temperature;

  dt = huge_number;
  
  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Find Multi-species fields and expansion factor. */
 
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  FLOAT a = 1.0, dadt;
 
  IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
			HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
 
  /* Find photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaHINum, gammaHeINum, gammaHeIINum;
  IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				  gammaHeINum, kphHeIINum, gammaHeIINum, 
				  kdissH2INum);


  /* Compute temperature field */

  temperature = new float[size];
  this->ComputeTemperatureField(temperature);  

  int tidx;
  float HIIdot, a3inv, a6inv, k1, k2;
  float logtem, logtem0, logtem9, dlogtem, t1, t2, tdef;

  a3inv = 1.0/(a*a*a);
  a6inv = a3inv * a3inv;

  logtem0 = log(RateData.TemperatureStart);
  logtem9 = log(RateData.TemperatureEnd);
  dlogtem = (logtem9 - logtem0) / float(RateData.NumberOfTemperatureBins-1);

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	/* Remember to convert densities from comoving to proper */
	
	if (BaryonField[kphHINum][index] > 0) {

	  logtem = min( max( log(temperature[index]), logtem0 ), logtem9 );
	  tidx = min( RateData.NumberOfTemperatureBins-1,
		      max(1, int((logtem - logtem0) / dlogtem)+1) );
	  t1 = logtem0 + (tidx - 1) * dlogtem;
	  t2 = logtem0 + (tidx    ) * dlogtem;
	  tdef = t2 - t1;

	  k1 = RateData.k1[tidx] + (logtem - t1) * 
	    (RateData.k1[tidx+1] - RateData.k1[tidx]) / tdef;
	  k2 = RateData.k2[tidx] + (logtem - t1) * 
	    (RateData.k2[tidx+1] - RateData.k2[tidx]) / tdef;

	  HIIdot = a6inv *
	    (k1 * BaryonField[HINum][index] * BaryonField[DeNum][index] -
	     k2 * BaryonField[HIINum][index] * BaryonField[DeNum][index]) +
	    a3inv * BaryonField[HINum][index] * BaryonField[kphHINum][index];

	  if (HIIdot > 0)
	    dt = min(dt, MAX_CHANGE * BaryonField[HIINum][index] / HIIdot);

	}

      } // ENDFOR i
    } // ENDFOR j


  return dt;
}
