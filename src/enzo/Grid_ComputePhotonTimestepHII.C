/***********************************************************************
/
/  GRID CLASS (COMPUTE TIME STEP FOR RADIATIVE TRANSFER)
/
/  written by: John Wise
/  date:       September, 2009
/  modified1:
/
/  PURPOSE: Calculates the photon timestep by restricting the change in
/           ionizing hydrogen to 50% in cells with radiation and optical
/           depths greater than 0.5, so we capture the I-front evolution
/           correctly.
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

#define MAX_CHANGE 0.1

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

float grid::ComputePhotonTimestepHII(float DensityUnits, float LengthUnits)
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return huge_number;

//  if (this->HasRadiation == FALSE)
//    return huge_number;

  int i, j, k, dim, index, size;
  float dt, this_dt, sigma_dx, *temperature;
  const double sigmaHI = 6.345e-18, mh = 1.673e-24;

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
  float HIIdot, a3inv, a6inv, k1, k2, tau;
  float logtem, logtem0, logtem9, dlogtem, t1, t2, tdef;

  a3inv = 1.0/(a*a*a);
  a6inv = a3inv * a3inv;

  logtem0 = log(RateData.TemperatureStart);
  logtem9 = log(RateData.TemperatureEnd);
  dlogtem = (logtem9 - logtem0) / float(RateData.NumberOfTemperatureBins-1);

  // Absorb all unit conversions into this factor, so we only have to
  // multiple by density.
//  sigma_dx = float(sigmaHI * double(LengthUnits) * double(DensityUnits)/mh * 
//		   double(CellWidth[0][0]));

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	/* Remember to convert densities from comoving to proper.  We
	   also only consider cells with cumulative optical depths >
	   0.1, i.e. near the I-front.  This is quantified by a
	   maximum kph in the grid and is calculated in
	   WalkPhotonPackage. */
	
	//tau = sigma_dx * BaryonField[HINum][index];
	if (BaryonField[kphHINum][index] < this->MaximumkphIfront &&
	    BaryonField[kphHINum][index] > 0) {

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

	  if (HIIdot > 0) {
	    this_dt = MAX_CHANGE * BaryonField[HIINum][index] / HIIdot;
//	    printf("kph=%g/%g, dt=%g, T=%g, HI=%g, HII=%g, k1=%g, k2=%g\n",
//		   BaryonField[kphHINum][index], MaximumkphIfront,
//		   this_dt, temperature[index], 
//		   BaryonField[HINum][index] / BaryonField[0][index],
//		   BaryonField[HIINum][index] / BaryonField[0][index],
//		   k1, k2);
	    dt = min(dt, this_dt);
	  }

	} // ENDIF radiation > max(kph) in I-front (tau>~0.1)

      } // ENDFOR i
    } // ENDFOR j

  return dt;
}
