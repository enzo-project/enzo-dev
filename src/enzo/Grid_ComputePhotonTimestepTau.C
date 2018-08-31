#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (COMPUTE TIME STEP FOR RADIATIVE TRANSFER)
/
/  written by: John Wise
/  date:       September, 2009
/  modified1:  TA+JHW (July, 2010) -- compute maximum change in 
/              intensity, which is proportional to the I-front velocity.
/
/  PURPOSE: Calculates the photon timestep by restricting the change in
/           intensity to 50% in cells, so we capture the I-front evolution
/           correctly.
/
/  RETURNS:
/    dt   - photon timestep
/
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "RadiativeTransferParameters.h"

#define MAX_CHANGE 0.5

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

float grid::ComputePhotonTimestepTau(float DensityUnits, float LengthUnits,
				     float VelocityUnits, float aye)
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

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
			     Vel3Num, TENum);

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
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum);


  /* Compute temperature field */

  temperature = new float[size];
  this->ComputeTemperatureField(temperature);  

  int tidx, nbins;
  float HIIdot, a3inv, a6inv, kr1, kr2, tau;
  float logtem, logtem0, logtem9, dlogtem, t1, t2, tdef;

  a3inv = 1.0/(a*a*a);
  a6inv = a3inv * a3inv;

#ifdef USE_GRACKLE
  if (grackle_data->use_grackle == TRUE) {
    logtem0 = log(grackle_data->TemperatureStart);
    logtem9 = log(grackle_data->TemperatureEnd);
    nbins = grackle_data->NumberOfTemperatureBins;
    dlogtem = (logtem9 - logtem0) / float(nbins-1);
  } else
#endif
  {
    logtem0 = log(RateData.TemperatureStart);
    logtem9 = log(RateData.TemperatureEnd);
    nbins = RateData.NumberOfTemperatureBins;
    dlogtem = (logtem9 - logtem0) / float(nbins-1);
  }

  // Absorb all unit conversions into this factor, so we only have to
  // multiple by density.
  sigma_dx = float(sigmaHI * double(LengthUnits) * double(DensityUnits)/mh * 
		   double(CellWidth[0][0]));

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	/* Remember to convert densities from comoving to proper. */
	
	if (BaryonField[kphHINum][index] > 0) {

	  // optical depth in a cell
	  tau = sigma_dx * BaryonField[HINum][index];

	  logtem = min( max( log(temperature[index]), logtem0 ), logtem9 );
	  tidx = min( nbins-1, max(1, int((logtem - logtem0) / dlogtem)+1) );
	  t1 = logtem0 + (tidx - 1) * dlogtem;
	  t2 = logtem0 + (tidx    ) * dlogtem;
	  tdef = t2 - t1;

#ifdef USE_GRACKLE
	  if (grackle_data->use_grackle == TRUE) {
	    kr1 = grackle_rates.k1[tidx] + (logtem - t1) * 
	      (grackle_rates.k1[tidx+1] - grackle_rates.k1[tidx]) / tdef;
	    kr2 = grackle_rates.k2[tidx] + (logtem - t1) * 
	      (grackle_rates.k2[tidx+1] - grackle_rates.k2[tidx]) / tdef;
	  } else
#endif
	  {
	    kr1 = RateData.k1[tidx] + (logtem - t1) * 
	      (RateData.k1[tidx+1] - RateData.k1[tidx]) / tdef;
	    kr2 = RateData.k2[tidx] + (logtem - t1) * 
	      (RateData.k2[tidx+1] - RateData.k2[tidx]) / tdef;
	  }

	  HIIdot = a6inv *
	    (kr1 * BaryonField[HINum][index] * BaryonField[DeNum][index] -
	     kr2 * BaryonField[HIINum][index] * BaryonField[DeNum][index]) +
	    a3inv * BaryonField[HINum][index] * BaryonField[kphHINum][index];

	  this_dt = MAX_CHANGE * BaryonField[HINum][index] / 
	    (HIIdot * min(tau, 1.0));
	  if (this_dt > 0) dt = min(dt, this_dt);
#ifdef UNUSED
	  if (this_dt > 0) {
	    if (this_dt < dt) {
	      dt = this_dt;
	      printf("dtPhoton: ijk = %d %d %d, dt=%g, tau=%g, HIIdot=%g, HI=%g, kph=%g\n",
		     i,j,k,this_dt,tau, HIIdot, 
		     BaryonField[HINum][index],
		     BaryonField[kphHINum][index]);
	    }
	  }
#endif /* UNUSED */

	} // ENDIF radiation

      } // ENDFOR i
    } // ENDFOR j
  
  delete [] temperature;

  return dt;
}
