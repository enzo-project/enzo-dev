#define DEBUG 0
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

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "RadiativeTransferParameters.h"


#define MAX_CHANGE 0.1

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

float grid::ComputePhotonTimestepHII(float DensityUnits, float LengthUnits,
				     float VelocityUnits, float aye,
				     float Ifront_kph)
{

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return huge_number;

//  if (this->HasRadiation == FALSE)
//    return huge_number;

  int i, j, k, dim, index, indexn, size;
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
  float dxinv, dyinv, dzinv;
  float HIIdot, a3inv, a6inv, kr1, kr2, cs;
  float logtem, logtem0, logtem9, dlogtem, t1, t2, tdef, cs_factor;
  float *alldt = new float[size];
  for (i = 0; i < size; i++) alldt[i] = -huge_number;

  a3inv = 1.0/(a*a*a);
  a6inv = a3inv * a3inv;
  dxinv = 1.0 / CellWidth[0][0];
  dyinv = 1.0 / CellWidth[1][0];
  dzinv = 1.0 / CellWidth[2][0];
  cs_factor = 9.082e3 / VelocityUnits / Mu;

#ifdef USE_GRACKLE
  logtem0 = log(grackle_data->TemperatureStart);
  logtem9 = log(grackle_data->TemperatureEnd);
  nbins = grackle_data->NumberOfTemperatureBins;
  dlogtem = (logtem9 - logtem0) / float(nbins-1);
#else
  logtem0 = log(RateData.TemperatureStart);
  logtem9 = log(RateData.TemperatureEnd);
  nbins = RateData.NumberOfTemperatureBins;
  dlogtem = (logtem9 - logtem0) / float(nbins-1);
#endif

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
	if (BaryonField[kphHINum][index] <= Ifront_kph &&
	    BaryonField[kphHINum][index] > 0) {

	  logtem = min( max( log(temperature[index]), logtem0 ), logtem9 );
	  tidx = min( nbins-1, max(1, int((logtem - logtem0) / dlogtem)+1) );
	  t1 = logtem0 + (tidx - 1) * dlogtem;
	  t2 = logtem0 + (tidx    ) * dlogtem;
	  tdef = t2 - t1;

#ifdef USE_GRACKLE
	  kr1 = grackle_rates.k1[tidx] + (logtem - t1) * 
	    (grackle_rates.k1[tidx+1] - grackle_rates.k1[tidx]) / tdef;
	  kr2 = grackle_rates.k2[tidx] + (logtem - t1) * 
	    (grackle_rates.k2[tidx+1] - grackle_rates.k2[tidx]) / tdef;
#else
	  kr1 = RateData.k1[tidx] + (logtem - t1) * 
	    (RateData.k1[tidx+1] - RateData.k1[tidx]) / tdef;
	  kr2 = RateData.k2[tidx] + (logtem - t1) * 
	    (RateData.k2[tidx+1] - RateData.k2[tidx]) / tdef;
#endif

	  HIIdot = a6inv *
	    (kr1 * BaryonField[HINum][index] * BaryonField[DeNum][index] -
	     kr2 * BaryonField[HIINum][index] * BaryonField[DeNum][index]) +
	    a3inv * BaryonField[HINum][index] * BaryonField[kphHINum][index];

	  alldt[index] = MAX_CHANGE * BaryonField[HIINum][index] / HIIdot;

	} // ENDIF radiation > max(kph) in I-front (tau>~0.1)

	/* If not in the I-front, use the normal Godunov formula */
	else {
	  cs = cs_factor * sqrt(Gamma * temperature[index]);
	  alldt[index] = CourantSafetyNumber * aye / 
	    ((cs + fabs(BaryonField[Vel1Num][index])) * dxinv +
	     (cs + fabs(BaryonField[Vel2Num][index])) * dyinv + 
	     (cs + fabs(BaryonField[Vel3Num][index])) * dzinv);
	}

      } // ENDFOR i
    } // ENDFOR j

  /* Use a weighted averaged of the timestep in a 3^3 cube */

  int i0, j0, k0, i1, j1, k1, ii, jj, kk;
  int imin, ikernel, nx, ny, nz, nn;
  float weight, kernel2;
  
  const float kernel[] = 
    {0.25, 0.33, 0.25, 0.33, 0.50, 0.33, 0.25, 0.33, 0.25,
     0.33, 0.50, 0.33, 0.50, 1.00, 0.50, 0.33, 0.50, 0.33,
     0.25, 0.33, 0.25, 0.33, 0.50, 0.33, 0.25, 0.33, 0.25};

  imin = INT_UNDEFINED;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    k0 = max(k-1, GridStartIndex[2]);
    k1 = min(k+1, GridEndIndex[2]);
    nz = k1-k0+1;
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      j0 = max(j-1, GridStartIndex[1]);
      j1 = min(j+1, GridEndIndex[1]);
      ny = j1-j0+1;
      index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {

	if (BaryonField[kphHINum][index] <= Ifront_kph &&
	    BaryonField[kphHINum][index] > 0) {

	  i0 = max(i-1, GridStartIndex[0]);
	  i1 = min(i+1, GridEndIndex[1]);
	  nx = i1-i0+1;
	  weight = 0.0;
	  this_dt = 0.0;
	  for (kk = k0; kk <= k1; kk++)
	    for (jj = j0; jj <= j1; jj++) {
	      nn = i0-(i-1) + 3*( jj-(j-1) + 3*(kk-(k-1)) );
	      for (ii = i0; ii <= i1; ii++, nn++) {
		indexn = GRIDINDEX_NOGHOST(ii,jj,kk);
		if (alldt[indexn] > 0 && BaryonField[kphHINum][indexn] > 0) {
		  kernel2 = kernel[nn] * kernel[nn];
		  this_dt += kernel2 * alldt[indexn];
		  weight += kernel2;
		  //this_dt += kernel[nn] * alldt[indexn];
		  //weight += kernel[nn];
		}
	      }
	    }

	  this_dt = (weight>0) ? this_dt/weight : huge_number;
	  if (DEBUG)
	    if (this_dt < dt) imin = GRIDINDEX_NOGHOST(i,j,k);
	  dt = min(dt, this_dt);

	} // ENDIF kph (I-front)
      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  if (DEBUG)
    if (imin != INT_UNDEFINED)
      printf("kph=%g/%g, dt=%g, avg(dt)=%g, T=%g, HI=%g, HII=%g\n",
	     BaryonField[kphHINum][imin], Ifront_kph,
	     alldt[imin], dt, temperature[imin], 
	     BaryonField[HINum][imin] / BaryonField[0][imin],
	     BaryonField[HIINum][imin] / BaryonField[0][imin]);
  
  delete [] temperature;
  delete [] alldt;

  return dt;
}
