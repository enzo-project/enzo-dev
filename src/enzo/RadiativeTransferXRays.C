/***********************************************************************
/
/  Calculate the effect of X-Ray radiation on the species
/  HI, HeI and HeII
/
/  written by: John Wise
/  date:     
/  modified1: John Regan
/             Moved to its own function and file from
/             WalkPhotonPackage
/
/  PURPOSE: Calculate the ionisations and xray photons absorbed by HI, HeI and HeII
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "CosmologyParameters.h"

int grid::RadiativeTransferXRays(PhotonPackageEntry **PP, FLOAT *dPXray, int cellindex,
				 int species, FLOAT ddr, FLOAT tau,  float geo_correction,
				 FLOAT photonrate, FLOAT *excessrate, float *ion2_factor,
				 float heat_factor, const int *kphNum, int gammaNum)
{  
  float dP1 = 0.0;
	
  // at most use all photons for photo-ionizations
  if (tau > 2.e1) 
    dPXray[species] = (1.0+BFLOAT_EPSILON) * (*PP)->Photons;
  else if (tau > 1.e-4) 
    dPXray[species] = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
  else
    dPXray[species] = min((*PP)->Photons*tau, (*PP)->Photons);
  
  dP1 = dPXray[species] * geo_correction;

  // contributions to the photoionization rate is over whole timestep
  // units are (1/LengthUnits^3)*(1/CodeTime)
  // This needs to be normalised - see Grid_FinalizeRadiationFields.C
  BaryonField[kphNum[species]][cellindex] += dP1 * photonrate * ion2_factor[species];
	
  // the heating rate is just the number of photo ionizations times
  // the excess energy; units are eV/CodeTime*((1/LengthUnits^3)); 
  // check Grid_FinalizeRadiationFields.C
  BaryonField[gammaNum][cellindex] += dP1 * excessrate[species] * heat_factor;
  
  return SUCCESS;
}

#define COMPTON 3

int grid::RadiativeTransferComptonHeating(PhotonPackageEntry **PP, FLOAT *dPXray, int cellindex,
					  float LengthUnits, float photonrate, 
					  int TemperatureField, FLOAT ddr, double dN, 
					  float geo_correction, int gammaNum)
{
  FLOAT xE = 0.0, ratioE = 0.0, dP1 = 0.0, xray_sigma = 0.0;
  FLOAT excess_heating = 0.0;
  const double k_b = 8.62e-5; // eV/K
  float tau = 0.0;
  // assume photon energy is much less than the electron rest mass energy 
  // nonrelativistic Klein-Nishina cross-section in Rybicki & Lightman (1979)
  xE = (*PP)->Energy/5.11e5;  // h*nu/(m_e*c^2)
  xray_sigma = 6.65e-25 * (1 - 2.*xE + 26./5.*xE*xE) * LengthUnits; //Equation 7.6a
  
  // also, nonrelativistic energy transfer in Ciotti & Ostriker (2001)
  excess_heating = photonrate * 4 * k_b * BaryonField[TemperatureField][cellindex] * xE;
  ratioE = 4 * k_b * BaryonField[TemperatureField][cellindex] * xE / (*PP)->Energy; 
  
  tau = dN*xray_sigma;
  
  // at most use all photons for Compton scattering
  if (tau > 2.e1) 
    dPXray[COMPTON] = (1.0+BFLOAT_EPSILON) * (*PP)->Photons;
  else if (tau > 1.e-4) 
    dPXray[COMPTON] = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
  else
    dPXray[COMPTON] = min((*PP)->Photons*tau, (*PP)->Photons);
  dP1 = dPXray[COMPTON] * geo_correction;

  // the heating rate by energy transfer during Compton scattering
  // [dP1] = 1/LengthUnits^3
  // [excess_heating] = eV/CodeTime
  // [BaryonField[gammaNum]] = eV/CodeTime/LengthUnits^3
  // This needs to be nomalised - see Grid_FinalizeRadiationFields.C
  BaryonField[gammaNum][cellindex] += dP1 * excess_heating; 
  
  // a photon loses only a fraction of photon energy in Compton scatering, 
  // and keeps propagating; to model this with monochromatic energy,
  // we instead subtract #photons (dPXray[COMPTON]_new) from PP
  // (photon energy absorbed) = dPXray[COMPTON]     * (4*k_B*T*xE) 
  //                          = dPXray[COMPTON]_new * (*PP)->Energy
  // See also Kim et al. (2011)
  dPXray[COMPTON] *= ratioE;
  
  
  return SUCCESS;
}
