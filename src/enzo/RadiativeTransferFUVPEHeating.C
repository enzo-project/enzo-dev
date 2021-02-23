#define DEBUG 0
/***********************************************************************
/
/  Calculate the propogation of FUV radiation associated specifically
/  with photoelectic heating on dust grains.
/
/  written by: Andrew Emerick - based on the previous work of John Wise + John Regan
/  date:     Oct 2019
/  modified1:
/
/  PURPOSE:
/  This accounts for the portion
/  of the spectrum NOT attenuated by H2 (5.6-11.2 eV). 11.2 - 13.6 eV is handled
/  in RadiativeTransferLW and RadiativeTransferH2II. Currently this only makes sense
/  if H2 is tracked (otherwise 11.2-13.6 eV photons should all go here), but this is
/  not currently coded this way
/
/  RETURNS: FAIL or SUCCESS
/
*********************************kdissH2IINum***************************************/
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
#include "phys_constants.h"

int FindField(int f, int farray[], int n);

int grid::RadiativeTransferFUVPEHeating(PhotonPackageEntry **PP,
                                        FLOAT &dP_FUV,
                                        const int cellindex, const float tau,
                                        const FLOAT photonrate,
                                        float geo_correction, int FUVRateNum)
{

  // for individual stars, make sure below is consistent with
  // Star_ComputePhotonRates energies (should probably just make this a param)
   const double FUV_energy = FUV_photon_energy * erg_eV; // 5.6 - 11.2 eV average
   // be careful with this value !!! See Star_ComputePhotonRates
  dP_FUV = 0.0;
  // attenuation of the 5.6 - 11.2 eV band
  if (tau > 2.e1){
    dP_FUV = (1.0+BFLOAT_EPSILON) * (*PP)->Photons;
  } else if (tau > 1.0E-4){
    dP_FUV = min((*PP)->Photons*(1-expf(-tau)), ((*PP)->Photons));
  } else { // optically thin
    dP_FUV = min((*PP)->Photons*tau, (*PP)->Photons);
  }

  // dP_FUV is the number of absorptions due to dust
  // this is used to reduce photon count for next cell along ray
  dP_FUV = dP_FUV * geo_correction;

  // for this cell, however, we want to know ALL photons. This is because
  // PE heating rate calculation is computed using full FUV ISRF flux
  // convert photon count to luminosity to energy flux density
  const FLOAT dx2 = this->CellWidth[0][0] * this->CellWidth[0][0];
  // Need to multiply FUVRate field by EnergyUnits in Grid_FinalizeRadiationField
  BaryonField[FUVRateNum][cellindex] +=
                ((*PP)->Photons * photonrate * (FUV_energy) * geo_correction)
                / (dx2);

  return SUCCESS;
}
