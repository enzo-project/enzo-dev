/***********************************************************************
/
/  Calculate the effect of Infrared radiation on the H^- fraction 
/
/  written by: John Regan - based on the previous work of John Wise
/  date:     September 2014
/  modified1: 
/
/  PURPOSE: Calculate the HM destruction rate
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

int grid::RadiativeTransferIR(PhotonPackageEntry **PP, FLOAT &dPIR, int cellindex, 
			      float tau, FLOAT photonrate, 
			      FLOAT *excessrate, float geo_correction,
			      int kphHMNum, int gammaNum)
{

  // at most use all photons for photo-ionizations
  if (tau > 2.e1) //Completely Optically Thick
    dPIR = (1.0+BFLOAT_EPSILON) * (*PP)->Photons;
  else if (tau > 1.e-4) //Exponential decline in photons
    dPIR = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
  else //Optically thin case
    dPIR = min((*PP)->Photons*tau, (*PP)->Photons);
  
  //dPIR is the number of absorptions
  dPIR = dPIR * geo_correction;

  // contributions to the photoionization rate is over whole timestep
  // Units = (1/CodeTime)*(1/LengthUnits**3)
  BaryonField[kphHMNum][cellindex] += dPIR*photonrate;
  // the heating rate is just the number of photo ionizations (Units = (1/LengthUnits**3))
  // times the excess energy units here are eV/CodeTime.
  // Units = Ev per time per LengthUnits^3 [Ev/CodeTime/LengthUnits**3]
  BaryonField[gammaNum][cellindex] += dPIR*excessrate[IR];
  
  return SUCCESS;
}
