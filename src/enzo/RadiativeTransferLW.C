#define DEBUG 0
/***********************************************************************
/
/  Calculate the effect of Lyman Werner Radiation on the H2I Density
/  This function is based on the H2I cross section and is not fit
/  by a shielding function
/
/  written by: John Regan - based on the previous work of John Wise
/  date:     September 2014
/  modified1: 
/
/  PURPOSE: Calculate the dissociation rate due to Lyman Werner radiation
/
/  RETURNS: FAIL or SUCCESS
/
*********************************kdissH2INum***************************************/
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

int grid::RadiativeTransferLW(PhotonPackageEntry **PP, FLOAT &dPLW, int cellindex, 
			      float tau, FLOAT photonrate, 
			      float geo_correction, int kdissH2INum)
{
  // at most use all photons for photo-ionizations
  if (tau > 2.e1) //Completely Optically Thick
    dPLW = (*PP)->Photons;
  else if (tau > 1.e-4) //Exponential decline in photons
    dPLW = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
  else //Optically thin case
    dPLW = min((*PP)->Photons*tau, (*PP)->Photons);
  
  //dPLW is the number of absorptions due to H2I 
  dPLW = dPLW * geo_correction;
  
  // contributions to the photoionization rate is over whole timestep
  // Units = (1/CodeTime)*(1/LengthUnits^3)
  // BaryonField[kdissH2INum] needs to be normalised - see 
  // Grid_FinalizeRadiationFields.C
  BaryonField[kdissH2INum][cellindex] += dPLW*photonrate;
 
  return SUCCESS;
}
