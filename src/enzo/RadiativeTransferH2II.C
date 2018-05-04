#define DEBUG 0
/***********************************************************************
/
/  Calculate the effect of Lyman Werner Radiation on the H2II Density
/  The cross section is taken from Stancil et al. (1994)
/
/  written by: John Regan - based on the previous work of John Wise
/  date:     April 2015
/  modified1: 
/
/  PURPOSE: Calculate the dissociation rate due to Lyman Werner radiation
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

int grid::RadiativeTransferH2II(PhotonPackageEntry **PP, int cellindex, 
				float tau, FLOAT photonrate, float geo_correction,
				int kdissH2IINum)
{
  FLOAT dPH2II = 0.0;
  // at most use all photons for photo-ionizations
  if (tau > 2.e1) //Completely Optically Thick
    dPH2II = (1.0+BFLOAT_EPSILON) * (*PP)->Photons;
  else if (tau > 1.e-4) //Exponential decline in photons
    dPH2II = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
  else //Optically thin case
    dPH2II = min((*PP)->Photons*tau, (*PP)->Photons);
  
  //dPH2II is the number of absorptions due to H2II
  dPH2II = dPH2II * geo_correction;
  // contributions to the photoionization rate is over whole timestep
  // Units = (1/CodeTime)*(1/LengthUnits^3)
  // BaryonField[kdissH2IINum] needs to be normalised - see 
  // Grid_FinalizeRadiationFields.C
  BaryonField[kdissH2IINum][cellindex] += dPH2II*photonrate;
  if(BaryonField[kdissH2IINum][cellindex] < tiny_number)
    {
      BaryonField[kdissH2IINum][cellindex] = tiny_number;
    }
      
  return SUCCESS;
}
