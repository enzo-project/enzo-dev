#define DEBUG 0
#define MYPROC MyProcessorNumber == ProcessorNumber
/***********************************************************************
/
/  Calculate the effect of LW radiation on H2 molecules
/
/  written by: John Wise
/  date:     
/  modified1: John Regan
/             Moved to its own function and file from
/             WalkPhotonPackage
/
/  PURPOSE: Calculate the photodissociation rate 
/           of H2 due to radiation in the Lyman Werner band. 
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

#define THRESHOLD_DENSITY_DB36 1e14
#define THRESHOLD_DENSITY_DB37 5e14

int grid::RadiativeTransferLWShielding(PhotonPackageEntry **PP, FLOAT &dP, 
				       FLOAT thisDensity, FLOAT ddr,
				       int cellindex, float LengthUnits, int kdissH2INum, 
				       int TemperatureField, float geo_correction)
{
  int H2Thin = 0;
  float shield1 = 0.0, shield2 = 0.0;
  float emission_dt_inv = 1.0 / (*PP)->EmissionTimeInterval;
  FLOAT dx = 0.0, dx2 = 0.0, dx3 = 0.0, Area_inv = 0.0;
  double cross_section = 0.0;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				  HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				  DIINum, HDINum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }

  cross_section =  3.71e-18 * LengthUnits; // H2I average cross-section
  dx = CellWidth[0][0];
  dx2 = dx*dx;
  dx3 = dx2*dx;
  Area_inv = ddr / dx3;
  //Calculate the dissociation rate per cell area
  //Units = cm^2*LengthUnits/(CodeLength^2*CodeTime)
  FLOAT dissrate = emission_dt_inv*cross_section*Area_inv; 
 
  if(RadiativeTransferH2ShieldType == 0) {
    /* We treat H2 dissociation with the shielding function from
       Draine & Bertoldi (1996) */
    if ((*PP)->ColumnDensity < THRESHOLD_DENSITY_DB36) {
      shield1 = 1;
      H2Thin = TRUE;
    } else {
      shield1 = pow((*PP)->ColumnDensity / THRESHOLD_DENSITY_DB36, -0.75);
      H2Thin = FALSE;
    }
    
    (*PP)->ColumnDensity += thisDensity * ddr * LengthUnits;
    if ((*PP)->ColumnDensity < THRESHOLD_DENSITY_DB36) {
      shield2 = 1;
    } else {
      shield2 = pow((*PP)->ColumnDensity / THRESHOLD_DENSITY_DB36, -0.75);
      H2Thin = FALSE;
    }
  }
  else if (RadiativeTransferH2ShieldType == 1) {
    /* We treat H2 dissociation with the shielding function from 
     * Equation 37 from Draine & Beltoldi with the exception that the 
     * power in the first term is 1.1 as per Wolcott-Green 2011
     */
    float b = 0;                      //doppler parameter
    float b5 = 0;
    float x = 0;
    float alpha = 1.1;
    float kb = 1.3806504e-16;        //erg K^-1
    float H2mass = 2.0*1.672623e-24; //grams
    float shield1_db = 1.0, shield2_db = 1.0;
    /*
     * The Wallcott-Green 2011 paper indicates that  a modified form
     * of the Draine & Bertoldi fitting formaula for H2 self shielding may be 
     * more appropriate 
     */
   
    
    b = sqrt(2.0*kb*BaryonField[TemperatureField][cellindex]/H2mass); // cm s^-1
    b5 = b/1e5;                                    // cm s^-1
    if ((*PP)->ColumnDensity < THRESHOLD_DENSITY_DB37) {
      shield1 = 1;
      H2Thin = TRUE;
    }
    else {
      x = (*PP)->ColumnDensity/(THRESHOLD_DENSITY_DB37);
      shield1 = 0.965/pow(1 + x/b5, alpha) + 
	(0.035/(pow(1 + x, 0.5)))*exp((-8.5e-4)*(pow(1 + x, 0.5)));
      H2Thin = FALSE;
    }
    
    (*PP)->ColumnDensity += thisDensity * ddr * LengthUnits;
    if ((*PP)->ColumnDensity < THRESHOLD_DENSITY_DB37) {
      shield2 = 1;
    } 
    else {
      x = (*PP)->ColumnDensity/(THRESHOLD_DENSITY_DB37);
      shield2 = 0.965/pow(1 + x/b5, alpha) 
	+ (0.035/(pow(1 + x, 0.5)))*exp((-8.5e-4)*(pow(1 + x, 0.5)));
      H2Thin = FALSE;
    }
  }
  else
    {
      fprintf(stderr, "%s. Bad bad - it appears that we have no shielding " \
	      "fitting function - abort\n", __FUNCTION__);
      return FAIL;
    }
  
  if (H2Thin == TRUE) {
    dP = 0;
  }
  else {	
    dP = (*PP)->Photons * (1 - shield2/shield1);
    if (MYPROC && DEBUG) 
      {
	fprintf(stdout, "%s: Column Density = %e\t dP(D&B) = %e\n", __FUNCTION__,
		(*PP)->ColumnDensity, dP);
	fprintf(stdout, "%s: shield = %g\n",  __FUNCTION__, shield2);
      }
  }
  
  /* [geo_correction] = None (unitless)
   * [(*PP)->Photons] = 1/(LengthUnits^3)
   * [dissrate] = cm^2*CodeLength/(CodeLength^2*CodeTime)
   /* Units = 1/(CodeTime) */
  
  BaryonField[kdissH2INum][cellindex] += geo_correction * (*PP)->Photons * 
    dissrate;
   if(BaryonField[kdissH2INum][cellindex] < tiny_number)
    {
#if DEBUG
      fprintf(stdout, "Changing kdissH2I  from %g to %g\n", BaryonField[kdissH2INum][cellindex], tiny_number);
      fprintf(stdout, "(*PP)->Photons = %g\t shield2 = %g\t dP = %g\n", (*PP)->Photons, shield2, dP );
#endif
      BaryonField[kdissH2INum][cellindex] = tiny_number;
    }
      
  return SUCCESS;
}
