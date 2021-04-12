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
/  modified2: Gen Chiaki
/             Add other molecules HD, CO, OH, and H2O
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
#include "phys_constants.h"

#define THRESHOLD_DENSITY_DB36 1e14
#define THRESHOLD_DENSITY_DB37 5e14
#define THRESHOLD_DENSITY_CO_CO   1e13
#define THRESHOLD_DENSITY_OH_OH   1e15
#define THRESHOLD_DENSITY_H2O_H2O 1e16
#define THRESHOLD_DENSITY_CO_H2   1e20
#define THRESHOLD_DENSITY_OH_H2   1e19

int grid::RadiativeTransferLWShielding(PhotonPackageEntry **PP, FLOAT &dP, 
				       FLOAT thisDensity, FLOAT thisDensityHDI, 
				       FLOAT thisDensityCO, FLOAT thisDensityOH, FLOAT thisDensityH2O, FLOAT ddr,
				       int cellindex, float LengthUnits, int kdissH2INum, int kdissHDINum,
				       int kdissCONum, int kdissOHNum, int kdissH2ONum,
				       int TemperatureField, float geo_correction)
{
  int H2Thin = 0;
  float shield1 = 0.0, shield2 = 0.0;
  float emission_dt_inv = 1.0 / (*PP)->EmissionTimeInterval;
  FLOAT dx = 0.0, dx2 = 0.0, dx3 = 0.0, Area_inv = 0.0;
  double cross_section = 0.0;
  double ColumnDensityPrim;

  float b = 0;                      //doppler parameter
  float b5 = 0;
  float x = 0;
  float alpha = 1.1;
  float H2mass = 2.0*mh; //grams

  float shield1HDI = 0.0, shield2HDI = 0.0;
  double cross_sectionHDI = 0.0;

  // fit from Heays et al. (2017) 
  int COThin, OHThin, H2OThin;
  float shield1CO  = 0.0, shield2CO  = 0.0;
  float shield1OH  = 0.0, shield2OH  = 0.0;
  float shield1H2O = 0.0, shield2H2O = 0.0;
  double cross_sectionCO = 0.0, cross_sectionOH = 0.0, cross_sectionH2O = 0.0;

  float N_CO_CO   = 9.48921e+14 ;//  +/- 3.906e+13    (4.116%)
  float a_CO_CO   = -0.822589   ;//  +/- 0.003422     (0.416%)
  float b_CO_CO   = -4.31591e-07;//  +/- 1.743e-08    (4.039%)
  float N_CO_H2   = 1.01533e+21 ;//  +/- 9.38e+19     (9.239%)
  float a_CO_H2   = -1.44968    ;//  +/- 0.06207      (4.282%)
  float b_CO_H2   = -0.171213   ;//  +/- 0.01413      (8.252%)

  float N_OH_OH   = 1.13614e+17 ;//  +/- 9.017e+15    (7.936%)
  float a_OH_OH   = -1.08514    ;//  +/- 0.01284      (1.183%)
  float b_OH_OH   = -0.00162352 ;//  +/- 0.0001285    (7.914%)
  float N_OH_H2   = 7.18749e+21 ;//  +/- 1.457e+21    (20.27%)
  float a_OH_H2   = -4.22208    ;//  +/- 0.6456       (15.29%)
  float b_OH_H2   = -0.901778   ;//  +/- 0.1155       (12.81%)

  float N_H2O_H2O = 1.02787e+17 ;//  +/- 1.365e+16    (13.28%)
  float a_H2O_H2O = -1.37727    ;//  +/- 0.02189      (1.589%)
  /* shielding by H2 is negligible for H2O */

  cross_section =  3.71e-18 * LengthUnits; // H2I average cross-section
  if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
    cross_sectionHDI =  4.22276e-18 * LengthUnits; // HDI average cross-section
  if(MetalChemistry && RadiativeTransferMetalDiss) {
    cross_sectionCO  =  1.74554e-17 * LengthUnits; // CO average cross-section
    cross_sectionOH  =  3.56213e-18 * LengthUnits; // OH average cross-section
    cross_sectionH2O =  1.01407e-17 * LengthUnits; // H2O average cross-section
  }
  dx = CellWidth[0][0];
  dx2 = dx*dx;
  dx3 = dx2*dx;
  Area_inv = ddr / dx3;
  //Calculate the dissociation rate per cell area
  //Units = cm^2*LengthUnits/(CodeLength^2*CodeTime)
  FLOAT dissrate = emission_dt_inv*cross_section*Area_inv; 

  FLOAT dissrateHDI = 0.0, dissrateCO = 0.0, dissrateOH = 0.0, dissrateH2O = 0.0;
  if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
    dissrateHDI = emission_dt_inv*cross_sectionHDI*Area_inv; 
  if(MetalChemistry && RadiativeTransferMetalDiss) {
    dissrateCO  = emission_dt_inv*cross_sectionCO *Area_inv; 
    dissrateOH  = emission_dt_inv*cross_sectionOH *Area_inv; 
    dissrateH2O = emission_dt_inv*cross_sectionH2O*Area_inv; 
  }
 
  ColumnDensityPrim = (*PP)->ColumnDensity; // H2 column density
  if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
    ColumnDensityPrim += (*PP)->ColumnDensityHDI; // HD column density
  if(RadiativeTransferH2ShieldType == 0) {
    /* We treat H2 dissociation with the shielding function from
       Draine & Bertoldi (1996) */
    if (ColumnDensityPrim < 0.01 * THRESHOLD_DENSITY_DB36) {
      shield1 = 1;
      if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
        shield1HDI = 1;
      H2Thin = TRUE;
    } else {
      shield1 = pow(ColumnDensityPrim / THRESHOLD_DENSITY_DB36, -0.75);
      if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
        shield1HDI = pow(ColumnDensityPrim / THRESHOLD_DENSITY_DB36, -0.75);
      H2Thin = FALSE;
    }
  }
  else if (RadiativeTransferH2ShieldType == 1) {
    /* We treat H2 dissociation with the shielding function from 
     * Equation 37 from Draine & Beltoldi with the exception that the 
     * power in the first term is 1.1 as per Wolcott-Green 2011
     */
    /*
     * The Wallcott-Green 2011 paper indicates that  a modified form
     * of the Draine & Bertoldi fitting formaula for H2 self shielding may be 
     * more appropriate 
     */
    b = sqrt(2.0*kboltz*BaryonField[TemperatureField][cellindex]/H2mass); // cm s^-1
    b5 = b/1e5;                                    // cm s^-1
    if (ColumnDensityPrim < 0.01 * THRESHOLD_DENSITY_DB37) {
      shield1 = 1;
      if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
        shield1HDI = 1;
      H2Thin = TRUE;
    }
    else {
      x = ColumnDensityPrim/(THRESHOLD_DENSITY_DB37);
      shield1 = 0.965/pow(1 + x/b5, alpha) + 
        (0.035/(pow(1 + x, 0.5)))*exp((-8.5e-4)*(pow(1 + x, 0.5)));
      if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
        shield1HDI = 0.965/pow(1 + x/(b5 * 0.8165), alpha) + 
          (0.035/(pow(1 + x, 0.5)))*exp((-8.5e-4)*(pow(1 + x, 0.5)));
            // multiply by sqrt(2/3)
      H2Thin = FALSE;
    }
//  printf("1 %13.5e %13.5e %13.5e %13.5e \n"
//     , BaryonField[TemperatureField][cellindex] 
//     , ColumnDensityPrim
//     , shield1, shield1HDI);
  }
  else
  {
    fprintf(stderr, "%s. Bad bad - it appears that we have no shielding " \
            "fitting function - abort\n", __FUNCTION__);
    return FAIL;
  }
  if(MetalChemistry && RadiativeTransferMetalDiss) {
    /* CO */
    if ((*PP)->ColumnDensityCO < THRESHOLD_DENSITY_CO_CO) {
      shield1CO = 1;
      COThin = TRUE;
    }
    else {
      x = (*PP)->ColumnDensityCO / N_CO_CO;
      shield1CO = pow(1.0 + x, a_CO_CO) * exp(b_CO_CO * x);
      COThin = FALSE;
    }
    if(MultiSpecies > 1) {
      if ((*PP)->ColumnDensity > THRESHOLD_DENSITY_CO_H2) {
        x = (*PP)->ColumnDensity / N_CO_H2;
        shield1CO *= pow(1.0 + x, a_CO_H2) * exp(b_CO_H2 * x);
        COThin = FALSE;
      }
    }

    /* OH */
    if ((*PP)->ColumnDensityOH < THRESHOLD_DENSITY_OH_OH) {
      shield1OH = 1;
      OHThin = TRUE;
    }
    else {
      x = (*PP)->ColumnDensityOH / N_OH_OH;
      shield1OH = pow(1.0 + x, a_OH_OH) * exp(b_OH_OH * x);
      OHThin = FALSE;
    }
    if(MultiSpecies > 1) {
      if ((*PP)->ColumnDensity > THRESHOLD_DENSITY_OH_H2) {
        x = (*PP)->ColumnDensity / N_OH_H2;
        shield1OH *= pow(1.0 + x, a_OH_H2) * exp(b_OH_H2 * x);
        OHThin = FALSE;
      }
    }

    /* H2O */
    if ((*PP)->ColumnDensityH2O < THRESHOLD_DENSITY_H2O_H2O) {
      shield1H2O = 1;
      H2OThin = TRUE;
    }
    else {
      x = (*PP)->ColumnDensityH2O / N_H2O_H2O;
      shield1H2O = pow(1.0 + x, a_H2O_H2O);
      H2OThin = FALSE;
    }
  }

  (*PP)->ColumnDensity += thisDensity * ddr * LengthUnits;
  if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
    (*PP)->ColumnDensityHDI += thisDensityHDI * ddr * LengthUnits;
  ColumnDensityPrim = (*PP)->ColumnDensity; // H2 column density
  if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
    ColumnDensityPrim += (*PP)->ColumnDensityHDI; // HD column density
  if(MetalChemistry && RadiativeTransferMetalDiss) {
    (*PP)->ColumnDensityCO  += thisDensityCO  * ddr * LengthUnits;
    (*PP)->ColumnDensityOH  += thisDensityOH  * ddr * LengthUnits;
    (*PP)->ColumnDensityH2O += thisDensityH2O * ddr * LengthUnits;
  }
  if(RadiativeTransferH2ShieldType == 0) {
    if (ColumnDensityPrim < 0.01 * THRESHOLD_DENSITY_DB36) {
      shield2 = 1;
      if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
        shield2HDI = 1;
    } else {
      shield2 = pow(ColumnDensityPrim / THRESHOLD_DENSITY_DB36, -0.75);
      if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
        shield2HDI = pow(ColumnDensityPrim / THRESHOLD_DENSITY_DB36, -0.75);
      H2Thin = FALSE;
    }
  }
  else if (RadiativeTransferH2ShieldType == 1) {
    if (ColumnDensityPrim < 0.01 * THRESHOLD_DENSITY_DB37) {
      shield2 = 1;
      if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
        shield2HDI = 1;
    } 
    else {
      x = ColumnDensityPrim/(THRESHOLD_DENSITY_DB37);
      shield2 = 0.965/pow(1 + x/b5, alpha) 
	+ (0.035/(pow(1 + x, 0.5)))*exp((-8.5e-4)*(pow(1 + x, 0.5)));
      if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
        shield2HDI = 0.965/pow(1 + x/(b5 * 0.8165), alpha) 
	  + (0.035/(pow(1 + x, 0.5)))*exp((-8.5e-4)*(pow(1 + x, 0.5)));
            // multiply by sqrt(2/3)
      H2Thin = FALSE;
    }
//  printf("2 %13.5e %13.5e %13.5e %13.5e \n"
//     , BaryonField[TemperatureField][cellindex] 
//     , ColumnDensityPrim
//     , shield2, shield2HDI);
  }
  else
  {
    fprintf(stderr, "%s. Bad bad - it appears that we have no shielding " \
            "fitting function - abort\n", __FUNCTION__);
    return FAIL;
  }
  if(MetalChemistry && RadiativeTransferMetalDiss) {
    /* CO */
    if ((*PP)->ColumnDensityCO < THRESHOLD_DENSITY_CO_CO) {
      shield2CO = 1;
      COThin = TRUE;
    }
    else {
      x = (*PP)->ColumnDensityCO / N_CO_CO;
      shield2CO = pow(1.0 + x, a_CO_CO) * exp(b_CO_CO * x);
      COThin = FALSE;
    }
    if(MultiSpecies > 1) {
      if ((*PP)->ColumnDensity > THRESHOLD_DENSITY_CO_H2) {
        x = (*PP)->ColumnDensity / N_CO_H2;
        shield2CO *= pow(1.0 + x, a_CO_H2) * exp(b_CO_H2 * x);
        COThin = FALSE;
      }
    }

    /* OH */
    if ((*PP)->ColumnDensityOH < THRESHOLD_DENSITY_OH_OH) {
      shield2OH = 1;
      OHThin = TRUE;
    }
    else {
      x = (*PP)->ColumnDensityOH / N_OH_OH;
      shield2OH = pow(1.0 + x, a_OH_OH) * exp(b_OH_OH * x);
      OHThin = FALSE;
    }
    if(MultiSpecies > 1) {
      if ((*PP)->ColumnDensity > THRESHOLD_DENSITY_OH_H2) {
        x = (*PP)->ColumnDensity / N_OH_H2;
        shield2OH *= pow(1.0 + x, a_OH_H2) * exp(b_OH_H2 * x);
        OHThin = FALSE;
      }
    }

    /* H2O */
    if ((*PP)->ColumnDensityH2O < THRESHOLD_DENSITY_H2O_H2O) {
      shield2H2O = 1;
      H2OThin = TRUE;
    }
    else {
      x = (*PP)->ColumnDensityH2O / N_H2O_H2O;
      shield2H2O = pow(1.0 + x, a_H2O_H2O);
      H2OThin = FALSE;
    }
  }


  if (H2Thin == TRUE) {
    dP = 0;
  }
  else {	
    dP = (*PP)->Photons * (1 - shield2/shield1);
    if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
      dP += (*PP)->Photons * (1 - shield2HDI/shield1HDI);
    if (MYPROC && DEBUG) 
      {
	fprintf(stdout, "%s: Column Density = %e\t dP(D&B) = %e\n", __FUNCTION__,
		ColumnDensityPrim, dP);
	fprintf(stdout, "%s: shield = %g\n",  __FUNCTION__, shield2);
      }
  }

  if(MetalChemistry && RadiativeTransferMetalDiss) {
    if (COThin == FALSE)
      dP += (*PP)->Photons * (1 - shield2CO /shield1CO );
    if (OHThin == FALSE)
      dP += (*PP)->Photons * (1 - shield2OH /shield1OH );
    if (H2OThin == FALSE)
      dP += (*PP)->Photons * (1 - shield2H2O/shield1H2O);
  }
  
  /* [geo_correction] = None (unitless)
   * [(*PP)->Photons] = 1/(LengthUnits^3)
   * [dissrate] = cm^2*CodeLength/(CodeLength^2*CodeTime)
   * Units = 1/(CodeTime) */
  
  BaryonField[kdissH2INum][cellindex] += geo_correction * (*PP)->Photons * 
    dissrate;
  if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
    BaryonField[kdissHDINum][cellindex] += geo_correction * (*PP)->Photons * 
      dissrateHDI;
  if(MetalChemistry && RadiativeTransferMetalDiss) {
    BaryonField[kdissCONum ][cellindex] += geo_correction * (*PP)->Photons * 
      dissrateCO;
    BaryonField[kdissOHNum ][cellindex] += geo_correction * (*PP)->Photons * 
      dissrateOH;
    BaryonField[kdissH2ONum][cellindex] += geo_correction * (*PP)->Photons * 
      dissrateH2O;
  }
  if(BaryonField[kdissH2INum][cellindex] < tiny_number)
    {
#if DEBUG
      fprintf(stdout, "Changing kdissH2I  from %g to %g\n", BaryonField[kdissH2INum][cellindex], tiny_number);
      fprintf(stdout, "(*PP)->Photons = %g\t shield2 = %g\t dP = %g\n", (*PP)->Photons, shield2, dP );
#endif
      BaryonField[kdissH2INum][cellindex] = tiny_number;
    }

  if(MultiSpecies > 2 && RadiativeTransferHDIDiss)
    if(BaryonField[kdissHDINum][cellindex] < tiny_number)
       BaryonField[kdissHDINum][cellindex] = tiny_number;
  if(MetalChemistry && RadiativeTransferMetalDiss) {
    if(BaryonField[kdissCONum][cellindex] < tiny_number)
       BaryonField[kdissCONum ][cellindex] = tiny_number;
    if(BaryonField[kdissOHNum][cellindex] < tiny_number)
       BaryonField[kdissOHNum ][cellindex] = tiny_number;
    if(BaryonField[kdissH2ONum][cellindex] < tiny_number)
       BaryonField[kdissH2ONum][cellindex] = tiny_number;
  }

//printf("%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e \n"
//   , (*PP)->ColumnDensity
//   , BaryonField[kdissH2INum][cellindex]
//   , (*PP)->ColumnDensityHDI
//   , BaryonField[kdissHDINum][cellindex]
//   , (*PP)->ColumnDensityCO
//   , BaryonField[kdissCONum][cellindex]
//   , (*PP)->ColumnDensityOH
//   , BaryonField[kdissOHNum][cellindex]
//   , (*PP)->ColumnDensityH2O
//   , BaryonField[kdissH2ONum][cellindex]);

  return SUCCESS;
}
