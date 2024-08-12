/***********************************************************************
/
/  READS GRACKLE PARAMETERS FROM INPUT FILE
/
/  written by: Andrew Emerick
/  date:       December, 2019
/  modified1:
/
/  PURPOSE:
/
/  NOTE: Handles reading of Grackle-specific parameters and mapping
/        shared parameters between Enzo and Grackle. Any new Grackle parameters
/        should be read in here.
/
************************************************************************/

#include "preincludes.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"


int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int GrackleSetDefaultParameters(FILE *fptr){

  // Check if Grackle is being used. If so, copy over Grackle's
  // default parameters to their Enzo equivalents

  char line[MAX_LINE_LENGTH];

  /* First, check if Grackle is actually being used in this problem */
  rewind(fptr);
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    int ret = 0;
    ret += sscanf(line, "use_grackle = %"ISYM, &use_grackle);
  }

#ifndef USE_GRACKLE
  if (use_grackle == TRUE) {
    ENZO_FAIL("Error: Enzo must be compiled with 'make grackle-yes' to run with use_grackle = 1.\n");
  }
#else

  // Construct the Grackle chemistry data structure
  //
  chemistry_data *my_chemistry;
  my_chemistry = new chemistry_data;
  if (set_default_chemistry_parameters(my_chemistry) == FAIL) {
    ENZO_FAIL("Error in grackle: set_default_chemistry_parameters\n");
  }

  grackle_data->use_grackle = use_grackle;

  /* If we are actually using Grackle, overwrite Enzo's default parameters with
     the Grackle quivalents. This overrides what is set in SetDefaultGlobalValues */
  if (use_grackle){
    // Map Grackle defaults to corresponding Enzo parameters
    Gamma                                 = (float) grackle_data->Gamma;
    MultiSpecies                          = (int) grackle_data->primordial_chemistry;
    MetalCooling                          = (int) grackle_data->metal_cooling;
    H2FormationOnDust                     = (int) grackle_data->h2_on_dust;
    CloudyCoolingData.CMBTemperatureFloor = (int) grackle_data->cmb_temperature_floor;
    ThreeBodyRate                         = (int) grackle_data->three_body_rate;
    CIECooling                            = (int) grackle_data->cie_cooling;
    H2OpticalDepthApproximation           = (int) grackle_data->h2_optical_depth_approximation;
    PhotoelectricHeating                  = max((int) grackle_data->photoelectric_heating, 0);  // Grackle default is < 0
    PhotoelectricHeatingRate              = (float) grackle_data->photoelectric_heating_rate;
    CoolData.NumberOfTemperatureBins      = (int) grackle_data->NumberOfTemperatureBins;
    RateData.CaseBRecombination           = (int) grackle_data->CaseBRecombination;
    CoolData.TemperatureStart             = (float) grackle_data->TemperatureStart;
    CoolData.TemperatureEnd               = (float) grackle_data->TemperatureEnd;
    RateData.NumberOfDustTemperatureBins  = (int) grackle_data->NumberOfDustTemperatureBins;
    RateData.DustTemperatureStart         = (float) grackle_data->DustTemperatureStart;
    RateData.DustTemperatureEnd           = (float) grackle_data->DustTemperatureEnd;
    CoolData.HydrogenFractionByMass       = (float) grackle_data->HydrogenFractionByMass;
    CoolData.DeuteriumToHydrogenRatio     = (float) grackle_data->DeuteriumToHydrogenRatio;
    CoolData.SolarMetalFractionByMass     = (float) grackle_data->SolarMetalFractionByMass;
  } // end use_grackle

#endif // end USE_GRACKLE

  return SUCCESS;
}

int GrackleReadParameters(FILE *fptr, FLOAT InitTime)
{
  /* Read Grackle specific parameters and copy over Enzo
     parameters to their Grackle equivalents. Concludes by
     setting up grackle_units */

  if (use_grackle == FALSE) {
    return SUCCESS;
  }

  char line[MAX_LINE_LENGTH];
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;


#ifdef USE_GRACKLE

  // Go back through parameter file to check for Grackle-specific
  // parameters that do not have Enzo equivalents
  rewind(fptr);
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    int ret = 0;
    ret += sscanf(line, "with_radiative_cooling = %d",
                  &grackle_data->with_radiative_cooling);
    ret += sscanf(line, "use_volumetric_heating_rate = %d",
                    &grackle_data->use_volumetric_heating_rate);
    ret += sscanf(line, "use_specific_heating_rate = %d",
                    &grackle_data->use_specific_heating_rate);
    ret += sscanf(line, "self_shielding_method = %d",
                    &grackle_data->self_shielding_method);
    ret += sscanf(line, "H2_self_shielding = %d",
                    &grackle_data->H2_self_shielding);

    if (sscanf(line, "grackle_data_file = %s", dummy) == 1) {
      grackle_data->grackle_data_file = dummy;
      ret++;
    }
    ret += sscanf(line, "UVbackground = %d", &grackle_data->UVbackground);
    ret += sscanf(line, "Compton_xray_heating = %d",
                    &grackle_data->Compton_xray_heating);
    ret += sscanf(line, "LWbackground_intensity = %lf",
                  &grackle_data->LWbackground_intensity);
    ret += sscanf(line, "LWbackground_sawtooth_suppression = %d",
                  &grackle_data->LWbackground_sawtooth_suppression);

    ret += sscanf(line, "local_dust_to_gas_ratio = %f",
                  &grackle_data->local_dust_to_gas_ratio);
    ret += sscanf(line, "interstellar_radiation_field = %f",
                  &grackle_data->interstellar_radiation_field);
    ret += sscanf(line, "dust_recombination_cooling = %d",
                  &grackle_data->dust_recombination_cooling);

    ret += sscanf(line, "dust_chemistry = %d",
                  &grackle_data->dust_chemistry);

    /* functionality for below two are not yet implemented but are
       involved in options for other Grackle settings. Read in
       here to do error checking to make sure these are not used */
    ret += sscanf(line, "use_isrf_field = %d",
                  &grackle_data->use_isrf_field);
    ret += sscanf(line, "use_dust_density_field = %d",
                  &grackle_data->use_dust_density_field);

    /* If the dummy char space was used, then make another. */
    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      dummy[0] = 0;
      ret++;
    }

  }
  /* clean up */
  rewind(fptr);

  // Finally, map all Enzo parameters to their Grackle equivalents. These
  // are read in ReadParameterFile
  //
  //
  // grackle_data->with_radiative_cooling already set
  // grackle_data->grackle_data_file already set
  // grackle_data->UVbackground already set
  // grackle_data->Compton_xray_heating already set
  // grackle_data->LWbackground_intensity already set
  // grackle_data->LWbackground_sawtooth_suppression already set
  grackle_data->use_grackle                    = (Eint32) use_grackle;
  grackle_data->Gamma                          = (double) Gamma;
  grackle_data->primordial_chemistry           = (Eint32) MultiSpecies;
  grackle_data->metal_cooling                  = (Eint32) MetalCooling;
  grackle_data->h2_on_dust                     = (Eint32) H2FormationOnDust;
  grackle_data->cmb_temperature_floor          = (Eint32) CloudyCoolingData.CMBTemperatureFloor;
  grackle_data->three_body_rate                = (Eint32) ThreeBodyRate;
  grackle_data->cie_cooling                    = (Eint32) CIECooling;
  grackle_data->h2_optical_depth_approximation = (Eint32) H2OpticalDepthApproximation;
  grackle_data->photoelectric_heating          = (Eint32) PhotoelectricHeating;
  grackle_data->photoelectric_heating_rate     = (double) PhotoelectricHeatingRate;
  grackle_data->NumberOfTemperatureBins        = (Eint32) CoolData.NumberOfTemperatureBins;
  grackle_data->CaseBRecombination             = (Eint32) RateData.CaseBRecombination;
  grackle_data->TemperatureStart               = (double) CoolData.TemperatureStart;
  grackle_data->TemperatureEnd                 = (double) CoolData.TemperatureEnd;
  grackle_data->NumberOfDustTemperatureBins    = (Eint32) RateData.NumberOfDustTemperatureBins;
  grackle_data->DustTemperatureStart           = (double) RateData.DustTemperatureStart;
  grackle_data->DustTemperatureEnd             = (double) RateData.DustTemperatureEnd;
  grackle_data->HydrogenFractionByMass         = (double) CoolData.HydrogenFractionByMass;
  grackle_data->DeuteriumToHydrogenRatio       = (double) CoolData.DeuteriumToHydrogenRatio;
  grackle_data->SolarMetalFractionByMass       = (double) CoolData.SolarMetalFractionByMass;
  grackle_data->UVbackground_redshift_on       = (double) CoolData.RadiationRedshiftOn;
  grackle_data->UVbackground_redshift_off      = (double) CoolData.RadiationRedshiftOff;
  grackle_data->UVbackground_redshift_fullon   = (double) CoolData.RadiationRedshiftFullOn;
  grackle_data->UVbackground_redshift_drop     = (double) CoolData.RadiationRedshiftDropOff;
  grackle_data->use_radiative_transfer         = (Eint32) RadiativeTransfer;
  // grackle_data->radiative_transfer_coupled_rate_solver set in RadiativeTransferReadParameters
  // grackle_data->radiative_transfer_hydrogen_only set in RadiativeTransferReadParameters


  // Error checking for behavior not implemented
  if ( (grackle_data->photoelectric_heating == 2) ||
       (grackle_data->use_isrf_field)){
    ENZO_FAIL("Photoelectric heating model 2, and ISRF field, in Grackle is not yet implemented.\n");
  }

  if ( grackle_data->use_dust_density_field ){
    ENZO_FAIL("Supplying dust density (use_dust_density_field) to Grackle is not yet implemented.\n");
  }

  // Initialize Grackle units structure.
  FLOAT a_value, dadt;
  a_value = 1.0;
  code_units grackle_units;
  grackle_units.a_units = 1.0;
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0,
    TimeUnits = 1.0, VelocityUnits = 1.0;
  double MassUnits = 1.0;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, &MassUnits, InitTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(InitTime, &a_value,
                                        &dadt) == FAIL) {
      ENZO_FAIL("Error in CosmologyComputeExpansionFactors.\n");
    }
    grackle_units.a_units            = (double) (1.0 / (1.0 + InitialRedshift));
  }
  grackle_units.comoving_coordinates = (Eint32) ComovingCoordinates;
  grackle_units.density_units        = (double) DensityUnits;
  grackle_units.length_units         = (double) LengthUnits;
  grackle_units.time_units           = (double) TimeUnits;
  grackle_units.velocity_units       = (double) VelocityUnits;
  grackle_units.a_value              = (double) a_value;

  // Initialize chemistry structure.
  if (initialize_chemistry_data(&grackle_units) == FAIL) {
    ENZO_FAIL("Error in Grackle initialize_chemistry_data.\n");
  }

  // Need to set these after initialize_chemistry_data since
  // that function sets them automatically based on the tables.
  if (FinalRedshift < grackle_data->UVbackground_redshift_off) {
    grackle_data->UVbackground_redshift_off = FinalRedshift;
    grackle_data->UVbackground_redshift_drop = FinalRedshift;
  }

#endif // end USE_GRACKLE

  delete [] dummy;

  return SUCCESS;
}
