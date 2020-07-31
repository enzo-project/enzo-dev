/***********************************************************************
/
/  WRITES GRACKLE PARAMETERS FROM INPUT FILE
/
/  written by: Andrew Emerick
/  date:       April, 2020
/  modified1:
/
/  PURPOSE:
/
/  NOTE: Handles writing of Grackle-specific parameters to file
/
************************************************************************/

#include "preincludes.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/* function prototypes */

int GrackleWriteParameters(FILE *fptr)
{

#ifndef USE_GRACKLE

  if (use_grackle == TRUE){
    ENZO_FAIL("Error: Enzo must be compiled with 'make grackle-yes' to run with use_grackle = 1");
  }

#else

  fprintf(fptr, "with_radiative_cooling      = %d\n", grackle_data->with_radiative_cooling);
  fprintf(fptr, "use_volumetric_heating_rate = %d\n", grackle_data->use_volumetric_heating_rate);
  fprintf(fptr, "use_specific_heating_rate   = %d\n", grackle_data->use_specific_heating_rate);
  fprintf(fptr, "self_shielding_method       = %d\n", grackle_data->self_shielding_method);
  fprintf(fptr, "H2_self_shielding           = %d\n", grackle_data->H2_self_shielding);
  fprintf(fptr, "grackle_data_file           = %s\n", grackle_data->grackle_data_file);
  fprintf(fptr, "UVbackground                = %d\n", grackle_data->UVbackground);
  fprintf(fptr, "Compton_xray_heating        = %d\n", grackle_data->Compton_xray_heating);
  fprintf(fptr, "LWbackground_intensity      = %lf\n", grackle_data->LWbackground_intensity);
  fprintf(fptr, "LWbackground_sawtooth_suppression = %d\n", grackle_data->LWbackground_sawtooth_suppression);
  fprintf(fptr, "dust_chemistry              = %d\n",  grackle_data->dust_chemistry);
  fprintf(fptr, "local_dust_to_gas_ratio     = %lf\n", grackle_data->local_dust_to_gas_ratio);
  fprintf(fptr, "use_isrf_field              = %d\n",  grackle_data->use_isrf_field);
  fprintf(fptr, "use_dust_density_field      = %d\n",  grackle_data->use_dust_density_field);
#ifdef GRACKLE_MD
  fprintf(fptr, "UseDustDensityField         = %d\n", grackle_data->use_dust_density_field);
  fprintf(fptr, "MetalChemistry              = %d\n", grackle_data->metal_chemistry);
  fprintf(fptr, "GrainGrowth                 = %d\n", grackle_data->grain_growth);
  fprintf(fptr, "MetalAbundances             = %d\n", grackle_data->metal_abundances);
  fprintf(fptr, "DustSpecies                 = %d\n", grackle_data->dust_species);
  fprintf(fptr, "DustTemperatureMulti        = %d\n", grackle_data->dust_temperature_multi);
  fprintf(fptr, "DustSublimation             = %d\n", grackle_data->dust_sublimation);
#endif

#endif

  return SUCCESS;
}
