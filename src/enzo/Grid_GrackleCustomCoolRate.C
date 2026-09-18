/***********************************************************************
/
/  GRID CLASS (COMPUTE THE COOLING TIME FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:  Elizabeth Harper-Clark, August 2009
/              added in CoolingModel parameter
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
// Compute the cooling time

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
 
/* This parameter controls whether the cooling function recomputes
   the metal cooling rates.  It is reset by RadiationFieldUpdate. */
 
extern int RadiationFieldRecomputeMetalRates;
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RadiationFieldCalculateRates(FLOAT Time);
int FindField(int field, int farray[], int numfields);

int grid::GrackleCustomCoolRate(int rank, int *dim, float *cool_rate,
				float *dens, float *thrmeng,
				float *velx, float *vely, float *velz,
				float *HIdens, float *HIIdens,
				float *HeIdens, float *HeIIdens, float *HeIIIdens,
				float *edens,
				float *HMdens, float *H2Idens, float *H2IIdens,
				float *DIdens, float *DIIdens, float *HDIdens,
				float *metaldens,
				float *kphHI, float *kphHeI, float *kphHeII,
				float *kdissH2I, float *gamma)
{

  // All passed fields should be in code units

#ifdef USE_GRACKLE 

  /* Return if this doesn't concern us. */

  if (grackle_data->use_grackle == FALSE) return SUCCESS;
  
  if (RadiativeCooling == 0) return SUCCESS;
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i;
  int size = 1;
  for (int d = 0; d < rank; d++)
    size *= dim[d];
 
  /* Compute the cooling time. */
 
  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
 
    aUnits = 1.0/(1.0 + InitialRedshift);
  } else if (RadiationFieldRedshift > -1){
    a       = 1.0 / (1.0 + RadiationFieldRedshift);
    aUnits  = 1.0;
  }
  float afloat = float(a);

  float *volumetric_heating_rate = NULL;
  float *specific_heating_rate   = NULL;

  // Double check if there's a metal field when we have metal cooling
  int metal_cooling = MetalCooling;
  if (metal_cooling && !metaldens) {
    if (debug)
      fprintf(stderr, "Warning: No metal field passed to GrackleCustomCoolRate. Not using metal cooling.\n");
    metal_cooling = FALSE;
  }

  Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
  g_grid_dimension = new Eint32[3];
  g_grid_start = new Eint32[3];
  g_grid_end = new Eint32[3];

  // Fortran code will act as if there are 3 dimensions regardless...
  for (i = 0; i < rank; i++) {
    g_grid_dimension[i] = (Eint32) dim[i];
    g_grid_start[i] = (Eint32) 0;
    g_grid_end[i] = (Eint32) dim[i]-1;
  }
  // ...so for any unused dimensions, set quantities to 0.
  // will do nothing if rank == 2 (3-dimensional problem)
  for (i = rank; i < 3; i++){
    g_grid_dimension[i] = (Eint32) 0;
    g_grid_start[i] = (Eint32) 0;
    g_grid_end[i] = (Eint32) 0;
  }

  /* Update units. */

  code_units grackle_units;
  grackle_units.comoving_coordinates = (Eint32) ComovingCoordinates;
  grackle_units.density_units        = (double) DensityUnits;
  grackle_units.length_units         = (double) LengthUnits;
  grackle_units.time_units           = (double) TimeUnits;
  grackle_units.velocity_units       = (double) VelocityUnits;
  grackle_units.a_units              = (double) aUnits;
  grackle_units.a_value              = (double) a;

  /* set up the my_fields */
  grackle_field_data my_fields;

  my_fields.grid_rank = (Eint32) rank;
  my_fields.grid_dimension = g_grid_dimension;
  my_fields.grid_start     = g_grid_start;
  my_fields.grid_end       = g_grid_end;
  my_fields.grid_dx        = this->CellWidth[0][0]; // CHANGE

  GrackleFieldBuffer buf_dens, buf_te, buf_v1, buf_v2, buf_v3;
  GrackleFieldBuffer buf_HI, buf_HII, buf_HeI, buf_HeII, buf_HeIII, buf_e;
  GrackleFieldBuffer buf_HM, buf_H2I, buf_H2II, buf_DI, buf_DII, buf_HDI;
  GrackleFieldBuffer buf_metal, buf_vol_heat, buf_spec_heat;
  GrackleFieldBuffer buf_cool_rate;
#ifdef TRANSFER
  GrackleFieldBuffer buf_rt_hi, buf_rt_hei, buf_rt_heii, buf_rt_h2, buf_rt_heat;
#endif

  /* now add in the baryon fields */
  my_fields.density         = buf_dens.prepare(dens, size);
  my_fields.internal_energy = buf_te.prepare(thrmeng, size);
  my_fields.x_velocity      = buf_v1.prepare(velx, size);
  my_fields.y_velocity      = buf_v2.prepare(vely, size);
  my_fields.z_velocity      = buf_v3.prepare(velz, size);

  if (MultiSpecies) {
    my_fields.HI_density      = buf_HI.prepare(HIdens, size);
    my_fields.HII_density     = buf_HII.prepare(HIIdens, size);
    my_fields.HeI_density     = buf_HeI.prepare(HeIdens, size);
    my_fields.HeII_density    = buf_HeII.prepare(HeIIdens, size);
    my_fields.HeIII_density   = buf_HeIII.prepare(HeIIIdens, size);
    my_fields.e_density       = buf_e.prepare(edens, size);

    if (MultiSpecies > 1) {
      my_fields.HM_density      = buf_HM.prepare(HMdens, size);
      my_fields.H2I_density     = buf_H2I.prepare(H2Idens, size);
      my_fields.H2II_density    = buf_H2II.prepare(H2IIdens, size);
    
      if (MultiSpecies > 2) {
	my_fields.DI_density      = buf_DI.prepare(DIdens, size);
	my_fields.DII_density     = buf_DII.prepare(DIIdens, size);
	my_fields.HDI_density     = buf_HDI.prepare(HDIdens, size);
      }
    }
  }
  
  my_fields.metal_density   = buf_metal.prepare(metaldens, size);
  
  my_fields.volumetric_heating_rate  = buf_vol_heat.prepare(volumetric_heating_rate, size);
  my_fields.specific_heating_rate    = buf_spec_heat.prepare(specific_heating_rate, size);

#ifdef TRANSFER

  /* unit conversion from Enzo RT units to CGS */
  const float ev2erg = 1.60217653E-12;
  float rtunits = ev2erg / TimeUnits;

  if ( RadiativeTransfer ){
    my_fields.RT_HI_ionization_rate = buf_rt_hi.prepare(kphHI, size);

    if (RadiativeTransferHydrogenOnly == FALSE){
      my_fields.RT_HeI_ionization_rate  = buf_rt_hei.prepare(kphHeI, size);
      my_fields.RT_HeII_ionization_rate = buf_rt_heii.prepare(kphHeII, size);
    }

    if (MultiSpecies > 1)
      my_fields.RT_H2_dissociation_rate = buf_rt_h2.prepare(kdissH2I, size);

    /* need to convert to CGS units */
    for( i = 0; i < size; i++) gamma[i] *= rtunits;

    my_fields.RT_heating_rate = buf_rt_heat.prepare(gamma, size);

  }
#endif // TRANSFER

  gr_float *cool_rate_ptr = buf_cool_rate.prepare(cool_rate, size, false);

  if (calculate_cooling_time(&grackle_units, &my_fields, cool_rate_ptr) == FAIL) {
    ENZO_FAIL("Error in Grackle calculate_cooling_time.\n");
  }
  buf_cool_rate.copy_back();

  // Code units
  for (i = 0; i < size; i++) {
    cool_rate[i] = thrmeng[i] / fabs(cool_rate[i]) / dens[i];
  }
    
#ifdef TRANSFER
  if (RadiativeTransfer){
    /* convert the RT units back to Enzo */
    for(i = 0; i < size; i ++) gamma[i] /= rtunits;

  }
#endif // TRANSFER

  delete [] g_grid_dimension;
  delete [] g_grid_start;
  delete [] g_grid_end;

#else

    printf("WARNING: Calling GrackleCustomCoolRate but USE_GRACKLE is False!\n");
    
#endif // USE_GRACKLE
    
    return SUCCESS;
}



