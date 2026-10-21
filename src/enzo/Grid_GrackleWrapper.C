/***********************************************************************
/
/  GRID CLASS (WRAP THE GRACKLE CHEMISTRY SOLVER)
/
/  written by: Britton Smith
/  date:       April, 2013
/  modified1:
/
/  PURPOSE: Solve chemistry and cooling with grackle.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "preincludes.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);

int grid::GrackleWrapper()
{

#ifdef USE_GRACKLE
  if (grackle_data->use_grackle == FALSE)
    return SUCCESS;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  LCAPERF_START("grid_GrackleWrapper");

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  double dt_cool = dtFixed;
#ifdef TRANSFER
  dt_cool = (grackle_data->radiative_transfer_intermediate_step == TRUE) ? dtPhoton : dtFixed;
#endif
  
  /* Compute the size of the fields. */
 
  int i;
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  Eint32 *g_grid_dimension, *g_grid_start, *g_grid_end;
  g_grid_dimension = new Eint32[GridRank];
  g_grid_start = new Eint32[GridRank];
  g_grid_end = new Eint32[GridRank];
  for (i = 0; i < GridRank; i++) {
    g_grid_dimension[i] = (Eint32) GridDimension[i];
    g_grid_start[i] = (Eint32) GridStartIndex[i];
    g_grid_end[i] = (Eint32) GridEndIndex[i];
  }
 
  /* Find fields: density, total energy, velocity1-3. */
 
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* Find Multi-species fields. */

  DeNum = HINum = HIINum = HeINum = HeIINum = HeIIINum = HMNum = H2INum = 
    H2IINum = DINum = DIINum = HDINum = 0;
 
  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }
 
  /* Get easy to handle pointers for each variable. */
 
  float *density     = BaryonField[DensNum];
  float *totalenergy = BaryonField[TENum];
  float *gasenergy   = BaryonField[GENum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];

  float *volumetric_heating_rate = NULL;
  float *specific_heating_rate   = NULL;

  /* Compute the cooling time. */

  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time+0.5*dt_cool, &a, &dadt);
 
    aUnits = 1.0/(1.0 + InitialRedshift);
  } else if (RadiationFieldRedshift > -1){
    a        = 1.0 / (1.0 + RadiationFieldRedshift);
    aUnits   = 1.0;
  }
  float afloat = float(a);

  /* Update units. */

  code_units grackle_units;
  grackle_units.comoving_coordinates = (Eint32) ComovingCoordinates;
  grackle_units.density_units        = (double) DensityUnits;
  grackle_units.length_units         = (double) LengthUnits;
  grackle_units.time_units           = (double) TimeUnits;
  grackle_units.velocity_units       = (double) VelocityUnits;
  grackle_units.a_units              = (double) aUnits;
  grackle_units.a_value              = (double) a;

  /* Metal cooling codes. */
 
  int MetalNum = 0, SNColourNum = 0;
  int MetalFieldPresent = FALSE;

  // First see if there's a metal field (so we can conserve species in
  // the solver)
  MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields);
  SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
  MetalFieldPresent = (MetalNum != -1 || SNColourNum != -1);

  // Double check if there's a metal field when we have metal cooling
  if (MetalCooling && MetalFieldPresent == FALSE) {
    ENZO_FAIL("Metal cooling is on, but no metal field present.");
  }

  /* If both metal fields (Pop I/II and III) exist, create a field
     that contains their sum */

  float *MetalPointer = NULL;
  float *TotalMetals = NULL;

  if (MetalNum != -1 && SNColourNum != -1) {
    TotalMetals = new float[size];
    for (i = 0; i < size; i++)
      TotalMetals[i] = BaryonField[MetalNum][i] + BaryonField[SNColourNum][i];
    MetalPointer = TotalMetals;
  } // ENDIF both metal types
  else {
    if (MetalNum != -1)
      MetalPointer = BaryonField[MetalNum];
    else if (SNColourNum != -1)
      MetalPointer = BaryonField[SNColourNum];
  } // ENDELSE both metal types

  /*
    Cap the metal density at 90% of the gas density.
    This was historically in Enzo, then Grackle and was removed in
    Grackle PR #121 (https://github.com/grackle-project/grackle/pull/121)
  */
  if (MetalPointer != NULL) {
    for (i = 0; i < size; i++) {
      MetalPointer[i] = min(MetalPointer[i], 0.9 * BaryonField[DensNum][i]);
    }
  }
 
  int temp_thermal = FALSE;
  float *thermal_energy;
  if ( UseMHD ){
    iBx = FindField(Bfield1, FieldType, NumberOfBaryonFields);
    iBy = FindField(Bfield2, FieldType, NumberOfBaryonFields);
    iBz = FindField(Bfield3, FieldType, NumberOfBaryonFields);  
  }

  if (HydroMethod==Zeus_Hydro) {
    thermal_energy = BaryonField[TENum];
  }
  else if (DualEnergyFormalism) {
    thermal_energy = BaryonField[GENum];
  }
  else {
    temp_thermal = TRUE;
    thermal_energy = new float[size];
    for (i = 0; i < size; i++) {
      thermal_energy[i] = BaryonField[TENum][i] - 
        0.5 * POW(BaryonField[Vel1Num][i], 2.0);
      if(GridRank > 1)
        thermal_energy[i] -= 0.5 * POW(BaryonField[Vel2Num][i], 2.0);
      if(GridRank > 2)
        thermal_energy[i] -= 0.5 * POW(BaryonField[Vel3Num][i], 2.0);

      if( UseMHD ) {
        thermal_energy[i] -= 0.5 * (POW(BaryonField[iBx][i], 2.0) + 
                                    POW(BaryonField[iBy][i], 2.0) + 
                                    POW(BaryonField[iBz][i], 2.0)) / 
          BaryonField[DensNum][i];
      }
    } // for (int i = 0; i < size; i++)
  }

  //
  // Put code here to assign fields to volumetric or specific
  // heating rate pointers
  //

  /* set up grackle fields object */
  grackle_field_data my_fields;

  my_fields.grid_rank = (Eint32) GridRank;
  my_fields.grid_dimension = g_grid_dimension;
  my_fields.grid_start     = g_grid_start;
  my_fields.grid_end       = g_grid_end;
  my_fields.grid_dx        = this->CellWidth[0][0];

  GrackleFieldBuffer buf_dens, buf_te, buf_v1, buf_v2, buf_v3;
  GrackleFieldBuffer buf_HI, buf_HII, buf_HeI, buf_HeII, buf_HeIII, buf_e;
  GrackleFieldBuffer buf_HM, buf_H2I, buf_H2II, buf_DI, buf_DII, buf_HDI;
  GrackleFieldBuffer buf_metal, buf_vol_heat, buf_spec_heat;
#ifdef TRANSFER
  GrackleFieldBuffer buf_rt_hi, buf_rt_hei, buf_rt_heii, buf_rt_h2, buf_rt_heat;
#endif

  /* now add in the baryon fields */
  my_fields.density         = buf_dens.prepare(density, size);
  my_fields.internal_energy = buf_te.prepare(thermal_energy, size);
  my_fields.x_velocity      = buf_v1.prepare(velocity1, size);
  my_fields.y_velocity      = buf_v2.prepare(velocity2, size);
  my_fields.z_velocity      = buf_v3.prepare(velocity3, size);
  my_fields.HI_density      = buf_HI.prepare(BaryonField[HINum], size);
  my_fields.HII_density     = buf_HII.prepare(BaryonField[HIINum], size);
  my_fields.HeI_density     = buf_HeI.prepare(BaryonField[HeINum], size);
  my_fields.HeII_density    = buf_HeII.prepare(BaryonField[HeIINum], size);
  my_fields.HeIII_density   = buf_HeIII.prepare(BaryonField[HeIIINum], size);
  my_fields.e_density       = buf_e.prepare(BaryonField[DeNum], size);

  my_fields.HM_density      = buf_HM.prepare(BaryonField[HMNum], size);
  my_fields.H2I_density     = buf_H2I.prepare(BaryonField[H2INum], size);
  my_fields.H2II_density    = buf_H2II.prepare(BaryonField[H2IINum], size);

  my_fields.DI_density      = buf_DI.prepare(BaryonField[DINum], size);
  my_fields.DII_density     = buf_DII.prepare(BaryonField[DIINum], size);
  my_fields.HDI_density     = buf_HDI.prepare(BaryonField[HDINum], size);

  my_fields.metal_density   = buf_metal.prepare(MetalPointer, size);

  my_fields.volumetric_heating_rate = buf_vol_heat.prepare(volumetric_heating_rate, size);
  my_fields.specific_heating_rate   = buf_spec_heat.prepare(specific_heating_rate, size);

#ifdef TRANSFER
  /* Find RT fields */
  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum,
        gammaNum, kphHMNum, kdissH2IINum;

  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum,
                                  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  /* unit conversion from Enzo RT units to CGS */
  float rtunits = erg_eV / TimeUnits;

  if( RadiativeTransfer ){
    my_fields.RT_HI_ionization_rate   = buf_rt_hi.prepare(BaryonField[kphHINum], size);

    if (RadiativeTransferHydrogenOnly == FALSE){
      my_fields.RT_HeI_ionization_rate  = buf_rt_hei.prepare(BaryonField[kphHeINum], size);
      my_fields.RT_HeII_ionization_rate = buf_rt_heii.prepare(BaryonField[kphHeIINum], size);
    }

    if (MultiSpecies > 1)
      my_fields.RT_H2_dissociation_rate = buf_rt_h2.prepare(BaryonField[kdissH2INum], size);

    /* need to convert to CGS units */
    for( i = 0; i < size; i++) BaryonField[gammaNum][i] *= rtunits;

    my_fields.RT_heating_rate = buf_rt_heat.prepare(BaryonField[gammaNum], size);
    

  }
#endif // TRANSFER

  /* Call the chemistry solver. */

  if (solve_chemistry(&grackle_units, &my_fields, (double) dt_cool) == FAIL){
    fprintf(stderr, "Error in Grackle solve_chemistry.\n");
    return FAIL;
  }

  buf_dens.copy_back();
  buf_te.copy_back();
  buf_v1.copy_back();
  buf_v2.copy_back();
  buf_v3.copy_back();
  buf_HI.copy_back();
  buf_HII.copy_back();
  buf_HeI.copy_back();
  buf_HeII.copy_back();
  buf_HeIII.copy_back();
  buf_e.copy_back();
  buf_HM.copy_back();
  buf_H2I.copy_back();
  buf_H2II.copy_back();
  buf_DI.copy_back();
  buf_DII.copy_back();
  buf_HDI.copy_back();
  buf_metal.copy_back();

  if (HydroMethod != Zeus_Hydro) {
    for (i = 0; i < size; i++) {
      BaryonField[TENum][i] = thermal_energy[i] +
        0.5 * POW(BaryonField[Vel1Num][i], 2.0);
      if(GridRank > 1)
        BaryonField[TENum][i] += 0.5 * POW(BaryonField[Vel2Num][i], 2.0);
      if(GridRank > 2)
        BaryonField[TENum][i] += 0.5 * POW(BaryonField[Vel3Num][i], 2.0);

      if( UseMHD ) {
        BaryonField[TENum][i] += 0.5 * (POW(BaryonField[iBx][i], 2.0) + 
                                        POW(BaryonField[iBy][i], 2.0) + 
                                        POW(BaryonField[iBz][i], 2.0)) / 
          BaryonField[DensNum][i];
      }

    } // for (int i = 0; i < size; i++)
  } // if (HydroMethod != Zeus_Hydro)

  if (temp_thermal == TRUE) {
    delete [] thermal_energy;
  }

#ifdef TRANSFER
  if (RadiativeTransfer){
    /* convert the RT units back to Enzo */
    for(i = 0; i < size; i ++) BaryonField[gammaNum][i] /= rtunits;

  }
#endif TRANSFER


  delete [] TotalMetals;
  delete [] g_grid_dimension;
  delete [] g_grid_start;
  delete [] g_grid_end;

  LCAPERF_STOP("grid_GrackleWrapper");

#endif
  return SUCCESS;
}
