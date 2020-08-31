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
#ifdef GRACKLE_MD
  int HeHIINum, DMNum   , HDIINum
    , CINum   , CIINum  , CONum     , CO2Num   , OINum   , OHNum
    , H2ONum  , O2Num   , SiINum    , SiOINum  , SiO2INum
    , CHNum   , CH2Num  , COIINum   , OIINum   , OHIINum , H2OIINum, H3OIINum, O2IINum
    , MgNum   , AlNum   , SNum      , FeNum
    , SiMNum  , FeMNum  , Mg2SiO4Num, MgSiO3Num, Fe3O4Num
    , ACNum   , SiO2DNum, MgONum    , FeSNum   , Al2O3Num
    , DustNum ;
#endif

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
#ifdef GRACKLE_MD 
    HeHIINum = DMNum     = HDIINum
  = CINum    =  CIINum   = CONum      = CO2Num    = OINum    = OHNum
  = H2ONum   =  O2Num    = SiINum     = SiOINum   = SiO2INum
  = CHNum    =  CH2Num   = COIINum    = OIINum    = OHIINum  = H2OIINum = H3OIINum = O2IINum
  = MgNum    =  AlNum    = SNum       = FeNum
  = SiMNum   = FeMNum    = Mg2SiO4Num = MgSiO3Num = Fe3O4Num
  = ACNum    =  SiO2DNum = MgONum     = FeSNum    = Al2O3Num
  = DustNum  = 0;
#endif
 
  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
		      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }
#ifdef GRACKLE_MD
    if (IdentifySpeciesFieldsMD( HeHIINum, DMNum   , HDIINum
                               , CINum   , CIINum  , CONum     , CO2Num   , OINum   , OHNum
                               , H2ONum  , O2Num   , SiINum    , SiOINum  , SiO2INum
                               , CHNum   , CH2Num  , COIINum   , OIINum   , OHIINum , H2OIINum,  H3OIINum,  O2IINum
                               , MgNum   , AlNum   , SNum      , FeNum
                               , SiMNum  , FeMNum  , Mg2SiO4Num, MgSiO3Num, Fe3O4Num
                               , ACNum   , SiO2DNum, MgONum    , FeSNum   , Al2O3Num
                               , DustNum ) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFieldsMD.\n");
    }
#endif
 
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
    if (debug)
      fprintf(stderr, "Warning: No metal field found.  Turning OFF MetalCooling.\n");
    MetalCooling = FALSE;
    MetalNum = 0;
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

  int ExtraType0Num, ExtraType1Num, ExtraType2Num, ExtraType3Num, ExtraType4Num, ExtraType5Num
    , ExtraType6Num, ExtraType7Num, ExtraType8Num, ExtraType9Num, ExtraType10Num,ExtraType11Num;
  if (MetalFieldPresent && MultiMetals)
    if (this->IdentifyExtraTypeFields(
      ExtraType0Num, ExtraType1Num, ExtraType2Num, ExtraType3Num, ExtraType4Num, ExtraType5Num,
      ExtraType6Num, ExtraType7Num, ExtraType8Num, ExtraType9Num, ExtraType10Num,ExtraType11Num
               ) == FAIL)
      ENZO_FAIL("Error in grid->IdentifyExtraTypeFields.\n");
 
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

  /* now add in the baryon fields */
  my_fields.density         = density;
  my_fields.internal_energy = thermal_energy;
  my_fields.x_velocity      = velocity1;
  my_fields.y_velocity      = velocity2;
  my_fields.z_velocity      = velocity3;

  if(MultiSpecies > 0) {
    my_fields.e_density       = BaryonField[DeNum];
    my_fields.HI_density      = BaryonField[HINum];
    my_fields.HII_density     = BaryonField[HIINum];
    my_fields.HeI_density     = BaryonField[HeINum];
    my_fields.HeII_density    = BaryonField[HeIINum];
    my_fields.HeIII_density   = BaryonField[HeIIINum];
  }

  if(MultiSpecies > 1) {
    my_fields.HM_density      = BaryonField[HMNum];
    my_fields.H2I_density     = BaryonField[H2INum];
    my_fields.H2II_density    = BaryonField[H2IINum];
  }

  if(MultiSpecies > 2) {
    my_fields.DI_density      = BaryonField[DINum];
    my_fields.DII_density     = BaryonField[DIINum];
    my_fields.HDI_density     = BaryonField[HDINum];
  }

  my_fields.metal_density   = MetalPointer;
  if(MultiMetals) {
    my_fields.metal_loc = BaryonField[ExtraType0Num];
    my_fields.metal_C13 = BaryonField[ExtraType1Num];
    my_fields.metal_C20 = BaryonField[ExtraType2Num];
    my_fields.metal_C25 = BaryonField[ExtraType3Num];
    my_fields.metal_C30 = BaryonField[ExtraType4Num];
    my_fields.metal_F13 = BaryonField[ExtraType5Num];
    my_fields.metal_F15 = BaryonField[ExtraType6Num];
    my_fields.metal_F50 = BaryonField[ExtraType7Num];
    my_fields.metal_F80 = BaryonField[ExtraType8Num];
    my_fields.metal_P170= BaryonField[ExtraType9Num];
    my_fields.metal_P200= BaryonField[ExtraType10Num];
    my_fields.metal_Y19 = BaryonField[ExtraType11Num];
  }

  if(UseDustDensityField)
    my_fields.dust_density = BaryonField[DustNum];

#ifdef GRACKLE_MD
  if(MultiSpecies > 3) {
    my_fields.     DM_density = BaryonField[     DMNum];
    my_fields.   HDII_density = BaryonField[   HDIINum];
    my_fields.  HeHII_density = BaryonField[  HeHIINum];
  }

  if(MetalChemistry > 0) {
    my_fields.     CI_density = BaryonField[     CINum];
    my_fields.    CII_density = BaryonField[    CIINum];
    my_fields.     CO_density = BaryonField[     CONum];
    my_fields.    CO2_density = BaryonField[    CO2Num];
    my_fields.     OI_density = BaryonField[     OINum];
    my_fields.     OH_density = BaryonField[     OHNum];
    my_fields.    H2O_density = BaryonField[    H2ONum];
    my_fields.     O2_density = BaryonField[     O2Num];
    my_fields.    SiI_density = BaryonField[    SiINum];
    my_fields.   SiOI_density = BaryonField[   SiOINum];
    my_fields.  SiO2I_density = BaryonField[  SiO2INum];
    my_fields.     CH_density = BaryonField[     CHNum];
    my_fields.    CH2_density = BaryonField[    CH2Num];
    my_fields.   COII_density = BaryonField[   COIINum];
    my_fields.    OII_density = BaryonField[    OIINum];
    my_fields.   OHII_density = BaryonField[   OHIINum];
    my_fields.  H2OII_density = BaryonField[  H2OIINum];
    my_fields.  H3OII_density = BaryonField[  H3OIINum];
    my_fields.   O2II_density = BaryonField[   O2IINum];
    if (GrainGrowth || DustSublimation) {
      if (DustSpecies > 0) {
        my_fields.   Mg_density = BaryonField[     MgNum];
      }
      if (DustSpecies > 1) {
        my_fields.   Al_density = BaryonField[     AlNum];
        my_fields.    S_density = BaryonField[      SNum];
        my_fields.   Fe_density = BaryonField[     FeNum];
      }
    }
  }

  if (GrainGrowth || DustSublimation) {
    if (DustSpecies > 0) {
      my_fields. MgSiO3_density = BaryonField[ MgSiO3Num];
      my_fields.     AC_density = BaryonField[     ACNum];
    }
    if (DustSpecies > 1) {
      my_fields.    SiM_density = BaryonField[    SiMNum];
      my_fields.    FeM_density = BaryonField[    FeMNum];
      my_fields.Mg2SiO4_density = BaryonField[Mg2SiO4Num];
      my_fields.  Fe3O4_density = BaryonField[  Fe3O4Num];
      my_fields.  SiO2D_density = BaryonField[  SiO2DNum];
      my_fields.    MgO_density = BaryonField[    MgONum];
      my_fields.    FeS_density = BaryonField[    FeSNum];
      my_fields.  Al2O3_density = BaryonField[  Al2O3Num];
    }
  }
#endif

  my_fields.volumetric_heating_rate = volumetric_heating_rate;
  my_fields.specific_heating_rate   = specific_heating_rate;

#ifdef TRANSFER
  /* Find RT fields */
  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum,
        gammaNum, kphHMNum, kdissH2IINum;

  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum,
                                  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);

  /* unit conversion from Enzo RT units to CGS */
  float rtunits = erg_eV / TimeUnits;

  if( RadiativeTransfer ){
    my_fields.RT_HI_ionization_rate   = BaryonField[kphHINum];

    if (RadiativeTransferHydrogenOnly == FALSE){
      my_fields.RT_HeI_ionization_rate  = BaryonField[kphHeINum];
      my_fields.RT_HeII_ionization_rate = BaryonField[kphHeIINum];
    }

    if (MultiSpecies > 1)
      my_fields.RT_H2_dissociation_rate = BaryonField[kdissH2INum];

    /* need to convert to CGS units */
    for( i = 0; i < size; i++) BaryonField[gammaNum][i] *= rtunits;

    my_fields.RT_heating_rate = BaryonField[gammaNum];
    

  }
#endif // TRANSFER

  int ISRFNum;

  if (grackle_data->use_isrf_field) {
    
    ISRFNum = FindField(ISRFHabing, FieldType, NumberOfBaryonFields);

    for( i = 0; i < size; i++) {
      BaryonField[ISRFNum][i] = BaryonField[kphHINum][i]
          / grackle_units.time_units            /* convert to CGS [1/s] */
          * (1.60217653e-12 * 13.6) / 6.30e-18  /* estimate flux */
          / 5.3e-3;                             /* convert to Habing units */
    }

    my_fields.isrf_habing = BaryonField[ISRFNum];

  }

  /* Call the chemistry solver. */

  if (solve_chemistry(&grackle_units, &my_fields, (double) dt_cool) == FAIL){
    fprintf(stderr, "Error in Grackle solve_chemistry.\n");
    return FAIL;
  }

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
#endif /* TRANSFER */


  delete [] TotalMetals;
  delete [] g_grid_dimension;
  delete [] g_grid_start;
  delete [] g_grid_end;

  LCAPERF_STOP("grid_GrackleWrapper");

#endif
  return SUCCESS;
}
