/******************************************************************
/
/ HANDLES FEEDBACK FOR INDIVIDUAL STAR PARTICLES
/
/ written by: Andrew Emerick
/ date:       August 2016
/ modified:
/
/ PURPOSE: Routine applies feedback from individual star particles
/          (stellar winds or supernovae) using the CIC stencil
/          interpolation methods from Simpson et. al. 2015. This
/          routine exists following the pre-existing Star Particle
/          class formalism in order to allow consistent feedback
/          across grids when particles are near grid boundaries.
/          This is a substantial improvement over previous method
/          which "shifted" feedback zone to avoid this.
/
/
/ OUTSTANDIG ISSUES: Type Ia supernovae rely on random number generator
/                    to go off... need to do this consistently (somehow)
/                    across all processors. Do in Set feedback flag
/
/
**********************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "performance.h"
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "phys_constants.h"
#include "Star.h"

#include "IndividualStarProperties.h"

#define DO_ZEROSOS 0

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);



int IsParticleFeedbackInGrid(float *pos, int ncell, LevelHierarchyEntry *Temp);
float ComputeOverlap(const int &i_shape, const float &radius,
                     const FLOAT &xc, const FLOAT &yc, const FLOAT &zc,
                     const FLOAT &xl, const FLOAT &yl, const FLOAT &zl,
                     const FLOAT &xr, const FLOAT &yr, const FLOAT &zr,
                     const int &nsample);

void ModifyStellarWindFeedback(float cell_mass, float T, float dx,
                               float MassUnits, float EnergyUnits, float &m_eject,
                               float &E_thermal, float * metal_mass,
                               float *grid_abundances);

float StellarYields_SolarAbundancesByNumber(const int &atomic_number);
float StellarYields_MMW(const int &atomic_number);
float StellarYields_ScaledSolarMassFractionByNumber(const float &metallicity,
                                                    const int &atomic_number);

int IndividualStarParticleAddFeedback(HierarchyEntry *Grids[],
                                      TopGridData *MetaData,
                                      LevelHierarchyEntry *LevelArray[],
                                      int level, Star* &AllStars,
                                      bool* &AddedFeedback){

  /* stars and hierarchy */
  Star *cstar;
  LevelHierarchyEntry *Temp;

  /* pos and vel for star */
  FLOAT *pos;
  float *vel;

  if (AllStars == NULL)
    return SUCCESS;

  TIMER_START("IndividualStarParticleAddFeedback");

  int count = 0;
  bool any_feedback_added = false;

  /* Loop over all stars, checking properties before doing feedback */
  for (cstar = AllStars; cstar; cstar = cstar->NextStar, count++){

    AddedFeedback[count] = false;

    if( ABS(cstar->ReturnType()) != IndividualStar &&
        ABS(cstar->ReturnType()) != IndividualStarWD &&
        ABS(cstar->ReturnType()) != IndividualStarRemnant &&
        ABS(cstar->ReturnType()) != IndividualStarPopIII  &&
        ABS(cstar->ReturnType()) != IndividualStarUnresolved ){
      continue; // This probably should never need to be checked
    }

    /* Check feedback flag - skip if particle isn't doing anything interesting */
    if( cstar->ReturnFeedbackFlag() < INDIVIDUAL_STAR_STELLAR_WIND ||
        cstar->ReturnFeedbackFlag() > INDIVIDUAL_STAR_POPIIISN){
      continue; // skip to next star
    }

    // This is fine IF the feedback is fully contained on that grid level
    //
    //    if( cstar->ReturnLevel() > level){
    //      continue; // only apply feedback on level of star
    //    }

    if(cstar->ReturnMass() < 0.0){
        cstar->PrintInfo();
        ENZO_FAIL("Particle Mass initially negative in IndividualStarParticleAddFeedback");
    }


    /* Get Units */
    float DensityUnits, LengthUnits, TemperatureUnits,
          TimeUnits, VelocityUnits, MassUnits, EnergyUnits;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                 &TimeUnits, &VelocityUnits, MetaData->Time) == FAIL){
        ENZO_FAIL("Error in GetUnits");
    }

    /* feedback is done in a cic interpolation. This is number of cells
       on eiher side of central cell (i.e. 3x3 CIC -> ncell = 1) */

    pos = cstar->ReturnPosition();
    vel = cstar->ReturnVelocity();

    double particle_mass;
    int gridnum=0;
    //
    // Check Stellar Winds - Apply if particle on grid and grid local
    //
    if( cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_STELLAR_WIND ||
        cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_WIND_AND_SN){
      for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){
        gridnum=0;
        for (Temp = LevelArray[l]; Temp; Temp = Temp ->NextGridThisLevel,gridnum++){

          int ncell = IndividualStarFeedbackStencilSize+1;
          if (IndividualStarFeedbackRadius > 0){
            ncell = (int) ceil(IndividualStarFeedbackRadius*pc_cm/LengthUnits / Temp->GridData->ReturnCellWidth());
          }

          if(Temp->GridData->isLocal() && IsParticleFeedbackInGrid(pos, ncell, Temp) ){
            // refresh mass every time to prevent double+ counting loss
            particle_mass = cstar->ReturnMass();

            Temp->GridData->IndividualStarAddFeedbackSphere(Temp->GridHierarchyEntry->NextGridNextLevel, cstar, &particle_mass, -1); // -1 = wind mode

            AddedFeedback[count] = TRUE;
          }
        }
      }

      if (AddedFeedback[count]){ // only if this particle did something
//        cstar->PrintInfo();
        float old_mass = cstar->ReturnMass();
        cstar->SetNewMass(particle_mass); // update mass (only once)
        cstar->AddToWindMassEjected(old_mass - particle_mass);
//        cstar->PrintInfo();
      }
    }

    //
    // Check Pop III Supernova
    if (cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_POPIIISN) {
      for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l++){
        gridnum=0;
        for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel, gridnum++){
          int ncell = IndividualStarFeedbackStencilSize+1;
          if (IndividualStarFeedbackRadius > 0){
            ncell = (int) ceil(IndividualStarFeedbackRadius*pc_cm/LengthUnits / Temp->GridData->ReturnCellWidth());
          }

          if(Temp->GridData->isLocal() && IsParticleFeedbackInGrid(pos, ncell, Temp)){

            particle_mass = cstar->ReturnMass();
            Temp->GridData->IndividualStarAddFeedbackSphere(Temp->GridHierarchyEntry->NextGridNextLevel, cstar, &particle_mass, 3); // 3 == popIII mode
            AddedFeedback[count] = TRUE;

          }
        }
      }

      if (AddedFeedback[count]){
        float old_mass = cstar->ReturnMass();
        cstar->SetFeedbackFlag(INDIVIDUAL_STAR_SN_COMPLETE);
        cstar->SetNewMass(particle_mass);
        cstar->AddToSNMassEjected(old_mass - particle_mass);
      }

    }
    //
    // Check Core Collapse Supernova
    if( cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_SNII ||
        cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_WIND_AND_SN){
      for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){
        gridnum=0;
        for (Temp = LevelArray[l]; Temp; Temp = Temp ->NextGridThisLevel,gridnum++){
          int ncell = IndividualStarFeedbackStencilSize+1;
          if (IndividualStarFeedbackRadius > 0){
            ncell = (int) ceil(IndividualStarFeedbackRadius*pc_cm/LengthUnits / Temp->GridData->ReturnCellWidth());
          }

          if(Temp->GridData->isLocal() && IsParticleFeedbackInGrid(pos, ncell, Temp)){
            // refresh mass every time to prevent double+ counting
            particle_mass = cstar->ReturnMass();
            Temp->GridData->IndividualStarAddFeedbackSphere(Temp->GridHierarchyEntry->NextGridNextLevel, cstar, &particle_mass, 1); // 1 == snII mode

            AddedFeedback[count] = TRUE;
          }
        }
      }

      if (AddedFeedback[count]){
        float old_mass = cstar->ReturnMass();
        AddedFeedback[count] = true;
        cstar->SetFeedbackFlag(INDIVIDUAL_STAR_SN_COMPLETE);
        cstar->SetNewMass(particle_mass); // update mass (only once)
        cstar->AddToSNMassEjected(old_mass - particle_mass);
      }
    }

    //
    // Check Type Ia Supernova
    //
    if( cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_SNIA){
      for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){
        gridnum = 0;
        for (Temp = LevelArray[l]; Temp; Temp = Temp ->NextGridThisLevel,gridnum++){
          int ncell = IndividualStarFeedbackStencilSize+1;
          if (IndividualStarFeedbackRadius > 0){
            ncell = (int) ceil(IndividualStarFeedbackRadius*pc_cm/LengthUnits / Temp->GridData->ReturnCellWidth());
          }

          if(Temp->GridData->isLocal() && IsParticleFeedbackInGrid(pos, ncell, Temp)){
            particle_mass = cstar->ReturnMass();

            Temp->GridData->IndividualStarAddFeedbackSphere(Temp->GridHierarchyEntry->NextGridNextLevel, cstar, &particle_mass, 2); // 2 == SNIa

          }
        }
      }

      AddedFeedback[count] = true;
      cstar->SetFeedbackFlag(INDIVIDUAL_STAR_SN_COMPLETE);
      cstar->AddToSNMassEjected(cstar->ReturnMass()); // not the actual mass ejection from SNIA !!!
      cstar->SetNewMass(0.0); // now a massless tracer
    }

    if(cstar->ReturnMass() < 0.0){

        if (IndividualStarIgnoreNegativeMass && cstar->ReturnFeedbackFlag() == INDIVIDUAL_STAR_SN_COMPLETE){
            /* This is experimental and obviously not physically valid. For use when using IMF averaged yields
               with individual particles which just scale mass ejection based on particle mass fraction which
               can lead to inconsistencies with mass ejection with winds + SN since they are IMF averaged
               independently without a consistency check. Fancier averaging would fix this.
               Setting remnant mass to 20% of birth mass, which is not an unreasonable estimate (maybe a bit high
               in some cases, but we're not doing a detailed dynamics sim here anyway). */
            cstar->SetNewMass(  0.2 * cstar->ReturnBirthMass());
        } else {
            cstar->PrintInfo();
            ENZO_FAIL("Particle Mass going negative in IndividualStarParticleAddFeedback");
        }
    }

    if (AddedFeedback[count]) any_feedback_added = true;

  } // end stars loop


  // debugging loop to ensure validity of mass ejection
  if (TRUE) {
    for (cstar = AllStars; cstar; cstar = cstar->NextStar){
      cstar->CheckMassEjectionValidity();
    }
  }

  /* Ensure injection is valid at all levels */
/*
//  if (any_feedback_added){
    for (int l = level; l > 0; l--){
      Temp = LevelArray[l];
      while (Temp != NULL) {
          if (Temp->GridData->ProjectSolutionToParentGrid(*Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL){
            fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid\n");
            return FAIL;
           }
            Temp = Temp->NextGridThisLevel;
      }
    }
  //}
*/

  TIMER_STOP("IndividualStarParticleAddFeedback");
  return SUCCESS;
}


int IsParticleFeedbackInGrid(float *pos, int ncell, LevelHierarchyEntry *Temp){
/* Check and see if particle feedback zone overlaps with any portion
   of grid before calling feedback functions. This is checked as well
   in feedback functions but reduces overhead somewhat by having
   a redundant check here
*/

  int Rank, Dims[MAX_DIMENSION];
  float CellWidth;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  Temp->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
  CellWidth = (RightEdge[0] - LeftEdge[0]) / (1.0 * Dims[0]);

  float fudge = 1.0;

  if( (pos[0] - (ncell + fudge)*CellWidth > RightEdge[0]) ||
      (pos[0] + (ncell + fudge)*CellWidth < LeftEdge[0])  ||
      (pos[1] - (ncell + fudge)*CellWidth > RightEdge[1]) ||
      (pos[1] + (ncell + fudge)*CellWidth < LeftEdge[1])  ||
      (pos[2] - (ncell + fudge)*CellWidth > RightEdge[2]) ||
      (pos[2] + (ncell + fudge)*CellWidth < LeftEdge[2])){
      // particle feedback zone is not on grid at all. skip
      // this check is performed also in actual feedback routines, but
      // redundancy is O.K. here
      return FALSE;
  }

  return TRUE;
}

int grid::IndividualStarAddFeedbackSphere(HierarchyEntry* SubgridPointer,
                                          Star *cstar, float *mp, const int mode){

/*
     General function to add feedback for a given star in a spherical region

     Thermal energy injection ONLY

     mode   :   integer, (Values)
                switches between stellar wind, core collapse SN, or type ia sn

*/
  if (this->NumberOfBaryonFields == 0 || !this->isLocal() )
    return SUCCESS;


  /* AJE Try something here */
  /* First, set under_subgrid field */
  HierarchyEntry *Subgrid;
  if (FALSE){
  this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  for (Subgrid = SubgridPointer; Subgrid; Subgrid = Subgrid->NextGridThisLevel){
    this->ZeroSolutionUnderSubgrid(Subgrid->GridData, ZERO_UNDER_SUBGRID_FIELD);
  }
  }
  

  float dx = this->CellWidth[0][0];

  float m_eject, E_thermal;

  const float mproj = cstar->ReturnBirthMass();
  const float lifetime = cstar->ReturnLifetime();


  float *metal_mass; // array of individual species masses

  FLOAT * pos;
  pos = cstar->ReturnPosition();

  /* Get Units */
  float DensityUnits, LengthUnits, TemperatureUnits,
        TimeUnits, VelocityUnits, MassUnits, EnergyUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits;
  EnergyUnits = MassUnits * VelocityUnits * VelocityUnits;

  /* If we are following yields, initialize array to hold ejecta masses */
  if(IndividualStarFollowStellarYields && MultiMetals == 2){

    metal_mass = new float[StellarYieldsNumberOfSpecies + 1];

    for (int i = 0; i < StellarYieldsNumberOfSpecies + 1; i ++){
      metal_mass[i] = 0.0;
    }

  } else { metal_mass = NULL;}


  float cgs_lifetime = lifetime * TimeUnits;

  int stellar_wind_mode = FALSE;

  if( mode < 0 ){  // compute properties for stellar wids

    // mproj needs to be in SolarMass - everything else in CGS
    IndividualStarSetStellarWindProperties(cstar, this->Time, this->dtFixed, TimeUnits,
                                           m_eject, E_thermal, metal_mass);
    stellar_wind_mode = TRUE;
  } else if (mode == 1){

    // core collapse supernova
    IndividualStarSetCoreCollapseSupernovaProperties(cstar, m_eject, E_thermal, metal_mass);

    stellar_wind_mode = FALSE;
  } else if (mode == 2){

    // Type Ia supernova properties
    IndividualStarSetTypeIaSupernovaProperties(m_eject, E_thermal, metal_mass);
    // printf("m_eject  for snia = %"FSYM"\n", m_eject);
    stellar_wind_mode = FALSE;

  } else if (mode == 3){

    IndividualStarSetPopIIISupernovaProperties(cstar, m_eject, E_thermal, metal_mass);

    stellar_wind_mode = FALSE;
  }

  /* Return surface abundances of stars */
  if (IndividualStarSurfaceAbundances){
    // This works under the assumption that yield tables provide just
    // the amount of each element produced (not production - ambient). By summing
    // here, we take the amount released for a given element (X) as:
    // ejected_mass_X = surface_abundance_X * total_ejecta_mass + yield_table_production_for_X
    //
    // --- likely this is never really a dominant effect to account for
    //
    double *abundances = cstar->ReturnAbundances();

    // this should vary with yield table assumptions of Z_solar and solar abund !!!
    const double z_solar = 0.0134; // Asplund+2009
    const double z_ratio = cstar->ReturnMetallicity() / z_solar;
    float dm_total = 0.0;

    const double Fe_H = StellarYields_SolarAbundancesByNumber(26);

    for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){

      if ((StellarYieldsAtomicNumbers[i] <= 2)) continue; // no H or He adjustments

      double a_solar;

      if (LimongiAbundances && Fe_H <= -1.0){
        double enhancement = 0.0;

        switch (StellarYieldsAtomicNumbers[i]){
          /* At Fe/H < solar, these abundances are enhanced by the following
             [X/Fe] values */
          case  6: enhancement = 0.18; break;
          case  8: enhancement = 0.47; break;
          case 12: enhancement = 0.27; break; // paper has 0.0.27 ...
          case 14: enhancement = 0.37; break;
          case 16: enhancement = 0.35; break;
          case 18: enhancement = 0.35; break;
          case 20: enhancement = 0.33; break;
          case 23: enhancement = 0.23; break;

          default:
            enhancement = 0.0;
        }

        a_solar = z_ratio * POW(10.0, enhancement + Fe_H) * 0.7381 * (StellarYields_MMW(StellarYieldsAtomicNumbers[i]) /
                                                          StellarYields_MMW(1));

      } else{
        a_solar = StellarYields_ScaledSolarMassFractionByNumber(cstar->ReturnMetallicity(),
                                                                       StellarYieldsAtomicNumbers[i]);
      }
      // abundances in stars are really mass fractions, so a_solar should be the scaled solar mass
      // fraction of the element (which is what is calculated above)
      double mass_change        = (cstar->abundances[i] - a_solar)*m_eject; // *z_ratio)*m_eject;
      if (TRUE) { // debugging
        if (ABS(mass_change) >= metal_mass[i+1]){
          printf("WARNING - SURFACE ABUNDANCES PRODUCING BIZZARE RESULTS\N");
          printf("Total mass change = %"ESYM" from initial = %"ESYM" and total ejection of %"ESYM"\n", mass_change, metal_mass[i+1], m_eject);
          printf("For element %"ISYM" with abundance %"ESYM" and scaled solar abundance %"ESYM"\n",StellarYieldsAtomicNumbers[i],cstar->abundances[i],a_solar);
          printf("NEED TO DOUBLE CHECK THESE VALUES");
        }
      }
      metal_mass[i+1] += mass_change;
      dm_total += mass_change;
    }

    /* strictly speaking this should be conserved if we follow ALL species
       but we don't. Unsure if I should leave as += 0 or += dm */
    metal_mass[0]     += dm_total ; // (cstar->ReturnMetallicity() - z_solar*z_ratio)*m_eject;
    m_eject += dm_total;

    if(TRUE){
      for (int i = 0; i < StellarYieldsNumberOfSpecies+1; i++){
        if ((metal_mass[i] <= 0) || (metal_mass[i] > metal_mass[0])){
            /* if any fail here, print all and exit */
          printf("Total mass = %"FSYM"\n", m_eject);
          for(int j=0;j<StellarYieldsNumberOfSpecies+1;j++){
            printf(" %"FSYM,metal_mass[j]);
          }
          printf("\n");
          ENZO_FAIL("Negative mass or too much mass in surface abundance return");
        }
      }
    }

  }

  /* convert computed parameters to code units */
  m_eject   = m_eject*SolarMass / MassUnits   / (dx*dx*dx);
  E_thermal = E_thermal      / EnergyUnits / (dx*dx*dx);

  if(IndividualStarFollowStellarYields && MultiMetals == 2){
    for(int i = 0; i < StellarYieldsNumberOfSpecies + 1; i++){
      // printf("metal mass species %"ISYM"   = %"ESYM"\n", i, metal_mass[i]);
      metal_mass[i] = metal_mass[i] * SolarMass / MassUnits / (dx*dx*dx);
    }
  }

  //
  // now that we've computed the explosion properties
  // find where we should go off
  //
  if( (m_eject > 0) || (E_thermal > 0)){ // can sometimes both be zero for stellar winds due to mass corrections
    this->IndividualStarInjectSphericalFeedback(cstar, pos[0], pos[1], pos[2], m_eject, E_thermal,
                                                metal_mass, stellar_wind_mode);
  }

  float new_mass = (*mp) - m_eject * (dx*dx*dx) * MassUnits / SolarMass; // update mass

  if( *mp < 0 && mode != 2){ // This can happen for Type 1a since using fixed mass model
    printf("new_mass = %"ESYM" mp = %"ESYM" m_eject =%"ESYM"\n", new_mass, *mp, m_eject*dx*dx*dx*MassUnits/SolarMass);
    ENZO_FAIL("IndividualStarFeedback: Ejected mass greater than current particle mass - negative particle mass!!!\n");
  } else if (mode == 2){
    *mp = 0.0;
  } else{
    *mp = new_mass;
  }

  delete [] metal_mass;

  /* Maybe? */
//  HierarchyEntry *Subgrid;
  if (FALSE){
  this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  for (Subgrid = SubgridPointer; Subgrid; Subgrid = Subgrid->NextGridThisLevel){
    this->ZeroSolutionUnderSubgrid(Subgrid->GridData, ZERO_UNDER_SUBGRID_FIELD);
  }
  }

  return SUCCESS;
}

int grid::IndividualStarInjectSphericalFeedback(Star *cstar,
                                                const FLOAT xp, const FLOAT yp, const FLOAT zp,
                                                const float m_eject, const float E_thermal,
                                                const float *metal_mass, const int stellar_wind_mode){

  float dx = float(this->CellWidth[0][0]); // for convenience

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;

  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  if ( MultiSpecies ) {
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  int PopIIIMetalNum, PopIIIPISNeMetalNum,
      AGBMetalNum, SNIaMetalNum, SNIIMetalNum, RProcMetalNum,
      ExtraMetalNum0, ExtraMetalNum1, ExtraMetalNum2;

  AGBMetalNum    = FindField(ExtraType0, FieldType, NumberOfBaryonFields);
  PopIIIMetalNum = FindField(ExtraType1, FieldType, NumberOfBaryonFields);
  PopIIIPISNeMetalNum = FindField(MetalPISNeDensity, FieldType, NumberOfBaryonFields);
  SNIaMetalNum   = FindField(MetalSNIaDensity, FieldType, NumberOfBaryonFields);
  SNIIMetalNum   = FindField(MetalSNIIDensity, FieldType, NumberOfBaryonFields);
  RProcMetalNum  = FindField(MetalRProcessDensity, FieldType, NumberOfBaryonFields);

  ExtraMetalNum0 = FindField(ExtraMetalField0, FieldType, NumberOfBaryonFields);
  ExtraMetalNum1 = FindField(ExtraMetalField1, FieldType, NumberOfBaryonFields);
  ExtraMetalNum2 = FindField(ExtraMetalField2, FieldType, NumberOfBaryonFields);

  if ( IndividualStarTrackAGBMetalDensity && (AGBMetalNum <= 0)){
    ENZO_FAIL("Error in finding AGB metal density field in individual_star_maker");
  }

  if ( IndividualStarPopIIIFormation && ((PopIIIMetalNum <= 0) || (PopIIIPISNeMetalNum <=0))){
    ENZO_FAIL("Error in finding Pop III metal density field in individual_star_maker");
  }

  if ( IndividualStarTrackSNMetalDensity && ( (SNIIMetalNum <= 0) || (SNIaMetalNum <=0))){
    ENZO_FAIL("Error in finding SNII and SNIa metal density field in individual_star_maker.");
  }

  if ( IndividualStarRProcessModel && (RProcMetalNum <= 0)) {
    ENZO_FAIL("Error in finding R Process model metal density field in individual_star_maker.");
  }

  if (IndividualStarSNIaModel == 2){
    if ( (ExtraMetalNum0 <= 0) || (ExtraMetalNum1 <= 0) || (ExtraMetalNum2 <= 0)){
      ENZO_FAIL("Error in finding all SNIa fields needed for model 2.");
    }
  }

  /* Get Units */
  float DensityUnits, LengthUnits, TemperatureUnits,
        TimeUnits, VelocityUnits, MassUnits, EnergyUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits;
  EnergyUnits = MassUnits * VelocityUnits * VelocityUnits;

  FLOAT xstart = this->CellLeftEdge[0][0],
        ystart = this->CellLeftEdge[1][0],
        zstart = this->CellLeftEdge[2][0];

  const int nx = *(this->GridDimension), ny = *(this->GridDimension+1), nz = *(this->GridDimension+2);
  const int size = nx*ny*nz;

  //
  // find position and index for grid zone
  // nearest to particle
  //
  const float xpos = (xp - xstart)/dx;
  const float ypos = (yp - ystart)/dx;
  const float zpos = (zp - zstart)/dx;

  const int ic   = ((int) floor(xpos ));
  const int jc   = ((int) floor(ypos ));
  const int kc   = ((int) floor(zpos ));

  float * temperature;
  temperature = new float[size];


  if(  this->ComputeTemperatureField(temperature) == FAIL ){
    ENZO_FAIL("Error in compute temperature called from PhotoelectricHeatingFromStar");
  }

  //
  //
  // now move outward from particle computing
  // fractional volume to inject
  //
  float radius = 0.0;
  int   r_int  = 0;

  if (IndividualStarFeedbackRadius > 0){
    radius = IndividualStarFeedbackRadius * pc_cm / LengthUnits; // code length
    r_int  = (int) ceil(radius / dx); // int of farthest cell in any dir.
  } else {
    radius = IndividualStarFeedbackStencilSize * dx; // code length
    r_int  = ceil(IndividualStarFeedbackStencilSize); // int of farthest cell in any dir.
  }

  const float sphere_volume = 4.0 * pi * radius * radius * radius / 3.0; // (code length)**3
  const float cell_volume_fraction = dx*dx*dx / sphere_volume;           // fraction of vol for each cell

  float * injected_metal_mass;
  injected_metal_mass = new float [StellarYieldsNumberOfSpecies+1];
  for(int i =0; i < StellarYieldsNumberOfSpecies+1; i++) injected_metal_mass[i] = 0.0;

  if (metal_mass == NULL && (cstar)){
    injected_metal_mass[0] = cstar->ReturnMetallicity() * m_eject;
  } else if (metal_mass == NULL){
    injected_metal_mass[0] = 0.0;
  }

  // for printing stats at the end
  float total_volume_fraction = 0.0, total_grid_mass = 0.0, total_mass_injected = 0.0;
  float total_energy_injected = 0.0, max_density_on_grid = 0.0, average_density_on_grid = 0.0;
  float total_metal_mass = 0.0;
  int   cells_this_grid = 0;

  const int num_factors = POW(r_int*2+1,3);
  float * injection_factors;
  injection_factors = new float[num_factors];
  for(int i =0; i <num_factors; i++) injection_factors[i]=0.0;

  // loop over all cells, compute fractional volume to deposit feedback
  int count = 0;
  float total_injection_volume=0.0;
  for(int k = kc - r_int; k <= kc + r_int; k++){
    if ( (k<0) || (k>=nz) ) continue;
    for(int j = jc - r_int; j <= jc + r_int; j++){
      if ( (j<0) || (j>=ny) ) continue;
      for(int i = ic - r_int; i <= ic + r_int; i++){
        if ( (i<0) || (i>=nx) ) continue;

        int index = i + (j + k*ny)*nx;

        float xc,yc,zc;
        xc = (i + 0.5) * dx + xstart;
        yc = (j + 0.5) * dx + ystart;
        zc = (k + 0.5) * dx + zstart;

        // do fractional overlap calculation
        float fractional_overlap = ComputeOverlap(1, radius, xp, yp, zp,
                                            xc - dx, yc - dx, zc - dx,
                                            xc + dx, yc + dx, zc + dx, IndividualStarFeedbackOverlapSample);

        injection_factors[count] = cell_volume_fraction * fractional_overlap;
        total_injection_volume   += injection_factors[count];
        count++;
      }
    }
  } // end fractional overlap computation

  // injection_fraction represents fraction of total volume per cell
  // do some error checking:

  /* AJE: This really only works if particle feedback is contained entirely 
          within a single grid. Otherwise this correction will over / undercorrect
          when particle ejects feedback across multiple grids (unless
          some MPI occurs to share the correction factors... which... woudl be 
          dumb... could maybe do a check to see if its contained fully in grid
          and then apply this, but otherwise don't do it (and accept the
          (~< 1 % error) */
  if (FALSE){
  if (total_injection_volume < 1.0){
    // this will always underestimte the total volume
    // need to scale up all of the injection fractions
    // to ensure sum is 1 (and mass + energy conservation)
    const float total_inv = 1.0 / total_injection_volume;

    if (total_injection_volume < 0.8){

      printf("total_injection_volume = %"ESYM"\n",total_injection_volume);
      ENZO_FAIL("Total injection volume is small... this is a large correction... make sure nohing is broken");
    }

    for (count = 0; count < num_factors; count++){
      if (injection_factors[count] < 0) {ENZO_FAIL("injection factor < 0");}

      injection_factors[count] *= total_inv;
    }

  } else if (total_injection_volume > 1.02){
    // here for testing
    // but if this is common, just switch to always dividing by the sum (above)
    printf("total_injection_volume = %"FSYM"\n", total_injection_volume);
    ENZO_FAIL("Error in computing feedback injection volume in individual stars. Greater than 1");
  }
  }


  // loop over cells and inject feedback
  count = 0;
  for(int k = kc - r_int; k <= kc + r_int; k++){
    if ( (k<0) || (k>=nz) ) continue;
    for(int j = jc - r_int; j <= jc + r_int; j++){
      if ( (j<0) || (j>=ny) ) continue;
      for(int i = ic - r_int; i <= ic + r_int; i++){
        if ( (i<0) || (i>=nx) ) continue;

        int index = i + (j + k*ny)*nx;

        float injection_factor = injection_factors[count];

        if (injection_factor <= 0){
          count++;
          continue; // no overlap, no feedback
        }

        float delta_mass  = m_eject   * injection_factor;
        float delta_therm = E_thermal * injection_factor;


        if (IndividualStarFollowStellarYields && (cstar || metal_mass)){
          for(int im = 0; im < StellarYieldsNumberOfSpecies+1; im++){

            // Hack  - remove contribution from stars to experiment fields
            //          remember, im = 0 is the total metal mass field
            if ((cstar) && (im > 0) && MetalMixingExperiment){
              for (int j = 0; j < StellarYieldsNumberOfSpecies; j++){
                if (StellarYieldsAtomicNumbers[im-1] == MixingExperimentData.anums[j] ){
                  injection_factor = 0.0; // zero out these fields to limit to just mixing experiments events ONLY
                }
              }
            }

            injected_metal_mass[im] = metal_mass[im]*injection_factor;
          }
        }

        if(stellar_wind_mode){
          // Stellar winds are challenging - to avoid both very high velocities
          // on the grid and superheated gas, modify feedback accordingly
          ModifyStellarWindFeedback(BaryonField[DensNum][index],
                                    temperature[index], dx, MassUnits, EnergyUnits,
                                    delta_mass, delta_therm, injected_metal_mass,
                                    this->AveragedAbundances);
        }

        float inv_dens = 1.0 / (BaryonField[DensNum][index] + delta_mass);

        /* record statistics accumulators */
        cells_this_grid++;
        total_volume_fraction += injection_factor;
        total_mass_injected   += delta_mass;
        total_energy_injected += delta_therm;
        total_grid_mass       += BaryonField[DensNum][index];
        max_density_on_grid     = fmax( BaryonField[DensNum][index], max_density_on_grid);

        BaryonField[TENum][index] = (BaryonField[TENum][index] * BaryonField[DensNum][index]
                                     + delta_therm) * inv_dens;
        if(DualEnergyFormalism){
          BaryonField[GENum][index] = (BaryonField[GENum][index] * BaryonField[DensNum][index]
                                       + delta_therm) * inv_dens;
        }
        float old_mass = BaryonField[DensNum][index];
        BaryonField[DensNum][index] += delta_mass;

        /* add metal species if we need to */
        if(MultiMetals == 2 && IndividualStarFollowStellarYields){
          int field_num;
          this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0); // gives metallicity field


          if (cstar){ // is feedback coming from a star particle?
            BaryonField[field_num][index] += injected_metal_mass[0]; // add to metallicity field

            // Add to separate source fields if they exist

            if (IndividualStarPopIIIFormation && cstar->IsPopIII() ){
              BaryonField[PopIIIMetalNum][index] += injected_metal_mass[0];

              if ((cstar->ReturnBirthMass() > PISNLowerMass) && (cstar->ReturnBirthMass() < PISNUpperMass)){
                BaryonField[PopIIIPISNeMetalNum][index] += injected_metal_mass[0];
              }

            } else if (IndividualStarTrackAGBMetalDensity &&
                      (cstar->ReturnBirthMass() < IndividualStarSNIIMassCutoff) &&
                      (stellar_wind_mode)){                                 // make sure we are an AGB wind
              BaryonField[AGBMetalNum][index] += injected_metal_mass[0];

            } else if (IndividualStarTrackSNMetalDensity && !(stellar_wind_mode)){

              if (IndividualStarRProcessModel){

                if ( (cstar->ReturnBirthMass() >= IndividualStarRProcessMinMass) &&
                     (cstar->ReturnBirthMass() <= IndividualStarRProcessMaxMass) ){

                  const float RProcModelMass = (1.0 * SolarMass) / MassUnits; // fully arbitrary

                  BaryonField[RProcMetalNum][index] += RProcModelMass * injection_factors[count]; //need to do correct value here

                }
              }

              if (cstar->ReturnBirthMass() > IndividualStarSNIIMassCutoff){
                BaryonField[SNIIMetalNum][index] += injected_metal_mass[0];

              } else if (cstar->ReturnBirthMass() < IndividualStarSNIaMaximumMass){

                if (IndividualStarSNIaModel == 1){

                  BaryonField[SNIaMetalNum][index] += injected_metal_mass[0];

                } else if (IndividualStarSNIaModel == 2){

                  // 0 : DDS, 1: sCh, 2: SDS, 3 :HeRS
                  const int fields[4] = {SNIaMetalNum, ExtraMetalNum0, ExtraMetalNum1, ExtraMetalNum2};

                  BaryonField[fields[cstar->ReturnSNIaType()]][index] += injected_metal_mass[0];

                }

              }

            } // end if tracking yield source modes

          } else if (metal_mass){
            BaryonField[field_num][index] += injected_metal_mass[0];

          } else {
            // keep same fraction if using artificial SN generaotr
            BaryonField[field_num][index] += delta_mass *
                                             BaryonField[field_num][index] / old_mass;
          }

            total_metal_mass += BaryonField[field_num][index];

          for(int im = 0; im < StellarYieldsNumberOfSpecies; im++){
            this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num,
                                                              StellarYieldsAtomicNumbers[im]);
            if (cstar || metal_mass){
              BaryonField[field_num][index] += injected_metal_mass[1 + im];
            } else { // keep same fraction if using artificial SN generator
              BaryonField[field_num][index] += delta_mass *
                                               BaryonField[field_num][index] / old_mass;
            }
          }

        } else{
          int field_num;
          this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0); // gives metallicity field

          if (cstar || metal_mass){
            BaryonField[field_num][index] += injected_metal_mass[0];
          } else{
            BaryonField[field_num][index] += delta_mass *
                                             BaryonField[field_num][index]/ old_mass;
          }
        } // end yields check

        count++; //increment injection factor index
      }
    }
  }


  // print SN stats to check if resolved if desired
  if (IndividualStarPrintSNStats && (!stellar_wind_mode) && (cstar || metal_mass)){
    // Column order: Grid ID, Particle ID, M_now, M_eject, Sphere Volume

    average_density_on_grid = total_grid_mass / (1.0 * cells_this_grid); // Sum Density / # cells

    /* convert to CGS */
    total_grid_mass         *= dx*dx*dx * MassUnits;
    total_mass_injected     *= dx*dx*dx * MassUnits;
    total_energy_injected   *= dx*dx*dx * EnergyUnits;
    total_metal_mass        *= dx*dx*dx * MassUnits;
    const float sphere_volume_cgs   = sphere_volume / (LengthUnits * LengthUnits * LengthUnits);
    max_density_on_grid     *= DensityUnits;
    average_density_on_grid *= DensityUnits;
    const float m_eject_cgs  = m_eject * dx*dx*dx * MassUnits;
    const float average_metallicity      = total_metal_mass / (total_mass_injected + total_grid_mass);

    /* compute Sedov-Taylor phase radius (R_dps) */

    if (cstar){
      printf("IndividualStarSNStats: %"ISYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",
            this->ID, cstar->ReturnID(), this->Time, cstar->ReturnMass(), cstar->ReturnBirthMass(), cstar->ReturnMetallicity(), m_eject_cgs,
            cells_this_grid, sphere_volume_cgs, total_volume_fraction, total_mass_injected, total_energy_injected,
            total_grid_mass, max_density_on_grid, average_density_on_grid, total_metal_mass, average_metallicity);
    } else {
      printf("IndividualStarSNStats: %"ISYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",
            this->ID, -1.0, this->Time, -1.0, -1.0, -1.0, m_eject_cgs,
            cells_this_grid, sphere_volume_cgs, total_volume_fraction, total_mass_injected, total_energy_injected,
            total_grid_mass, max_density_on_grid, average_density_on_grid, total_metal_mass, average_metallicity);

    }
  }


  delete [] temperature;
  delete [] injection_factors;
  delete [] injected_metal_mass;
  // done with spherical injection feedback

  return SUCCESS;
}

void ModifyStellarWindFeedback(float cell_mass, float T, float dx,
                               float MassUnits, float EnergyUnits, float &m_eject,
                               float &E_thermal, float * metal_mass,
                               float *grid_abundances){

/* ----------------------------------------------------------------------------
 * ModifyStellarWindFeedback
 * -----------------------------------------------------------------------------
 * A. Emerick - Sep 2016
 *
 * Stellar winds are very challenging to get right, even at very high resolution.
 * We assume feedback here is used on ~1-5 pc resolution, galaxy scale, Gyr sims,
 * so feedback needs to be modified to be tractable. Inject thermal energy at
 * 100% thermalization of wind KE first (E_thermal_max) using velocities from
 * full wind model. However, if cell temperature gets very hot (~10^6 K, as set
 * by IndividualStarWindTemperature) then employ maximum wind velocity cutoff,
 * changing energy injection to E_thermal_min. If E_thermal_min will large large
 * temperatures above IndividualStarWindTemperature (as happens if soruce region
 * is devoid of gas), we use an ISM mass loading model to account for shell
 * mixing, which will be heineously unresolved at 1 pc resolution (need ~0.01 pc).
 * ----------------------------------------------------------------------------- */

  if (!IndividualStarUseWindMixingModel){
     return;
  }

  const float est_mu  = 0.5; // estimated - this is approximate anyway

  float m_ism = 0.0;

  E_thermal = E_thermal * dx *dx *dx * EnergyUnits;
  m_eject   = m_eject * dx *dx *dx * MassUnits;
  cell_mass = cell_mass *dx*dx*dx*MassUnits;

  float T_final = (E_thermal + 1.5 * cell_mass * kboltz * T / (est_mu*mh)) *
                  (2.0 * est_mu * mh/(3.0 * kboltz * (cell_mass + m_eject)));

    if(T_final > IndividualStarWindTemperature || T > IndividualStarWindTemperature){
      /* Compute the mass that needs to be injected */
      float E_final = (3.0/2.0) * (cell_mass/(est_mu*mh))*kboltz * T + E_thermal;
      T_final = fmax(T, IndividualStarWindTemperature);
      m_ism   = fmax( (E_final * (2.0 * est_mu * mh)/(3.0*kboltz * T_final)) - cell_mass - m_eject, 0.0);

      /* modify metal abundances here */
      m_ism = m_ism / (dx*dx*dx) / MassUnits;
      if(MultiMetals == 2 && IndividualStarFollowStellarYields && m_ism > 0.0){
        for (int im = 0; im < StellarYieldsNumberOfSpecies + 1; im++){

          metal_mass[im] = metal_mass[im] + m_ism * grid_abundances[im]; // ISM abundances

          if(metal_mass[im] < 0.0){
            printf("metal_mass %"ESYM" %"ISYM" %"ESYM"\n", metal_mass[im], im, grid_abundances[im]);
            ENZO_FAIL("IndividualStarFeedback: Metal mass correction < 0 and m_ism >0");
          }
        }
      } else{
        for (int im = 0; im < StellarYieldsNumberOfSpecies + 1; im++){
          if(metal_mass[im] < 0.0){
            printf("metal_mass %"ESYM" %"ISYM"\n", metal_mass[im], im);
            ENZO_FAIL("IndividualStarFeedback: Metal mass correction < 0 and m_ism < 0");
          }
        }
      }

  }

  /* make sure things aren't whacky */
  if (E_thermal < 0.0 || m_eject < 0.0 || m_ism < 0.0){
    printf("Error in stellar wind calculation. E_thermal = %"ESYM" m_eject = %"ESYM" m_ism (code) = %"ESYM"\n",E_thermal, m_eject, m_ism);
    ENZO_FAIL("IndividualStarFeedback: Negative injection values in stellar wind feedback modification\n");
  }

  /* convert back into code units */
  E_thermal  = E_thermal / (dx*dx*dx) / EnergyUnits;
  m_eject    = (m_eject)  / (dx*dx*dx) / MassUnits + m_ism;

}

float ComputeOverlap(const int &i_shape, const float &radius,
                     const FLOAT &xc, const FLOAT &yc, const FLOAT &zc,
                     const FLOAT &xl, const FLOAT &yl, const FLOAT &zl,
                     const FLOAT &xr, const FLOAT &yr, const FLOAT &zr,
                     const int &nsample){
 /* -------------------------------------------------------------------------
  * ComputeVolumeOverlap
  * -------------------------------------------------------------------------
  * A. Emerick - 9/21/16
  *
  * Computes overlap between a given rectangular grid cell and a spherical or
  * cylindrical volume using a Monte Carlo sampling method.
  *
  * Adopted from Joshua Wall's version of David Clarke's method in ZEUS-MP
  *
  * INPUTS:
  *   i_shape    - switch for shape (1 == sphere, 2 == right cylinder)
  *   xc,yc,zc   - center coordinates of sphere
  *   xl,yl,zl   - coordinates of lower left corner of grid cell
  *   xr,yr,zr   - coordinates of upper right corner of grid cell
  *   nsample    - number of Monte Carlo sample points
  *
  * OUTPUTS:
  *   overal     - fractional volume overlap for given cell
  * ------------------------------------------------------------------------- */

  float xsq[nsample], ysq[nsample], zsq[nsample];


  float factor = 1.0;
  if (i_shape == 2) factor = 0.0; // cylinder

  //
  float dx = (xr - xl) / ((float) nsample);
  float dy = (yr - yl) / ((float) nsample);
  float dz = (zr - zl) / ((float) nsample);

  for( int i = 0; i < nsample; i++){
    xsq[i] = (xl + (0.5 + i)*dx - xc)*(xl + (0.5+i)*dx - xc);
    ysq[i] = (yl + (0.5 + i)*dy - yc)*(yl + (0.5+i)*dy - yc);
    zsq[i] = (zl + (0.5 + i)*dz - zc)*(zl + (0.5+i)*dz - zc);
  }

  int inside_count = 0;

  for (int k = 0; k < nsample; k++){
    for (int j = 0; j < nsample; j++){
      for (int i = 0; i < nsample; i++){
       float r = sqrt( factor * xsq[i] + ysq[j] + zsq[k] );

       if ( r <= radius) inside_count++;
      }
    }
  }


  return ((float) inside_count) / ((float) nsample*nsample*nsample);
}
