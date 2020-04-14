
/****************************************************************************
/
/ STAR FORMATION AND FEEDBACK ALGORITHMS FOR INDIVIDUAL STARS
/
/ written by: Andrew Emerick
/ date:       February, 2016
/ modified1:
/
/ Controls star formation for individual stars as sampled from an IMF. Stars
/ are formed stochastically following an adaptation of Goldbaum et. al. 2015
/ as in star_maker_ssn.F
/ First use case of these particles is to tie to galaxy scale chemodynamics.
/ Particle creation tags these stars with local chemical abundances and feedback
/ is tied to yield tables to deposit elements in the ISM.
*****************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
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
#include "StellarYieldsRoutines.h"


/* Following Grid_ComputeTemperatureField.C */
#ifndef MU_METAL
# define MU_METAL 16.0
#endif

/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
int FindField(int f, int farray[], int n);
unsigned_long_int mt_random(void);

int search_lower_bound(float *arr, float value, int low, int high, int total);

extern "C" void FORTRAN_NAME(pop3_properties)(FLOAT *mass, FLOAT* luminosity,
                                              FLOAT *lifetime);


/* Internal */
float SampleIMF(void);
float SampleIMF(float * data, const float & lower_mass, const float & upper_mass);
float SamplePopIII_IMF(void);






int grid::individual_star_maker(float *dm, float *temp, int *nmax, float *mu, int *np,
                                float *ParticleMass,
                                int *ParticleType, FLOAT *ParticlePosition[],
                                float *ParticleVelocity[], float *ParticleAttribute[],
                                float *StellarAbundances[]){
/*-----------------------------------------------------------------------------
  INPUTS:
    dm          - dark matter density field (computed in Grid_StarParticleHandler)
    temp        - temperature field (computed in Grid_StarParticleHandler)
    nmax        - Maximum allowed number of stars that can form on a single grid
    mu          - global Mean Molecular weight of gas
    ctype       - number for desired particle type assignment

  OUTPUTS: SUCCESS or FAIL
    Creates star particle and updates all particle arrays
    modifies baryon fields during star formation
    np - number of particles created
    ParticleMass - particle masses on grid
    ParticleType - particle types on grid
    ParticlePosition - particle positions on grid
    ParticleVelocity - particle velocities
    ParticleAttribute - particle attributes
-----------------------------------------------------------------------------*/

  const double sndspdC = 1.3095e8;

  int i, j, k, index, ii=0, istar=0, index_presf=0;
  int xo, yo, zo, rsign=1;
  float bmass, div, star_mass=0.0, sum_mass=0.0, metal_mass=0.0, H2mass=0.0;
  float pstar, mass_to_stars, mass_available, tdyn;
  float dtot, isosndsp2, jeansmass, star_fraction;
  float umean, vmean, wmean, px, py, pz, px_excess, py_excess, pz_excess;
  float rnum;
  float inv_metal_mol = 1.0 /MU_METAL;

  const int max_random = (1<<16);

  int form_method = -1; // tracker for debugging purposes

  /* for convenience, rename some grid properties - will likely get optimized out */
  int  nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  FLOAT xstart = CellLeftEdge[0][0], ystart = CellLeftEdge[1][0], zstart = CellLeftEdge[2][0];
  float   dx   = CellWidth[0][0];

  if (! this->isLocal()) return SUCCESS;

  if ( this->dtFixed == 0.0){
    printf("DT EQUAL TO ZERO\N");
    return FAIL;
  }

  /* obtain baryon field indexes */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* identify species fields if they exist for proper computation of Mu */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if ( MultiSpecies ){
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  /* get metallicity tracer field number */
  int MetalNum;
  MetalNum   = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);

 int PopIIIMetalNum, PopIIIPISNeMetalNum, AGBMetalNum,
     SNIaMetalNum, SNIIMetalNum, RProcMetalNum;

  AGBMetalNum    = FindField(ExtraType0, FieldType, NumberOfBaryonFields);
  PopIIIMetalNum = FindField(ExtraType1, FieldType, NumberOfBaryonFields);
  PopIIIPISNeMetalNum = FindField(MetalPISNeDensity, FieldType, NumberOfBaryonFields);
  SNIaMetalNum   = FindField(MetalSNIaDensity, FieldType, NumberOfBaryonFields);
  SNIIMetalNum   = FindField(MetalSNIIDensity, FieldType, NumberOfBaryonFields);
  RProcMetalNum  = FindField(MetalRProcessDensity, FieldType, NumberOfBaryonFields);

  if ( IndividualStarTrackAGBMetalDensity && (AGBMetalNum <= 0)){
    ENZO_FAIL("Error in finding AGB metal density field in individual_star_maker");
  }

  if ( IndividualStarPopIIIFormation && ((PopIIIMetalNum <= 0) || (PopIIIPISNeMetalNum <= 0))){
    ENZO_FAIL("Error in finding Pop III metal density field in individual_star_maker");
  }

  if ( IndividualStarTrackSNMetalDensity && ( (SNIIMetalNum <= 0) || (SNIaMetalNum <=0))){
    ENZO_FAIL("Error in finding SNII and SNIa metal density field in individual_star_maker.");
  }

  if ( IndividualStarRProcessModel && ( RProcMetalNum <=0 )){
    ENZO_FAIL("Error in finding R process metal density field in individual_star_maker.");
  }

  /* get units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units


    // 3D -> 1D index
    xo = 1;
    yo = nx;
    zo = (nx) * (ny);

    // over density threshold in code units
    // if multispecies is off, assumes a value for MU
    float odthreshold, secondary_odthreshold;
    if (MultiSpecies == FALSE){
        odthreshold           = StarMakerOverDensityThreshold * mh * (*mu) / (DensityUnits); // code density
        secondary_odthreshold = IndividualStarSecondaryOverDensityThreshold * mh * (*mu) / (DensityUnits);
    }

    // loop over all cells, check condition, form stars stochastically
    ii = 0; index_presf = 0;

    int number_of_sf_cells = 0;
    int integer_sep = ((int) (IndividualStarCreationStencilSize + 1) / 2.0 - 1); // stencil size must be odd number

    float *ke_before;
    if (HydroMethod != 2){
      ke_before = new float[ ((int) POW(IndividualStarCreationStencilSize,3)) ];

      for (i = 0; i < ((int) POW(IndividualStarCreationStencilSize,3)); i++) ke_before[i] = 0.0;

    } else{
      ke_before = NULL;
    }

    for (k = ibuff; k < nz - ibuff; k++){
      for (j = ibuff; j < ny - ibuff; j++){
        for (i = ibuff; i < nx - ibuff; i++){

          index = i + ( j + k * (ny)) * nx;

          if (ProblemType == 30){
            if ( BaryonField[NumberOfBaryonFields][index] != 0) continue; // not highest refined zone of this region
          }

          /* if distributed star formation */
          // check center cell's SF condition, if met, do SF
          /* loop and sum over all*/

           bmass      = (BaryonField[DensNum][index]*(dx*dx*dx)) * MassUnits / SolarMass; // in solar masses
           metal_mass = (BaryonField[MetalNum][index] * bmass); // metalNum is converted to fraction in Grid_StarParticleHandler

               // perform the following easy checks for SF before proceeding
               // 1) Is density greater than the density threshold?
               // 2) Is temperature < the minimum temperature?
               // 3) Do not allow star formation if minimum star particle mass is
               //    above some fraction of cell mass. This is very unlikely to occur
               //    in intended use case:
               //    (i.e. baryon mass likely always going to be > ~10 solar masses)

          sum_mass = 0.0; index_presf = ii;

          /* Need to compute Mu exactly and recompute threshold */
          if ( MultiSpecies ){
            float mu_cell, number_density;

            number_density =
              0.25*(BaryonField[HeINum][index] + BaryonField[HeIINum][index] +
                    BaryonField[HeIIINum][index])                        +
                   BaryonField[HINum][index] + BaryonField[HIINum][index]    +
                   0.5 * BaryonField[DeNum][index];
            /* if H2 is present */
            if (MultiSpecies > 1){
              number_density += BaryonField[HMNum][index] +
                          0.5*(BaryonField[H2INum][index] + BaryonField[H2IINum][index]);

              H2mass = (BaryonField[H2INum][index] + BaryonField[H2IINum][index]) * bmass;
            }

            /* Metal field must be present in this star formation scheme */
            number_density += BaryonField[MetalNum][index] * inv_metal_mol;


            number_density *= BaryonField[DensNum][index] / mh ; // now actual n density

            mu_cell = BaryonField[DensNum][index] / (number_density * mh);
            odthreshold = StarMakerOverDensityThreshold * (mu_cell) * mh / (DensityUnits);

            secondary_odthreshold = IndividualStarSecondaryOverDensityThreshold * (mu_cell) * mh / DensityUnits;
          }


          if (   BaryonField[DensNum][index]      > odthreshold
              && temp[index] <= IndividualStarTemperatureThreshold){
              //&& IndividualStarMassFraction*bmass > IndividualStarIMFLowerMassCutoff
              //&& 0.5*bmass > IndividualStarIMFUpperMassCutoff){

            // star formation may be possible
            // compute values and check jeans mass unstable

            // AJE: Feb 2017 - Apparently the dm field is only computed for dm particles, and not
            //                 computed for the static background used in the isolated galaxy sims.
            //                 This *shouldn't* be an issue if the density threshold is high, as the SF
            //                 regions should be dom by self-gravity and the local DM density should be
            //                 much less than the local baryon density... this should be fixed
            //                 if used in low resolution simulations

            dtot = ( BaryonField[DensNum][index] + dm[index] ) * (DensityUnits);         // total density
            tdyn = sqrt(3.0 * pi / 32.0 / GravConst / dtot) / (TimeUnits);            // in code units

            if (StarMakerUseJeansMass){
              isosndsp2 = sndspdC * temp[index] ;
              jeansmass = pi / (6.0 * sqrt(BaryonField[DensNum][index]*DensityUnits) *
                              POW(pi * isosndsp2 / GravConst ,1.5)) / SolarMass; // in solar masses
            } else{
              jeansmass = -1.0; // so always < bmass
            }

            float vel_div = -1.0;
            if (IndividualStarCheckVelocityDiv){
              if (HydroMethod == 2){
                vel_div = BaryonField[Vel1Num][index + xo] - BaryonField[Vel1Num][index] +
                          BaryonField[Vel2Num][index + yo] - BaryonField[Vel2Num][index] +
                          BaryonField[Vel3Num][index + zo] - BaryonField[Vel3Num][index];
              } else{
                vel_div = BaryonField[Vel1Num][index + xo] - BaryonField[Vel1Num][index - xo] +
                          BaryonField[Vel2Num][index + yo] - BaryonField[Vel2Num][index - yo] +
                          BaryonField[Vel3Num][index + zo] - BaryonField[Vel3Num][index - zo];
              }
            }

            if (jeansmass <= bmass && vel_div < 0.0){
              float lowest_cell_mass = bmass;
              bmass = 0.0; number_of_sf_cells = 0;

              int istart, iend, jstart, jend, kstart, kend;

              // compute indeces for adjacent cells integer_sep away from center in each dir
              // stop at grid edge if cell near boundary
              istart = iend = jstart = jend = kstart = kend = 0;
              if (integer_sep > 0){
                istart   = min( i - ibuff             , integer_sep);
                iend     = min( (nx - ibuff - 1 ) - i, integer_sep);
                jstart   = min( j - ibuff             , integer_sep);
                jend     = min( (ny - ibuff - 1 ) - j, integer_sep);
                kstart   = min( k - ibuff             , integer_sep);
                kend     = min( (nz - ibuff - 1 ) - k, integer_sep);
              }

              // loop through cells and add up total amount of mass available for SF
              int l = 0;
              for (int k_loc = -kstart; k_loc <= kend; k_loc++){
                for(int j_loc = -jstart; j_loc <= jend; j_loc++){
                  for (int i_loc = -istart; i_loc <= iend; i_loc++){
                    int loc_index = (i + i_loc) + ( (j + j_loc) + (k + k_loc)*(ny))*(nx);

                    if(BaryonField[DensNum][loc_index] > secondary_odthreshold &&
                       temp[loc_index] < IndividualStarTemperatureThreshold       ){
                      float current_cell_mass = BaryonField[DensNum][loc_index]*(dx*dx*dx)*MassUnits/SolarMass;

                      bmass += current_cell_mass; // solar masses
                      metal_mass += BaryonField[MetalNum][loc_index] * current_cell_mass;

                      if (MultiSpecies > 1){
                        H2mass += (BaryonField[H2INum][loc_index] + BaryonField[H2IINum][loc_index])*\
                                     current_cell_mass;
                      }

                      number_of_sf_cells++;
                      if ( current_cell_mass < lowest_cell_mass){
                        lowest_cell_mass = current_cell_mass;
                      }
                    }

                    /* if PPM, need to be careful about energies */
                    if (HydroMethod != 2 && FALSE){
                      ke_before[ l ] = 0.5 * BaryonField[DensNum][loc_index] *
                                 ( BaryonField[Vel1Num][loc_index] * BaryonField[Vel1Num][loc_index] +
                                   BaryonField[Vel2Num][loc_index] * BaryonField[Vel2Num][loc_index] +
                                   BaryonField[Vel3Num][loc_index] * BaryonField[Vel3Num][loc_index]);
                    }
                    l++;
                  }
                }
              }


              // only allow star formation if IMF can fully sampled (safely) in the given region
              //   aka, make sure there is enough mass to form the most massive star + some fudge
              if ( bmass * IndividualStarMassFraction < IndividualStarIMFUpperMassCutoff){
                break;
              }

              int add_unresolved_star = FALSE; // only used when IMF mass floor is < imf lower limit
              int form_popIII_stars   = 0;

              /* Here is where star formation actually occurs */
              if( bmass*IndividualStarMassFraction > IndividualStarSFGasMassThreshold ){
                // if true, we can try and form stars. compute probability that this mass will
                // form stars this timestep
                star_fraction  = min(StarMakerMassEfficiency*(this->dtFixed)/tdyn, 1.0);
                mass_to_stars  = star_fraction * bmass;

                pstar          = mass_to_stars / IndividualStarSFGasMassThreshold;

                rnum           = (float) (random() % max_random) / ( (float) max_random);

                if ( rnum < pstar){ // form stars until mass runs out - keep star if too much is made

                  form_popIII_stars = ((IndividualStarPopIIIFormation) *\
                                       ((metal_mass/bmass) <= PopIIIMetalCriticalFraction) ); // critial metallicity

                  if (form_popIII_stars){

                    if ( (H2mass/bmass) > PopIIIH2CriticalFraction){ // must check this separately
                      float mass_counter    = IndividualStarSFGasMassThreshold;
                      float unresolved_mass = 0.0;

                      while( mass_counter > 0.0){
                        float temp_mass = SamplePopIII_IMF();

                        ParticleMass[ii] = temp_mass;
                        ParticleType[ii] = -PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII;
                        ii++;
                        sum_mass += temp_mass;
                        mass_counter -= temp_mass;
                      }
                    } // else, do not form any stars!!

                  } else {
                    float mass_counter = IndividualStarSFGasMassThreshold;
                    float unresolved_mass = 0.0;
                    while( mass_counter > 0.0){
                        float temp_mass = SampleIMF();

                        // if the mass is between the lower mass cutoff and the mass floor, sum the total mass of
                        // these particles and dump into a single particle at the very end.
                        if ( (temp_mass >= IndividualStarIMFLowerMassCutoff) && (temp_mass < IndividualStarIMFMassFloor)){
                          unresolved_mass  += temp_mass;
                        } else{
                          ParticleMass[ii]  = temp_mass;
                          ParticleType[ii]  = -PARTICLE_TYPE_INDIVIDUAL_STAR;
                          ii++;
                        }
                        sum_mass         += temp_mass;
                        mass_counter     -= temp_mass;
                    }

                    if (unresolved_mass > 0.0) { // we've formed tiny stars (should always happen).. account for this
                      add_unresolved_star = TRUE;
                      ParticleMass[ii] = unresolved_mass;
                      ParticleType[ii] = -PARTICLE_TYPE_INDIVIDUAL_STAR_UNRESOLVED;
                      ii++;
                    }
                  } // check if doing popIII or normal SF

                } // endif randum number draw check

              } // endif mass threshold check

              // prepare for assigning star properties by computing the local
              // gas velocity properties (this is for velocity assignment)
              // 2 = Zeus .. otherwise PPM
              // copied from pop3_maker.F
              if (HydroMethod == 2){
                umean = (
                       0.5 * (BaryonField[Vel1Num][index   ] + BaryonField[Vel1Num][index+xo])*BaryonField[DensNum][index] +
                       0.5 * (BaryonField[Vel1Num][index-xo] + BaryonField[Vel1Num][index   ])*BaryonField[DensNum][index-xo] +
                       0.5 * (BaryonField[Vel1Num][index+xo] + BaryonField[Vel1Num][index + xo + xo])*BaryonField[DensNum][index+xo] +
                       0.5 * (BaryonField[Vel1Num][index+yo] + BaryonField[Vel1Num][index + xo + yo])*BaryonField[DensNum][index+yo] +
                       0.5 * (BaryonField[Vel1Num][index-yo] + BaryonField[Vel1Num][index + xo - yo])*BaryonField[DensNum][index-yo] +
                       0.5 * (BaryonField[Vel1Num][index+zo] + BaryonField[Vel1Num][index + xo + zo])*BaryonField[DensNum][index+zo] +
                       0.5 * (BaryonField[Vel1Num][index-zo] + BaryonField[Vel1Num][index + xo - zo])*BaryonField[DensNum][index-zo]) /
                      ( BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                        BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                        BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo] ); //
                vmean = (
                         0.5 * (BaryonField[Vel2Num][index   ] + BaryonField[Vel2Num][index+xo])*BaryonField[DensNum][index] +
                         0.5 * (BaryonField[Vel2Num][index-xo] + BaryonField[Vel2Num][index   ])*BaryonField[DensNum][index-xo] +
                         0.5 * (BaryonField[Vel2Num][index+xo] + BaryonField[Vel2Num][index + xo + xo])*BaryonField[DensNum][index+xo] +
                         0.5 * (BaryonField[Vel2Num][index+yo] + BaryonField[Vel2Num][index + xo + yo])*BaryonField[DensNum][index+yo] +
                         0.5 * (BaryonField[Vel2Num][index-yo] + BaryonField[Vel2Num][index + xo - yo])*BaryonField[DensNum][index-yo] +
                         0.5 * (BaryonField[Vel2Num][index+zo] + BaryonField[Vel2Num][index + xo + zo])*BaryonField[DensNum][index+zo] +
                         0.5 * (BaryonField[Vel2Num][index-zo] + BaryonField[Vel2Num][index + xo - zo])*BaryonField[DensNum][index-zo]) /
                        ( BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                         BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                          BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo] ); //
                wmean = (
                         0.5 * (BaryonField[Vel3Num][index   ] + BaryonField[Vel3Num][index+xo])*BaryonField[DensNum][index] +
                         0.5 * (BaryonField[Vel3Num][index-xo] + BaryonField[Vel3Num][index   ])*BaryonField[DensNum][index-xo] +
                         0.5 * (BaryonField[Vel3Num][index+xo] + BaryonField[Vel3Num][index + xo + xo])*BaryonField[DensNum][index+xo] +
                         0.5 * (BaryonField[Vel3Num][index+yo] + BaryonField[Vel3Num][index + xo + yo])*BaryonField[DensNum][index+yo] +
                         0.5 * (BaryonField[Vel3Num][index-yo] + BaryonField[Vel3Num][index + xo - yo])*BaryonField[DensNum][index-yo] +
                         0.5 * (BaryonField[Vel3Num][index+zo] + BaryonField[Vel3Num][index + xo + zo])*BaryonField[DensNum][index+zo] +
                         0.5 * (BaryonField[Vel3Num][index-zo] + BaryonField[Vel3Num][index + xo - zo])*BaryonField[DensNum][index-zo]) /
                        ( BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                          BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                          BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo] ); //
              }
              else{ // PPM case
                umean = (BaryonField[Vel1Num][index]*BaryonField[DensNum][index] +
                              BaryonField[Vel1Num][index-xo]*BaryonField[DensNum][index-xo] +
                              BaryonField[Vel1Num][index+xo]*BaryonField[DensNum][index+xo] +
                              BaryonField[Vel1Num][index-yo]*BaryonField[DensNum][index-yo] +
                              BaryonField[Vel1Num][index+yo]*BaryonField[DensNum][index+yo] +
                              BaryonField[Vel1Num][index+zo]*BaryonField[DensNum][index+zo] +
                              BaryonField[Vel1Num][index-zo]*BaryonField[DensNum][index-zo] ) /
                              (BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                               BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                               BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo]);
                vmean = (BaryonField[Vel2Num][index]*BaryonField[DensNum][index] +
                              BaryonField[Vel2Num][index-xo]*BaryonField[DensNum][index-xo] +
                              BaryonField[Vel2Num][index+xo]*BaryonField[DensNum][index+xo] +
                              BaryonField[Vel2Num][index-yo]*BaryonField[DensNum][index-yo] +
                              BaryonField[Vel2Num][index+yo]*BaryonField[DensNum][index+yo] +
                              BaryonField[Vel2Num][index+zo]*BaryonField[DensNum][index+zo] +
                              BaryonField[Vel2Num][index-zo]*BaryonField[DensNum][index-zo] ) /
                              (BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                               BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                               BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo]);

                wmean = (BaryonField[Vel3Num][index]*BaryonField[DensNum][index] +
                              BaryonField[Vel3Num][index-xo]*BaryonField[DensNum][index-xo] +
                              BaryonField[Vel3Num][index+xo]*BaryonField[DensNum][index+xo] +
                              BaryonField[Vel3Num][index-yo]*BaryonField[DensNum][index-yo] +
                              BaryonField[Vel3Num][index+yo]*BaryonField[DensNum][index+yo] +
                              BaryonField[Vel3Num][index+zo]*BaryonField[DensNum][index+zo] +
                              BaryonField[Vel3Num][index-zo]*BaryonField[DensNum][index-zo] ) /
                              (BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                               BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                               BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo]);
              } // imethod velocity computation


              // now assign particle properties, loop over every star
              px = 0.0; py = 0.0; pz =0.0; // initialize momentum counters
              for (istar = index_presf; istar < ii; istar++){

                ParticleAttribute[0][istar]    = this->Time;                        // formation time
                ParticleAttribute[2][istar]    = metal_mass / bmass ; //BaryonField[MetalNum][index]; // metal fraction (conv from density in Grid_StarParti$

                if (ParticleType[istar] == -PARTICLE_TYPE_INDIVIDUAL_STAR){

                  if(IndividualStarInterpolateLifetime(ParticleAttribute[1][istar], ParticleMass[istar],
                                                                                    ParticleAttribute[2][istar], 1) == FAIL){
                    printf(" %"ESYM"  %"ESYM"  %"ESYM"\n",ParticleAttribute[1][istar], ParticleMass[istar], ParticleAttribute[2][istar]);
                  ENZO_FAIL("Error in stellar lifetime interpolation");
                  }

                } else if (ParticleType[istar] == -PARTICLE_TYPE_INDIVIDUAL_STAR_UNRESOLVED){

                  ParticleAttribute[1][istar] = huge_number * TimeUnits; // make sure its VERY large

                } else if (ParticleType[istar] == -PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII){

                  float temp_mass, temp_lifetime, temp_luminosity;
                  temp_mass = 1.0 * ParticleMass[istar];

                  FORTRAN_NAME(pop3_properties)(&temp_mass, &temp_luminosity, &temp_lifetime);

                  ParticleAttribute[1][istar] = temp_lifetime * yr_s; // in seconds

                } // end check for particle type for assigning lifetimes

                ParticleAttribute[1][istar] /= TimeUnits; // convert from s to code units


                ParticleAttribute[3][istar]    = ParticleMass[istar]; //progenitor mass in solar (main sequence mass)
                ParticleMass[istar]            = ParticleMass[istar] * SolarMass / MassUnits;   // mass in code (not yet dens)

               // give the star particle a position chosen at random
                // within the cell ( so they are not all at cell center )

                rnum =  (float) (random() % max_random) / (float) (max_random);
                ParticlePosition[0][istar] = this->CellWidth[0][i]*rnum + this->CellLeftEdge[0][i];

                rnum =  (float) (random() % max_random) / (float) (max_random);
                ParticlePosition[1][istar] = this->CellWidth[1][j]*rnum + this->CellLeftEdge[1][j];

                rnum =  (float) (random() % max_random) / (float) (max_random);
                ParticlePosition[2][istar] = this->CellWidth[2][k]*rnum + this->CellLeftEdge[2][k];

                // assume velocity dispersion is isotropic in each velocity component. Multiply disp by
                // sqrt(1/3) to get disp in each component... taking above velocities as the mean
                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[0][istar] = umean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*km_cm*(TimeUnits)/(LengthUnits);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[1][istar] = vmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*km_cm*(TimeUnits)/(LengthUnits);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[2][istar] = wmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*km_cm*(TimeUnits)/(LengthUnits);

                // ENSURE MOMENTUM CONSERVATION!!!!!
                // make running total of momentum in each direction
                px += ParticleVelocity[0][istar]*ParticleMass[istar];
                py += ParticleVelocity[1][istar]*ParticleMass[istar];
                pz += ParticleVelocity[2][istar]*ParticleMass[istar];

                // We did the metallicity tagging already, but now loop through and
                // do individual chemical tagging for each species tracked in the simulation
                // these are stored as particle attributes starting with attr number 5 (index 4)
                if (MultiMetals == 2){
//                  if (IndividualStarOutputChemicalTags){
                    // output the particle formation time, birth mass (in SolarMass), and metallicity
//                    printf(" %"ISYM" %"ESYM" %"ESYM" %"ESYM, ParticleType[istar],
//                                                     ParticleAttribute[0][istar], ParticleAttribute[3][istar],
//                                                     ParticleAttribute[2][istar]);
//                  }
                  int iyield = 0;
                  for(iyield = 0; iyield < StellarYieldsNumberOfSpecies; iyield++){
                    double temp_fraction = 0.0;

                    if(StellarYieldsAtomicNumbers[iyield] > 2){
                      int field_num;
                      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[iyield]);

                      if (IndividualStarOutputChemicalTags){
//                        printf(" %"ESYM, BaryonField[field_num][index]);
                        StellarAbundances[iyield][istar]     = BaryonField[field_num][index];
                      } else {
                        ParticleAttribute[4 + iyield][istar] = BaryonField[field_num][index];
                      }

                    } else if (StellarYieldsAtomicNumbers[iyield] == 1){
                      // Take H and He fractions as TOTAL amount of H in the cell
                      // this is probably not needed since it should all be HI in a star forming region anyway
                      temp_fraction = BaryonField[HINum][index] + BaryonField[HIINum][index];

                      if (MultiSpecies > 1){
                        temp_fraction += BaryonField[HMNum][index] +
                                         BaryonField[H2INum][index] + BaryonField[H2IINum][index];
                      }

                      if (IndividualStarOutputChemicalTags){
//                        printf(" %"ESYM, temp_fraction);
                        StellarAbundances[iyield][istar]     = temp_fraction;
                      } else {
                        ParticleAttribute[4 + iyield][istar] = temp_fraction;
                      }
                    } else if (StellarYieldsAtomicNumbers[iyield] == 2){
                      // Again, total amount of Helium - probably not necessary, should all be HeI anyway
                      temp_fraction = BaryonField[HeINum][index] + BaryonField[HeIINum][index] +
                                      BaryonField[HeIIINum][index];

                      if (IndividualStarOutputChemicalTags){
//                        printf(" %"ESYM, temp_fraction);
                        StellarAbundances[iyield][istar]     = temp_fraction;
                      } else{
                        ParticleAttribute[4 + iyield][istar] = temp_fraction;
                      }
                    }
                  }// loop over yields

                  if (IndividualStarOutputChemicalTags){
                    int offset = 0;
                    // add two to the stellar abundances for popIII and AGB metals (if tracked)
                    if (IndividualStarTrackAGBMetalDensity){
                      if (IndividualStarOutputChemicalTags){
                        StellarAbundances[StellarYieldsNumberOfSpecies + 0][istar] = BaryonField[AGBMetalNum][index];
                      } else {
                        ParticleAttribute[4 + iyield + offset][istar] = BaryonField[AGBMetalNum][index];
                      }
                      offset++;
                    }

                    if (IndividualStarPopIIIFormation){
                      if (IndividualStarOutputChemicalTags){
                        StellarAbundances[StellarYieldsNumberOfSpecies + offset][istar] = BaryonField[PopIIIMetalNum][index];
                        offset++;
                        StellarAbundances[StellarYieldsNumberOfSpecies + offset][istar] = BaryonField[PopIIIPISNeMetalNum][index];
                      } else {
                        ParticleAttribute[4 + iyield + offset][istar] = BaryonField[PopIIIMetalNum][index];
                        offset++;
                        ParticleAttribute[4 + iyield + offset][istar] = BaryonField[PopIIIPISNeMetalNum][index];
                      }
                      offset++;
                    }

                    if (IndividualStarTrackSNMetalDensity) {
                      if (IndividualStarOutputChemicalTags){
                        StellarAbundances[StellarYieldsNumberOfSpecies + offset][istar] = BaryonField[SNIaMetalNum][index];
                        offset++;
                        StellarAbundances[StellarYieldsNumberOfSpecies + offset][istar] = BaryonField[SNIIMetalNum][index];
                      } else {
                        ParticleAttribute[4 + iyield + offset][istar] = BaryonField[SNIaMetalNum][index];
                        offset++;
                        ParticleAttribute[4 + iyield + offset][istar] = BaryonField[SNIIMetalNum][index];
                      }
                      offset++;
                    }

                    if (IndividualStarRProcessModel) {
                      if (IndividualStarOutputChemicalTags){
                        StellarAbundances[StellarYieldsNumberOfSpecies + offset][istar] = BaryonField[RProcMetalNum][index];
                      } else {
                        ParticleAttribute[4 + iyield + offset][istar] = BaryonField[RProcMetalNum][index];
                      }
                    }

                  } // if output

//                  if (IndividualStarOutputChemicalTags) printf("\n");
                } // check multimetals

                if (IndividualStarSaveTablePositions && (ParticleType[istar] == -IndividualStar)){
                  int tstart = ParticleAttributeTableStartIndex;

                  // stellar evolution table (attr 3 = birth mass, attr 2 = metallicity)
                  int t_i = -1, t_j = -1, t_k = -1;
                  IndividualStarGetSETablePosition(t_i, t_j,
                                                   ParticleAttribute[3][istar], ParticleAttribute[2][istar]);
                  ParticleAttribute[tstart    ][istar] = t_i;
                  ParticleAttribute[tstart + 1][istar] = t_j;
                  // radiation properties table (only do if particle can radiate - saves time)
                  if( ParticleAttribute[3][istar] >= IndividualStarRadiationMinimumMass){
                    float Teff, R;
                    IndividualStarInterpolateProperties(Teff, R, (int)ParticleAttribute[tstart][istar],
                                                       (int)ParticleAttribute[tstart+1][istar],
                                                       ParticleAttribute[3][istar], ParticleAttribute[2][istar]);
                    float g = IndividualStarSurfaceGravity(ParticleAttribute[3][istar], R);

                    t_i = -1; t_j = -1; t_k = -1;
                    IndividualStarGetRadTablePosition(t_i, t_j, t_k,
                                                    Teff, g, ParticleAttribute[2][istar]);
                    ParticleAttribute[tstart + 2][istar] = t_i;
                    ParticleAttribute[tstart + 3][istar] = t_j;
                    ParticleAttribute[tstart + 4][istar] = t_k;
                  } else {
                    ParticleAttribute[tstart + 2][istar] = -1;
                    ParticleAttribute[tstart + 3][istar] = -1;
                    ParticleAttribute[tstart + 4][istar] = -1;
                  }

                  // yields table position
                  t_i = -1 ; t_j = -1;
                  StellarYieldsGetYieldTablePosition(t_i, t_j,
                                                   ParticleAttribute[3][istar], ParticleAttribute[2][istar]);
                  ParticleAttribute[tstart + 5][istar] = t_i;
                  ParticleAttribute[tstart + 6][istar] = t_j;

                } // end check for saving table positions

                /* Keeping this as a particle attribute */
                ParticleAttribute[NumberOfParticleAttributes-2][istar] = 0.0; // wind mass ejected
                ParticleAttribute[NumberOfParticleAttributes-1][istar] = 0.0; // sn mass ejected


              } // end while loop for assigning particle properties

              // ensure zero net momentum from mean velocity
              // momentum of gas converted into stars (sum_mass * umean)
              // should be equal to the total momentum of the stars
              // compute excess momentum and modify star velocity evenly (mass weighted)
              // this is not completely physical, as pre-SF and post-SF gas vel is the same
              sum_mass = sum_mass * SolarMass / MassUnits; // in code units

              px_excess = px - umean * sum_mass;
              py_excess = py - vmean * sum_mass;
              pz_excess = pz - wmean * sum_mass;

              px_excess = py_excess = pz_excess = 0.0; // TURNED OFF

              // remove or add momentum evenly from each star if needed
/*
              if ( abs(px_excess) > tiny_number) {
                for (istar = index_presf; istar < ii; istar++){
                  ParticleVelocity[0][istar] += (-1.0 * px_excess) / (ParticleMass[istar] * (float) (ii-index_presf));
                }
              }
              if ( abs(py_excess) > tiny_number){
                for (istar = index_presf; istar < ii; istar++){
                  ParticleVelocity[1][istar] = (-1.0 * py_excess) / (ParticleMass[istar] * (float) (ii-index_presf));
                }
              }
              if ( abs(pz_excess) > tiny_number){
                for (istar = index_presf; istar < ii; istar++){
                  ParticleVelocity[2][istar] = (-1.0 * pz_excess) / (ParticleMass[istar] * (float) (ii-index_presf));
                }
              }
*/
              // now remove mass from grid - do not need to do this for tracer fields since they are kept as fractions
              // and will be modified accordingly when converted back to densities in Grid_StarParticleHandler
              l = 0;
              for (int k_loc = -kstart; k_loc <= kend; k_loc++){
                for(int j_loc = -jstart; j_loc <= jend; j_loc++){
                  for (int i_loc = -istart; i_loc <= iend; i_loc++){
                    int loc_index = (i + i_loc) + ( (j + j_loc) + (k + k_loc)*(ny))*(nx);

                    if(BaryonField[DensNum][loc_index] > secondary_odthreshold &&
                       temp[loc_index] < IndividualStarTemperatureThreshold    ){
                      // mass is removed as weighted by the previous cell mass (more mass is
                      // taken out of higher density regions). M_new = M_old - M_sf * (M_old / M_tot)
                      // where M_tot is mass of cells that meet above SF conditions. Simplifies to below eq:
                      float old_density = BaryonField[DensNum][loc_index];

                      BaryonField[DensNum][loc_index] *= (1.0 - sum_mass / bmass);

                      // adjust total energy if we are using PPM
                      if (HydroMethod != 2 && FALSE){
                          float ke_after, delta_ke;

                          ke_after = 0.5 * BaryonField[DensNum][loc_index] *
                                   ( BaryonField[Vel1Num][loc_index] * BaryonField[Vel1Num][loc_index] +
                                     BaryonField[Vel2Num][loc_index] * BaryonField[Vel2Num][loc_index] +
                                     BaryonField[Vel3Num][loc_index] * BaryonField[Vel3Num][loc_index]);

                          delta_ke = ke_after - ke_before[ l ];

                           // or TE_new = TE_old/(1.0 -sum_mass / bmass) + delta_ke / new_density;
                          BaryonField[TENum][loc_index] = (BaryonField[TENum][loc_index]*old_density +
                                                                  delta_ke) / BaryonField[DensNum][loc_index];
                      }
                      l++;
                    }
                  }
                }
              }

            } // if jeans mass unstable
          } // resolution and density


        } // enx x loop
      } // end y loop
    } // end z loop

    // Done forming stars!!! Output and exit
    if (ii > 0){
      printf("P(%"ISYM"): individual_star_maker[add]: %"ISYM" new star particles\n", MyProcessorNumber, ii);
    }
    if (ii >= *nmax){
      fprintf(stdout, "individual_star_maker: reached max new particle count!! Available: %"ISYM". Made: %"ISYM"\n", *nmax, ii);

    }

  // star masses are recorded as densities (mass / cell volume)
  // set progenitor masses in solar
  for (int counter = 0; counter < ii; counter++){
    ParticleMass[counter]   = ParticleMass[counter] / (dx*dx*dx); // code units / cell volume
  }

  *np = ii; // number of stars formed : AJE 2/29 check if this is a bug with the -1


  delete[] ke_before;

  return SUCCESS;
}


float SamplePopIII_IMF(void){

  return SampleIMF(SecondaryIMFData,
                   PopIIILowerMassCutoff, PopIIIUpperMassCutoff);
}

float SampleIMF(void){
// to retain default behavoir

  return SampleIMF(IMFData,
                   IndividualStarIMFLowerMassCutoff,
                   IndividualStarIMFUpperMassCutoff);

}

float SampleIMF(float * data, const float & lower_mass, const float & upper_mass){
/*-----------------------------------------------------------------------------
  SampleIMF

  Samples the tabulated initial mass function a la cumulative probablity density
  as created in StarParticleIndividual_IMFInitialize.

  INPUTS -
    None

  OUTPUTS -
    Mass of randomly selected star in solar masses

-----------------------------------------------------------------------------*/

  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  float x = (float) (random_int%max_random) / (float) (max_random);
  float dm = log10(upper_mass / lower_mass )/ ((float) (IMF_TABLE_ENTRIES-1));
  float m;

  int bin_number;

  if (x <= data[0] ){
    bin_number = 0;
  } else if (x >= data[IMF_TABLE_ENTRIES-1]){
    bin_number = IMF_TABLE_ENTRIES - 1;
  } else{
    bin_number = search_lower_bound(data, x, 0, IMF_TABLE_ENTRIES, IMF_TABLE_ENTRIES);
  }

  m = lower_mass * POW(10.0, bin_number * dm);

  IndividualStarIMFCalls++;

  return m;
}
