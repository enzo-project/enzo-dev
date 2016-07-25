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

float SampleIMF(void);

float ComputeSnIaProbability(const float &current_time, const float &formation_time, const float &lifetime, const float &TimeUnits);
unsigned_long_int mt_random(void);

void ComputeStellarWindVelocity(const float &mproj, const float &metallicity,
                                const float &lifetime,  float *v_wind);

void ComputeStellarWindMassLossRate(const float &mproj, const float &metallicity,
                                   const float & lifetime, float *dMdt);


void CheckFeedbackCellCenter(const FLOAT &xp, const FLOAT &yp, const FLOAT &zp,
                             const FLOAT &xstart, const FLOAT &ystart, const FLOAT &zstart,
                             const FLOAT &dx,
                             const int &nx, const int &ny, const int &nz, const int &ibuff,
                             FLOAT *xfc, FLOAT *yfc, FLOAT *zfc);

void Momentum(float *u, float *v, float *w, float *d,
              const float &up, const float &vp, const float &wp,
              const int &nx, const int &ny, const int &nz,
              const int &ic, const int &jc, const int &kc,
              const int &iface, const int &jface, const int &kface,
              const int &stencil, int idir);

void MetalConversion(float *metal_field, float *d, const float &dx,
                     const int &nx, const int &ny, const int &nz,
                     const int &ic, const int &jc, const int &kc,
                     const int &stencil, int idir);

void SumMassEnergy(float *pu, float *pv, float *pw, float *d, float *ge, float *te,
                   const int &nx, const int &ny, const int &nz,
                   const int &iface, const int &jface, const int &kface,
                   const int &ic, const int &jc, const int &kc, const int &stencil,
                   float *mass_sum, float *energy_sum, float *kin_energy_sum);

void ComputeAbcCoefficients(float *pu, float *pv, float *pw, float *d,
                            float *ge, float *pu_l, float *pv_l, float *pw_l,
                            float *d_l, const int &nx, const int &ny, const int &nz,
                            const int &iface, const int &jface, const int &kface,
                            const int &ic, const int &jc, const int &kc,
                            const int &stencil,
                            float &A, float &B, float &C);

void AddFeedbackToGridCells(float *pu, float *pv, float *pw, 
                            float *d, float *ge, float *te,
                            const int &nx, const int &ny, const int &nz,
                            const int &ic, const int &jc, const int &kc,
                            const int &iface, const int &jface, const int &kface,
                            const float &dxf, const float &dyf, const float &dzf,
                            const float &dxc, const float &dyc, const float &dzc,
                            const float &mass_per_cell, const float &mom_per_cell,
                            const float &therm_per_cell, const int &stencil);

void AddMetalSpeciesToGridCells(float *m, const float &mass_per_cell,
                                const int &nx, const int &ny, const int &nz,
                                const int &ic, const int &jc, const int &kc,
                                const float &dxc, const float &dyc, const float &dzc,
                                const int &stencil);




int grid::individual_star_maker(float *dm, float *temp, int *nmax, float *mu, int *np,
                                float *ParticleMass,
                                int *ParticleType, FLOAT *ParticlePosition[],
                                float *ParticleVelocity[], float *ParticleAttribute[]){
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

  const double msolar  = 1.989e33;
  const double sndspdC = 1.3095e8;
  const double myr     = 3.1536e13;

  int i, j, k, index, ii=0, istar=0, index_presf=0;
  int xo, yo, zo, rsign=1;
  float bmass, div, star_mass=0.0, sum_mass=0.0;
  float pstar, mass_to_stars, mass_available, tdyn;
  float dtot, isosndsp2, jeansmass, star_fraction;
  float umean, vmean, wmean, px, py, pz, px_excess, py_excess, pz_excess;
  float rnum;
  const double m_h = 1.673e-24;
  float inv_metal_mol = 1.0 /MU_METAL;

  const int max_random = (1<<16);

  int form_method = -1; // tracker for debugging purposes

  /* for convenience, rename some grid properties - will likely get optimized out */
  int  nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  FLOAT xstart = CellLeftEdge[0][0], ystart = CellLeftEdge[1][0], zstart = CellLeftEdge[2][0];
  float   dx   = CellWidth[0][0];


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
  if (MultiSpecies == TRUE){
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  /* get metallicity tracer field number */
  int MetalNum;
  MetalNum   = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);

  /* get units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units


  if (ProblemType == 260){ // place a star by hand and exit

    // Only form on first call of individual_star_maker
    if(ChemicalEvolutionTestStarFormed){
      return SUCCESS;
    } else {
      FLOAT xx, yy, zz;
      xx = ChemicalEvolutionTestStarPosition[0];
      yy = ChemicalEvolutionTestStarPosition[1];
      zz = ChemicalEvolutionTestStarPosition[2];

      // make sure particle position is on this grid / processor
      if( !( (xx > this->CellLeftEdge[0][ibuff - 1]) && (xx < this->CellLeftEdge[0][nx - ibuff ] )) ||
          !( (yy > this->CellLeftEdge[1][ibuff - 1]) && (yy < this->CellLeftEdge[1][ny - ibuff ] )) ||
          !( (zz > this->CellLeftEdge[2][ibuff -1 ]) && (zz < this->CellLeftEdge[2][nz - ibuff ] )) ) {
        ChemicalEvolutionTestStarFormed = TRUE; // setting this here to avoid doing MPI communication
                                                // on whatever processor the star actually gets placed
        printf("P(%"ISYM") individual_star_maker: Particle not on this grid. Leaving\n", MyProcessorNumber);
        return SUCCESS;
      }

      // deposit the star by hand
      ParticleMass[0] = ChemicalEvolutionTestStarMass * msolar / MassUnits / (dx*dx*dx);
      ParticleType[0] = - PARTICLE_TYPE_INDIVIDUAL_STAR;
      ParticleAttribute[0][0] = this->Time;

      // allow user to set lifetime artificially
      if(ChemicalEvolutionTestStarLifetime > 0){
        ParticleAttribute[1][0] = ChemicalEvolutionTestStarLifetime * myr / (TimeUnits);
      } else{
        // last arg tells function to return total stellar lifetime
        if(IndividualStarInterpolateLifetime(ParticleAttribute[1][0], ChemicalEvolutionTestStarMass,
                                                                      ChemicalEvolutionTestStarMetallicity, 1) == FAIL){
            ENZO_FAIL("Failure in stellar lifetime interpolation");
        }

        ParticleAttribute[1][0] /= TimeUnits; // convert from s to code units

      }

      ParticleAttribute[3][0] = ChemicalEvolutionTestStarMass; // in solar!!!
      ParticleAttribute[2][0] = ChemicalEvolutionTestStarMetallicity;

      ParticlePosition[0][0] = ChemicalEvolutionTestStarPosition[0];
      ParticlePosition[1][0] = ChemicalEvolutionTestStarPosition[1];
      ParticlePosition[2][0] = ChemicalEvolutionTestStarPosition[2];

      ParticleVelocity[0][0] = ChemicalEvolutionTestStarVelocity[0]*1.0E5 / VelocityUnits;
      ParticleVelocity[1][0] = ChemicalEvolutionTestStarVelocity[1]*1.0E5 / VelocityUnits;
      ParticleVelocity[2][0] = ChemicalEvolutionTestStarVelocity[2]*1.0E5 / VelocityUnits;

      // find grid cell and assign chemical tags
      int ip, jp, kp, n;
      ip = int ( (ParticlePosition[0][0] - (xstart)) / (dx));
      jp = int ( (ParticlePosition[1][0] - (ystart)) / (dx));
      kp = int ( (ParticlePosition[2][0] - (zstart)) / (dx));

      n  = ip + (jp + kp * (ny)) * (nx);

      /* Metal fields are all in fractions, as set in Grid_StarParticleHandler */
      if(TestProblemData.MultiMetals == 2){
        for( int ii = 0; ii < StellarYieldsNumberOfSpecies; ii++){
          if(StellarYieldsAtomicNumbers[ii] > 2){
            int field_num;

            this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[ii]);

            ParticleAttribute[4 + ii][0] = BaryonField[field_num][n];

          } else if (StellarYieldsAtomicNumbers[ii] == 1){
            /* Take H and He fractions as TOTAL amount of H and He species in the cell */
            ParticleAttribute[4 + ii][0] = BaryonField[HINum][n] + BaryonField[HIINum][n];

          } else if (StellarYieldsAtomicNumbers[ii] == 2){

            ParticleAttribute[4 + ii][0] = BaryonField[HeINum][n]  +
                                           BaryonField[HeIINum][n] + BaryonField[HeIIINum][n];

          }
        }
      }

      *np = 1;
      ChemicalEvolutionTestStarFormed = TRUE;
      printf("individual_star_maker: Formed star ChemicalEvolutionTest. M =  %"FSYM" and Z = %"FSYM". tau = %"ESYM"\n", ParticleMass[0]*(dx*dx*dx)*MassUnits/msolar, ParticleAttribute[2][0], ParticleAttribute[1][0]); 


      return SUCCESS;
    }
  } // end ChemicalEvolutionTest ProblemType check


    // 3D -> 1D index
    xo = 1;
    yo = nx;
    zo = (nx) * (ny);

    // over density threshold in code units
    // if multispecies is off, assumes a value for MU
    float odthreshold, secondary_odthreshold;
    if (MultiSpecies == FALSE){
        odthreshold           = StarMakerOverDensityThreshold * m_h * (*mu) / (DensityUnits); // code density
        secondary_odthreshold = IndividualStarSecondaryOverDensityThreshold * m_h * (*mu) / (DensityUnits);
    }

    // loop over all cells, check condition, form stars stochastically
    ii = 0; index_presf = 0;

    int number_of_sf_cells = 0;
    int integer_sep = ((int) (IndividualStarCreationStencilSize + 1) / 2.0 - 1); // stencil size must be odd number

    float *ke_before;
    if (HydroMethod != 2){
      ke_before = new float[ ((int) POW(IndividualStarCreationStencilSize,3)) ];
    } else{
      ke_before = NULL;
    }

    for (k = ibuff; k < nz - ibuff; k++){
      for (j = ibuff; j < ny - ibuff; j++){
        for (i = ibuff; i < nx - ibuff; i++){

          index = i + ( j + k * (ny)) * nx;

          /* if distributed star formation */
          // check center cell's SF condition, if met, do SF
          /* loop and sum over all*/

           bmass = (BaryonField[DensNum][index]*(dx*dx*dx)) * MassUnits / msolar; // in solar masses
               // perform the following easy checks for SF before proceeding
               // 1) Is density greater than the density threshold?
               // 2) Is temperature < the minimum temperature?
               // 3) Do not allow star formation if minimum star particle mass is
               //    above some fraction of cell mass. This is very unlikely to occur
               //    in intended use case:
               //    (i.e. baryon mass likely always going to be > ~10 solar masses)

          sum_mass = 0.0; index_presf = ii;

          /* Need to compute Mu exactly and recompute threshold */
          if (MultiSpecies == TRUE){
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
            }

            /* Metal field must be present in this star formation scheme */
            number_density += BaryonField[MetalNum][index] * inv_metal_mol;


            number_density *= BaryonField[DensNum][index] / m_h ; // now actual n density

            mu_cell = BaryonField[DensNum][index] / (number_density * m_h);
            odthreshold = StarMakerOverDensityThreshold * (mu_cell) * m_h / (DensityUnits);

            secondary_odthreshold = IndividualStarSecondaryOverDensityThreshold * (mu_cell) * m_h / DensityUnits;
          }


          if (   BaryonField[DensNum][index]      > odthreshold
              && temp[index] <= IndividualStarTemperatureThreshold){
              //&& IndividualStarMassFraction*bmass > IndividualStarIMFLowerMassCutoff
              //&& 0.5*bmass > IndividualStarIMFUpperMassCutoff){



            // star formation may be possible
            // compute values and check jeans mass unstable
            dtot = ( BaryonField[DensNum][index] + dm[index] ) * (DensityUnits);         // total density
            tdyn = sqrt(3.0 * pi / 32.0 / GravConst / dtot) / (TimeUnits);            // in code units
            isosndsp2 = sndspdC * temp[index] ;
            jeansmass = pi / (6.0 * sqrt(BaryonField[DensNum][index]*DensityUnits) *
                            POW(pi * isosndsp2 / GravConst ,1.5)) / msolar; // in solar masses

            if (jeansmass <= bmass){
              float lowest_cell_mass = bmass;
              bmass = 0.0; number_of_sf_cells = 0;

              int istart, iend, jstart, jend, kstart, kend;

              istart = iend = jstart = jend = kstart = kend = 0;
              if (integer_sep > 0){
                istart   = min( i - ibuff             , integer_sep);
                iend     = min( (nx - ibuff - 1 ) - i, integer_sep);
                jstart   = min( j - ibuff             , integer_sep);
                jend     = min( (ny - ibuff - 1 ) - j, integer_sep);
                kstart   = min( k - ibuff             , integer_sep);
                kend     = min( (nz - ibuff - 1 ) - k, integer_sep);
              }

              int l = 0;
              for (int k_loc = -kstart; k_loc <= kend; k_loc++){
                for(int j_loc = -jstart; j_loc <= jend; j_loc++){
                  for (int i_loc = -istart; i_loc <= iend; i_loc++){
                    int loc_index = (i + i_loc) + ( (j + j_loc) + (k + k_loc)*(ny))*(nx);

                    if(BaryonField[DensNum][loc_index] > secondary_odthreshold){
                      float current_cell_mass = BaryonField[DensNum][loc_index]*(dx*dx*dx)*MassUnits/msolar;
                      bmass += current_cell_mass; // solar masses
                      number_of_sf_cells++;
                      if ( current_cell_mass < lowest_cell_mass){
                        lowest_cell_mass = current_cell_mass;
                      }
                    }

                    /* if PPM, need to be careful about energies */
                    if (HydroMethod != 2){
                      ke_before[ l ] = 0.5 * BaryonField[DensNum][index] *
                                 ( BaryonField[Vel1Num][index] * BaryonField[Vel1Num][index] +
                                   BaryonField[Vel2Num][index] * BaryonField[Vel2Num][index] +
                                   BaryonField[Vel3Num][index] * BaryonField[Vel3Num][index]);
                    }
                    l++;
                  }
                }
              }


              float M_max_star;
              if (IndividualStarAllowTruncatedIMF){
                  /* allow IMF to be truncated to safely allow SF in less dense regions. This may be an
                     unphysical thing to do, as it won't necessarily reproduce the input IMF. In reality
                     SF clouds will do this (low mass clouds are less likely to make very massive stars), but
                     this is all convolved together in the input IMF (IMF may vary cloud-by-cloud)........
                     use this with caution */
                  M_max_star =  min(bmass * IndividualStarMassFraction, IndividualStarIMFUpperMassCutoff);
              } else {
                  /* only allow star formation if IMF can fully sampled (safely) in the given region */
                  if ( bmass * IndividualStarMassFraction < IndividualStarIMFUpperMassCutoff){
                    break;
                  }
                  M_max_star = IndividualStarIMFUpperMassCutoff;
              }


              if(IndividualStarSFAlgorithm == 0){
                // calculate mass in cell that can be converted to stars in timestep
                // generally this should be small (comparable to or less than the lower mass
                // cutoff of the IMF)

                star_fraction  = min(StarMakerMassEfficiency*(this->dtFixed)/tdyn, 1.0);
                mass_to_stars  = star_fraction * bmass;
                mass_available = IndividualStarMassFraction * bmass; // AJE - fixed 6/4
                mass_to_stars  = min(mass_to_stars, mass_available);

                // If mass_to_stars greater than available mass, convert
                // all of available mass into stars
                // Frankly this is very unlikely to occur...
                // Tests as of 2/22/16 show NO SF here for at least 10^5 stars in a LMC dwarf galaxy
                if(mass_to_stars >= mass_available){
                  mass_to_stars = mass_available;
                  while( ii < *nmax && mass_to_stars > M_max_star){
                    ParticleMass[ii] = SampleIMF();
                    sum_mass        += ParticleMass[ii]; // counter for mass formed in this cell
                    mass_to_stars   -= ParticleMass[ii]; // reduce available mass
                    ii++;
                  }
                }

                // Tests (as of 2/22/16) show NO SF here for at least the first 10^5 stars
                if (mass_to_stars > M_max_star){
                  while (ii < *nmax && mass_to_stars > M_max_star){
                    ParticleMass[ii]  = SampleIMF();
                    sum_mass         += ParticleMass[ii];
                    mass_to_stars    -= ParticleMass[ii];
                    ii++;
                  }
                }

                // If mass is above IMF lower limit, star formation will happen.
                // Just form stars randomly over IMF until mass dips below lower cutoff
                if(mass_to_stars > IndividualStarIMFLowerMassCutoff){

                  // loop until mass to stars is less than 10% of smallest star particle size
                  while( ii < *nmax && mass_to_stars > 1.1 * IndividualStarIMFLowerMassCutoff){
                    float tempmass;
                    tempmass = SampleIMF();

                    if (tempmass < M_max_star){
                        ParticleMass[ii]  = SampleIMF();
                        sum_mass         += ParticleMass[ii];
                        mass_to_stars    -= ParticleMass[ii];
                        ii++;
                    } // else redraw
  
                    if (mass_to_stars < 0.0){
                      mass_to_stars = 0.0;
                    }
                  }
                } // end mass above individual star cutoff

                // now we are in the Goldbaum et. al. 2015 regime (star_maker_ssn.F)
                // Calculate probability of star forming and form stars stochastically
                if (mass_to_stars < IndividualStarIMFLowerMassCutoff && mass_to_stars > tiny_number){
                  star_mass = SampleIMF();
                  pstar     = mass_to_stars / star_mass;
                  rnum =  (float) (random() % max_random) / ((float) max_random);
                  if (rnum < pstar){
                    ParticleMass[ii]  = star_mass;
                    sum_mass         += ParticleMass[ii];
                    ii++;
                  }
                }

              } else if (IndividualStarSFAlgorithm == 1){
                /* take chunks of mass */

                if( bmass*IndividualStarMassFraction > IndividualStarSFGasMassThreshold ){ // set to ~ 2 x M_max_star
                  // if true, we can try and form stars. compute probability that this mass will
                  // form stars this timestep
                  star_fraction  = min(StarMakerMassEfficiency*(this->dtFixed)/tdyn, 1.0);
                  mass_to_stars  = star_fraction * bmass;

                  pstar          = mass_to_stars / IndividualStarSFGasMassThreshold;

                  rnum           = (float) (random() % max_random) / ( (float) max_random);

                  if ( rnum < pstar){ // form stars until mass runs out - keep star if too much is made
                      float mass_counter = IndividualStarSFGasMassThreshold;
                      while( mass_counter > 0.0){
                          ParticleMass[ii]  = SampleIMF();
                          sum_mass         += ParticleMass[ii];
                          mass_counter     -= ParticleMass[ii];
                          ii++;
                      }

                  } // endif randum number draw check

                } // endif mass threshold check

              } // endif sf algorithm for star formation



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

                ParticleType[istar]            = - PARTICLE_TYPE_INDIVIDUAL_STAR;   // negative is a "new" star
                ParticleAttribute[0][istar]    = this->Time;                        // formation time
                ParticleAttribute[2][istar]    = BaryonField[MetalNum][index]; // metal fraction (conv from density in Grid_StarParticleHandler)


                if(IndividualStarInterpolateLifetime(ParticleAttribute[1][istar], ParticleMass[istar],
                                                                                  ParticleAttribute[2][istar], 1) == FAIL){
                  ENZO_FAIL("Error in stellar lifetime interpolation");
                }
                ParticleAttribute[1][istar] /= TimeUnits; // convert from s to code units

                ParticleAttribute[3][istar]    = ParticleMass[istar]; //progenitor mass in solar (main sequence mass)
                ParticleMass[istar]            = ParticleMass[istar] * msolar / MassUnits;   // mass in code (not yet dens)

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
                ParticleVelocity[0][istar] = umean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0E5*(TimeUnits)/(LengthUnits);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[1][istar] = vmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0e5*(TimeUnits)/(LengthUnits);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[2][istar] = wmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0E5*(TimeUnits)/(LengthUnits);

                // ENSURE MOMENTUM CONSERVATION!!!!!
                // make running total of momentum in each direction
                px += ParticleVelocity[0][istar]*ParticleMass[istar];
                py += ParticleVelocity[1][istar]*ParticleMass[istar];
                pz += ParticleVelocity[2][istar]*ParticleMass[istar];

                // We did the metallicity tagging already, but now loop through and
                // do individual chemical tagging for each species tracked in the simulation
                // these are stored as particle attributes starting with attr number 5 (index 4)
                if(TestProblemData.MultiMetals == 2){
                  for( int iyield = 0; iyield < StellarYieldsNumberOfSpecies; iyield++){
                    if(StellarYieldsAtomicNumbers[ii] > 2){
                      int field_num;

                      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[iyield]);

                      ParticleAttribute[4 + iyield][istar] = BaryonField[field_num][index];

                    } else if (StellarYieldsAtomicNumbers[iyield] == 1){
                      // Take H and He fractions as TOTAL amount of H in the cell
                      // this is probably not needed since it should all be HI in a star forming region anyway
                      ParticleAttribute[4 + iyield][istar] = BaryonField[HINum][index] + BaryonField[HIINum][index];

                    } else if (StellarYieldsAtomicNumbers[iyield] == 2){
                      // Again, total amount of Helium - probably not necessary, should all be HeI anyway
                      ParticleAttribute[4 + iyield][istar] = BaryonField[HeINum][index]  +
                                           BaryonField[HeIINum][index] + BaryonField[HeIIINum][index];

                    }
                  }
                }


              } // end while loop for assigning particle properties
              // ---------------------------------------------------

              // ensure zero net momentum from mean velocity
              // momentum of gas converted into stars (sum_mass * umean)
              // should be equal to the total momentum of the stars
              // compute excess momentum and modify star velocity evenly (mass weighted)
              // this is not completely physical, as pre-SF and post-SF gas vel is the same
              sum_mass = sum_mass * msolar / MassUnits; // in code units


              px_excess = umean * sum_mass + px;
              py_excess = vmean * sum_mass + py;
              pz_excess = wmean * sum_mass + pz;

              // remove or add momentum evenly from each star if needed
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

              // now remove mass from grid - do not need to do this for tracer fields since they are kept as fractions
              // and will be modified accordingly when converted back to densities in Grid_StarParticleHandler
              l = 0;
              for (int k_loc = -kstart; k_loc <= kend; k_loc++){
                for(int j_loc = -jstart; j_loc <= jend; j_loc++){
                  for (int i_loc = -istart; i_loc <= iend; i_loc++){
                    int loc_index = (i + i_loc) + ( (j + j_loc) + (k + k_loc)*(ny))*(nx);

                    if(BaryonField[DensNum][loc_index] > secondary_odthreshold){
                      // mass is removed as weighted by the previous cell mass (more mass is
                      // taken out of higher density regions). M_new = M_old - M_sf * (M_old / M_tot)
                      // where M_tot is mass of cells that meet above SF conditions. Simplifies to below eq:
                      BaryonField[DensNum][loc_index] *= (1.0 - sum_mass / bmass);

                      // adjust total energy if we are using PPM
                      if (HydroMethod != 2){
                          float ke_after, delta_ke;

                          ke_after = 0.5 * BaryonField[DensNum][index] *
                                   ( BaryonField[Vel1Num][index] * BaryonField[Vel1Num][index] +
                                     BaryonField[Vel2Num][index] * BaryonField[Vel2Num][index] +
                                     BaryonField[Vel3Num][index] * BaryonField[Vel3Num][index]);

                          delta_ke = ke_after - ke_before[ l ];

                          BaryonField[TENum][index] += delta_ke / BaryonField[DensNum][index];
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

/*-----------------------------------------------------------------------------
  SampleIMF

  Samples the tabulated initial mass function a la cumulative probablity density
  as created in StarParticleIndividual_IMFInitialize.

  INPUTS -
    None

  OUTPUTS -
    Mass of randomly selected star in solar masses

-----------------------------------------------------------------------------*/
float SampleIMF(void){
  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  float x = (float) (random_int%max_random) / (float) (max_random);
  float dm = log10(IndividualStarIMFUpperMassCutoff / IndividualStarIMFLowerMassCutoff)/ ((float) (IMF_TABLE_ENTRIES-1));
  float m;


  int width = IMF_TABLE_ENTRIES/2;
  int bin_number = IMF_TABLE_ENTRIES/2;

  while (width > 1) {
    width /= 2;
    if (x > IMFData[bin_number]){
      bin_number += width;
    } else if (x < IMFData[bin_number]){
      bin_number -= width;
    } else{
      break;
    }
  } // found the bin

  m = IndividualStarIMFLowerMassCutoff * POW(10.0, bin_number * dm);

  IndividualStarIMFCalls++;

  if ( m > IndividualStarIMFUpperMassCutoff || m < IndividualStarIMFLowerMassCutoff){
    printf("individual_star_maker: IMF sampling not working m %"FSYM" %"FSYM"\n",m, IndividualStarIMFUpperMassCutoff);
  }

  return m;
}

int grid::individual_star_feedback(int *np,
                                   float *ParticleMass, int *ParticleType, FLOAT *ParticlePosition[],
                                   float *ParticleVelocity[], float *ParticleAttribute[]){
/*-----------------------------------------------------------------------------
  Handles the feedback for the indivual stars formed. This includes mechanical
  feedback from stellar winds, supernovae, and (if enabled) chemical yield 
  deposition.

  INPUTS
    nx, ny, nz   - size of grid in each dimension
    dx           - current grid size (code units)
    dt           - current timestep  (code units)
    current_time - time (code units)
    DensityUnits,LengthUnits,VelocityUnits,TimeUnits  - conversion between code units and cgs
    x/y/z start  - start position of grid in each dimension
    ibuff        - size of ghost zones in each dimension
    np           - number of particles to loop over
    ParticleMass -
    ParticlePosition -
    ParticleVelocity -
    ParticleAttribute -
-----------------------------------------------------------------------------*/
  /* copy some grid parameters for convenience */
  int  nx = this->GridDimension[0];
  int  ny = this->GridDimension[1];
  int  nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  FLOAT xstart = CellLeftEdge[0][0], ystart = CellLeftEdge[1][0], zstart = CellLeftEdge[2][0];
  float   dx   = CellWidth[0][0];



  /* Get Units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units

  float mp, distmass, energy, dratio;
  const double msolar = 1.989e33;                 // solar mass in cgs
  const double speed_of_light = 2.99792458e10 ;

  FLOAT particle_age, lifetime;

  const int max_random = (1<<16);

  int ip, jp, kp, index; // particle location

  int do_stellar_winds, go_supernova;

  const float fudge_factor = 1.0E10; // see comment below - hack to get WD and SNIa to work

  int IndividualStar        = PARTICLE_TYPE_INDIVIDUAL_STAR,
      IndividualStarWD      = PARTICLE_TYPE_INDIVIDUAL_STAR_WD,
      IndividualStarRemnant = PARTICLE_TYPE_INDIVIDUAL_STAR_REMNANT; // convenience


  // loop over all star particles
  for(int i = 0; i < (*np); i++){

    // where is the particle?
    ip = int ( (ParticlePosition[0][i] - (xstart)) / dx);
    jp = int ( (ParticlePosition[1][i] - (ystart)) / dx);
    kp = int ( (ParticlePosition[2][i] - (zstart)) / dx);

    float birth_mass = ParticleAttribute[3][i];

    mp = ParticleMass[i] * dx*dx*dx; // mass in code units
    mp = mp * (MassUnits) / msolar ;        // Msun

    // warning if outside current grid
    if( ip < 0 || ip > nx || jp < 0 || jp > ny || kp < 0 || kp > nz){
      printf("Warning: star particle is outside of grid\n");
      printf(" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" \n", ip, jp, kp, nx, ny, nz);
      printf(" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n", dx, xstart, ystart, zstart, ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i]);
      return FAIL;
    }


    particle_age = (this->Time) - ParticleAttribute[0][i];
    lifetime     = ParticleAttribute[1][i];

    do_stellar_winds           = FALSE;
    go_supernova               = FALSE;

    if(ParticleType[i] == IndividualStar){
      if(IndividualStarStellarWinds){

        float wind_start_age = 0.0;
        if(birth_mass < IndividualStarAGBThreshold){
          // star has AGB phase, do winds only in last fraction of lifetime
          if(IndividualStarInterpolateLifetime(wind_start_age, birth_mass,
                                               ParticleAttribute[2][i], 2) == FAIL){
            ENZO_FAIL("IndividualStarFeedback: Failure in main sequence lifetime interpolation");
          }
        }

        // dt addition here so stars w/ unresolved agb phases dump all winds at last timestep
        if (particle_age  < lifetime && particle_age + this->dtFixed > wind_start_age){
          do_stellar_winds = TRUE;
        }

      } // end winds check

      // move all attribute changes until AFTER feedback has occured!!!!!
      if ( birth_mass >= IndividualStarSNIIMassCutoff && ((particle_age + this->dtFixed) > lifetime)){
        go_supernova = TRUE;
        ParticleAttribute[1][i] = lifetime*1.0E10; // lifetime set to arbitrarily large number

      } else if ( birth_mass >= IndividualStarWDMinimumMass && birth_mass <= IndividualStarWDMaximumMass
                  && particle_age + (this->dtFixed) > lifetime && IndividualStarUseSNIa){

        float wd_mass; // compute WD mass using linear fit from Salaris et. al. 2009
        if(      birth_mass  < 4.0){ wd_mass = 0.134 * birth_mass + 0.331;}
        else if (birth_mass >= 4.0){ wd_mass = 0.047 * birth_mass + 0.679;}

        ParticleMass[i] = wd_mass * (msolar / MassUnits) / (dx*dx*dx);
        ParticleType[i] = IndividualStarWD;
        /* fudge factor makes lifetime very long so particle is not deleted,
           however, SNIa scheme needs to know the main sequence lifetime. This
           confusing, but a little bit more efficient than making a new particle attribute */
        ParticleAttribute[1][i] = lifetime*fudge_factor; // make a big number

      } else if (particle_age + (this->dtFixed) > lifetime){
        ParticleMass[i] = 0.0; // kill silently
      }


    } else if (ParticleType[i] == IndividualStarWD){
      /* White Dwarf Feedback - Make SNIa */

      /* Does the progenitor mass (main sequence mass) fit within range */
      if( (birth_mass > IndividualStarSNIaMinimumMass) &&
          (birth_mass < IndividualStarSNIaMaximumMass) ){

        float formation_time = ParticleAttribute[0][i]; // formation time of main sequence star

        float PSNIa;
        float rnum;

        /* Probability that the star will explode as SNIa in this timestep */
        PSNIa  = ComputeSnIaProbability( this->Time, formation_time, lifetime/fudge_factor, TimeUnits); // units of t^((beta)) / s
        PSNIa *= this->dtFixed;

        rnum =  (float) (random() % max_random) / ((float) (max_random));

//        printf("individual_star_feedback: SNIa - M_proj, PSNIa, rnum = %"ESYM" %"ESYM"\n", PSNIa, rnum);

        if (rnum < PSNIa){
          go_supernova = TRUE;
        }
      } // end SNIa progenitor + WD check

    } // end type check

    if(do_stellar_winds || go_supernova){
        float sum_dens = 0.0;

        // call feedback function to do stellar winds
        if(do_stellar_winds){


          printf("ISF: Calling feedback general to do stellar winds\n");
          printf("ISF: Current Mass = %"ESYM" Particle aatribute 3 = %"ESYM" mproj = %"ESYM"\n", mp,ParticleAttribute[3][i], birth_mass);

          /* Apply stellar wind feedback. Determined by setting last arguemtn to -1 */
          this->IndividualStarAddFeedbackGeneral(ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i],
                                                 ParticleVelocity[0][i], ParticleVelocity[1][i], ParticleVelocity[2][i],
                                                 birth_mass, ParticleAttribute[1][i], particle_age, ParticleAttribute[2][i], &mp, -1);

          ParticleMass[i] = mp* (msolar/MassUnits) / (dx*dx*dx); // update particle mass and put in code units
        }

        // call feedback function to do supernova feedback (either SNIa or core collapse)
        if(go_supernova){

          if( ParticleType[i] != IndividualStarWD){
            printf("Calling feedback to do cc supernova");
            /* do core collapse supernova feedback - set by last value == 1 */
            this->IndividualStarAddFeedbackGeneral(ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i],
                                                   ParticleVelocity[0][i], ParticleVelocity[1][i], ParticleVelocity[2][i],
                                                   birth_mass, ParticleAttribute[1][i], particle_age, ParticleAttribute[2][i], &mp, 1);

            ParticleMass[i] = mp * (msolar/MassUnits) / (dx*dx*dx); // update particle mass and put in code units
            ParticleType[i] = IndividualStarRemnant; // change type
            ParticleAttribute[1][i] = 1.0E10 * ParticleAttribute[1][i];
          } else{
            printf("calling feedback to do supernova 1a\n");
            /* do SNIa supernova feedback - set by last value == 1 */
            this->IndividualStarAddFeedbackGeneral(ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i],
                                             ParticleVelocity[0][i], ParticleVelocity[1][i], ParticleVelocity[2][i],
                                             birth_mass, ParticleAttribute[1][i], particle_age, ParticleAttribute[2][i], &mp, 2);

            ParticleMass[i]     = 0.0;                             // make particle mass zero - now a masless tracer
            ParticleAttribute[1][i] = 1.0E10 * ParticleAttribute[1][i];    // make lifetime infinite  -

          }






          // put attriute changes here
        }

        // if we went supernova, check if radius is resolved:
/*
        if(go_supernova){
          const double m_proton = 1.6726E-24;
          const double mu       = 1.31; // O.K. assumption since we are just doing this approx
          float n, r_rad;
          sum_dens /= (POW(IndividualStarFeedbackStencilSize,3)); // now average mass density
          n         = sum_dens*(*DensityUnits) / (mu * m_proton);     // in cgs
          r_rad = (7.32E19) * POW(energy * MassUnits * (*VelocityUnits)*(*VelocityUnits) / 1.0E51,0.29) / POW(n,0.42);

          if( (*dx) > 0.25* (r_rad/(*LengthUnits)) ){ // unresolved

            printf("IndividualStarFeedback: Unresolved supernova with <n> = %"ESYM" <E> = %"ESYM" r_rad = %"ESYM" and dx = %"ESYM"\n",n,energy*MassUnits*(*VelocityUnits)*(*VelocityUnits),r_rad/3.086E18,(*dx)*(*LengthUnits)/3.086E18);
          } else if( (*dx) > 0.125 * (r_rad/(*LengthUnits)) ){ // moderately resolved
            printf("IndividualStarFeedback: Moderately resolved SN with <n> = %"ESYM" <E> = %"ESYM" r_rad = %"ESYM" and dx = %"ESYM"\n",n,energy*MassUnits*(*VelocityUnits)*(*VelocityUnits),r_rad/3.086E18,(*dx) * (*LengthUnits) / 3.086E18);
          } else{ // yay
            printf("IndividualStarFeedback: Resolved supernova with <n> = %"ESYM" <E> = %"ESYM" r_rad = %"ESYM" and dx = %"ESYM"\n",n,energy*MassUnits*(*VelocityUnits)*(*VelocityUnits),r_rad/3.086E18, (*dx)*(*LengthUnits)/3.086E18);
          }

        }
*/
    } // if do feedback


  } // loop over particles

  return SUCCESS;
}

float ComputeSnIaProbability(const float &current_time, const float &formation_time, const float &lifetime, const float &TimeUnits){
 /* -----------------------------------------------------------------
  * ComputeSnIaProbability
  *
  * Computes dPdt for a given white dwarf that might go supernova. The
  * probability is such that the integral over dP from the time the WD
  * was born for a hubble time afterwards is a free parameter on order of
  * a few percent. This is IndividualStarSNIaFraction, or fraction of WD's
  * over a certain progenitor mass range that will go supernova in a hubble
  * time.
  * ------------------------------------------------------------------- */


 float dPdt;
 const float hubble_time = 4.382E17; // need to allow for cosmology units AJE TO DO
                                     // currently assumes H_o = 70.4

 dPdt = IndividualStarSNIaFraction;

 /* conmpute normalized probability - normalized by integral over WD formation time to hubble time */
 if (IndividualStarDTDSlope == 1.0){
   dPdt /= log( ((hubble_time / TimeUnits) + lifetime) / lifetime );
 } else{
   dPdt *= (-IndividualStarDTDSlope + 1.0);
   dPdt /= ( POW( (hubble_time / TimeUnits) + lifetime   , -IndividualStarDTDSlope + 1) -
             POW( (lifetime)                      , -IndividualStarDTDSlope + 1));
 }

 dPdt = dPdt * POW( ((current_time) - (formation_time)), -IndividualStarDTDSlope);


 return dPdt;
}

void ComputeStellarWindMassLossRate(const float &mproj, const float &metallicity,
                                   const float & lifetime, float *dMdt){
 /* ------------------------------------------------------------------
  * ComputeStellarWindEjectaMass
  * -------------------------------------------------------------------
  * A. Emerick - 4/22/16
  *
  * Model for stellar wind mass loss rate taken from Leitherer et. al. 1992.
  * This is the same model used in STARBURST 99 stellar wind models.
  * -------------------------------------------------------------------- */

  float L, Teff, Z, R;

  const double solar_z = 0.02; // as assumed in Leithener et. al. 1992
  const double yr      = 3.16224E7; // number of seconds in a year
  const double solar_L = 3.9E33;

 
  /* get properties */
  if(IndividualStarInterpolateLuminosity(L, mproj, metallicity) == FAIL){
    ENZO_FAIL("ComputeStellarWindMassLossRate: Failed to interpolate luminosity");
  }

  if(IndividualStarInterpolateProperties(Teff, R, mproj, metallicity) == FAIL){
    ENZO_FAIL("ComputeStellarWindMassLossRate: Failed to interpolate stellar properties");
  }

  /* compute logged mass loss rate */
  *dMdt = -24.06 + 2.45 * log10(L/solar_L) - 1.10 * log10(mproj) + 1.31 * log10(Teff)
                                   + 0.80 * log10(metallicity / solar_z);

  *dMdt = POW(10.0, *dMdt) / yr ; // Msun / yr -> Msun / s
}

void ComputeStellarWindVelocity(const float &mproj, const float &metallicity,
                                const float &lifetime, float *v_wind){
 /* ------------------------------------------------------------------
  * ComputeStellarWindVelocity
  * -------------------------------------------------------------------
  * A. Emerick - 4/22/16
  *
  * Model for stellar wind velocities taken from Leitherer et. al. 1992.
  * This is the same model used in STARBURST 99 stellar wind models.
  * The mass loss rate is computed elsewhere from stellar yields tables,
  * but velocity is set below using the fit function in luminosity,
  * stellar mass, effective temperature, and metallicity
  * -------------------------------------------------------------------- */

  float L, Teff, Z, R;

  const double solar_z = 0.02; // as assumed in Leithener et. al. 1992
  const double solar_L = 3.9E33;

  /* get properties */
  if (IndividualStarInterpolateLuminosity(L, mproj, metallicity) == FAIL){
    ENZO_FAIL("ComputeStellarWindVelocity: failure in interpolating luminosity");
  }
  if( IndividualStarInterpolateProperties(Teff, R, mproj, metallicity) == FAIL){
    ENZO_FAIL("ComputeStellarWindVelocity: failure in interpolating stellar properties");
  }

  // wind is in units of km / s
  // L - solar units
  // T - Kelvin
  // M - solar units
  // Z - solar units
  *v_wind = 1.23 - 0.30*log10(L/solar_L) + 0.55 * log10(mproj) + 0.64 * log10(Teff) + 0.13*log10(metallicity/solar_z);
  *v_wind = POW(10.0, *v_wind);
}


int grid::IndividualStarAddFeedbackGeneral(const FLOAT &xp, const FLOAT &yp, const FLOAT &zp,
                       const float &up, const float &vp, const float &wp,
                       const float &mproj, const float &lifetime, const float &particle_age,
                       const float &metallicity, float *mp, int mode     ){

  float m_eject, E_thermal, E_kin, f_kinetic, v_wind, p_feedback;
  const double c_light = 2.99792458E10; const double msolar = 1.989E33;

  float *metal_mass;

  /* Get Units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units


  /* rename some grid parameters for convenience */
  int  nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  FLOAT xstart = CellLeftEdge[0][0], ystart = CellLeftEdge[1][0], zstart = CellLeftEdge[2][0];
  float   dx   = CellWidth[0][0];


  printf("ISF: Assiging initital values to metal mass\n");
  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){

    metal_mass = new float[StellarYieldsNumberOfSpecies + 1];

    for (int i = 0; i < StellarYieldsNumberOfSpecies; i ++){
      metal_mass[i] = 0.0;
    }

  } else { metal_mass = NULL; }

  /* General function call to handle feedback precomputing of numbers */

  /* handle mass removal from particle here */

  float wind_dt = 0.0, wind_lifetime = lifetime;
  if (mode < 0){ // computing stellar winds

    /* -------------------------------------------
     *
     * Compute the mass ejecta and energetics for
     * stellar winds
     *
     * -------------------------------------------
     */

    /* compute ejecta mass - use yield tables if yields are tracked, otherwise use model */
    if ( IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
      m_eject   = StellarYieldsInterpolateYield(1, mproj, metallicity, -1) *msolar/ MassUnits; /* first arg, 1 = wind ; last -1 = tot mass */

      wind_lifetime = lifetime;
      if (mproj < IndividualStarAGBThreshold){

        float agb_start_time ; // point in star's life where AGB phase occurs

        if(IndividualStarInterpolateLifetime(agb_start_time, mproj, metallicity, 2) == FAIL){
          ENZO_FAIL("Individual star add feedback general: failure in interpolating lifetime");
        }

        agb_start_time /= TimeUnits;
        wind_lifetime   = lifetime - agb_start_time; // wind lifetime is particle life - start time of AGB

        //
        // To ensure total (integrated) mass ejected is accurate, make sure we don't overinject
        // mass when winds should only be "ON" for part of a timestep, either at beginning or end
        // of AGB phase, or when AGB phase is unresolved (i.e. AGB time < dt)
        //
        if (particle_age > agb_start_time && particle_age + this->dtFixed > lifetime) {
          wind_dt = lifetime - particle_age; // wind only occurs for part of timestep + star dies

        } else if (particle_age < agb_start_time && particle_age + this->dtFixed > lifetime) {
          //
          // AGB phase is unresolved. Set wind timestep to lifetime to do all ejecta this timestep
          //
          wind_dt = wind_lifetime; 

        } else if (particle_age < agb_start_time && particle_age + this->dtFixed > agb_start_time){
          wind_dt = particle_age + this->dtFixed - agb_start_time; // wind only occurs for part of timestep
        }

      } else { // massive stars (constant winds over lifetime)

        //
        // Check timestep to make sure we don't overinject yields at end of life
        //
        wind_dt = fmin(lifetime - particle_age, this->dtFixed);
      }

      m_eject  /= wind_lifetime ; // average mass loss rate over entire wind lifetime

    } else{
      // gives m_eject as Msun / s
      ComputeStellarWindMassLossRate(mproj, metallicity, lifetime * TimeUnits, &m_eject);
      m_eject = m_eject * msolar / MassUnits * TimeUnits;  // convert to code mass / code time

      // make sure we don't over-inject mass (i.e. partial timestep)
      wind_dt = fmin( lifetime - particle_age, this->dtFixed);
    }

    // Finally, compute total amount of mass ejected this timestep
    m_eject = m_eject * wind_dt;

    E_thermal = 0.0;

    /* set wind velocity depending on mode */
    if ( IndividualStarStellarWindVelocity < 0){
      ComputeStellarWindVelocity(mproj, metallicity, lifetime * (LengthUnits/VelocityUnits), &v_wind); // compute wind velocity in km / s using model
    } else {
      v_wind = IndividualStarStellarWindVelocity; // wind velocity in km / s
    }

    printf("ISF: Stellar wind mass in Msun = %"ESYM" in code units %"ESYM"\n", m_eject *MassUnits/msolar, m_eject);
    printf("ISF: Total Expected momentum in cgs %"ESYM" and in code units %"ESYM"\n", v_wind*1.0E5*m_eject*MassUnits/msolar, m_eject * v_wind*1.0E5/VelocityUnits);
    printf("ISF: Stellar wind in km / s %"ESYM" in code %"ESYM" and code vel = %"ESYM"\n", v_wind , v_wind *1.0E5/ VelocityUnits, VelocityUnits);

    v_wind     = (v_wind*1.0E5) / VelocityUnits; // convert from km/s to cm/s then to code units

    p_feedback = m_eject * v_wind;
    E_kin      = 0.0;

  } else if (mode == 1) {

    /* -------------------------------------------
     *
     * Compute the mass ejecta and energetics for
     * core collapse supernova
     *
     * -------------------------------------------
     */

    /* use yield tables to compute supernova ejecta mass - otherwise just eject the entire star */
    if ( IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
      m_eject   = StellarYieldsInterpolateYield(0, mproj, metallicity, -1) * msolar / MassUnits;
    } else{
      m_eject   = StarMassEjectionFraction * mproj * msolar / MassUnits;
    }

    if( IndividualStarSupernovaEnergy < 0){
      E_thermal = m_eject * StarEnergyToThermalFeedback * (c_light * c_light/(VelocityUnits*VelocityUnits));
    } else {
      E_thermal = IndividualStarSupernovaEnergy * 1.0E51 / (MassUnits*VelocityUnits*VelocityUnits);
    }

    v_wind     = 0.0;
    p_feedback = 0.0;
    E_kin      = 0.0;
    printf("AddFeedbackGeneral: M_proj %"FSYM" Z = %"FSYM", M_eject = %"ESYM" E_thermal = %"ESYM"\n", mproj, metallicity, m_eject*MassUnits/msolar, E_thermal*VelocityUnits*VelocityUnits*MassUnits);
  } else if ( mode == 2){ // Type Ia supernova

    /* -------------------------------------------
     *
     * Compute the mass ejecta and energetics for
     * Type Ia supernova
     *
     * -------------------------------------------
     */

    m_eject = StellarYields_SNIaYieldsByNumber(-1) * msolar / MassUnits; // total ejected mass

    if( IndividualStarSupernovaEnergy < 0){
      E_thermal = m_eject * StarEnergyToThermalFeedback * (c_light * c_light/(VelocityUnits*VelocityUnits));
    } else {
      E_thermal = IndividualStarSupernovaEnergy * 1.0E51 / (MassUnits*VelocityUnits*VelocityUnits);
    }

    v_wind     = 0.0;
    p_feedback = 0.0;
    E_kin      = 0.0;
  }


  /* if we are tracking yeilds, interpolat the ejecta mass for each species */
  printf("Tabulating metal ejecta mass fields \n");
  if(TestProblemData.MultiMetals == 2 && IndividualStarFollowStellarYields){
    int interpolation_mode;     // switch for stellar winds vs cc SN interpolation

    //
    // Switch around modes for interpolating either winds or supernova
    // to wrapper functions
    //
    if (mode == 1){
      interpolation_mode = 0;
    } else if (mode < 0){
      interpolation_mode = 1;
    }

    // for each metal species, compute the total metal mass ejected depending on supernova type
    if (mode == 0){

      metal_mass[0] = StellarYieldsInterpolateYield(1, mproj, metallicity, 0) *msolar / MassUnits / (dx*dx*dx);

      for (int i = 0; i < StellarYieldsNumberOfSpecies; i++){
        metal_mass[1 + i] = StellarYieldsInterpolateYield(1, mproj, metallicity, StellarYieldsAtomicNumbers[i]) * msolar / MassUnits / (dx*dx*dx);
      }

      // metal_mass now contains total mass ejected over wind lifetime. Adjust using wind loss rate and 
      // finite timestep check performed above
      for (int i = 0; i < StellarYieldsNumberOfSpecies + 1; i++){
        metal_mass[1+i] /= wind_dt / wind_lifetime;
      }

    } else if (mode == 1){
      // Core collapse supernova

      // First argument teslls interpolations to look at supernova yield tables
      // last argument (atomic number) of zero means get total metal mass
      metal_mass[0] = StellarYieldsInterpolateYield(0, mproj, metallicity, 0) * msolar / MassUnits / (dx*dx*dx);

      for(int i = 0; i < StellarYieldsNumberOfSpecies; i ++){
        metal_mass[1 + i] = StellarYieldsInterpolateYield(0, mproj, metallicity, StellarYieldsAtomicNumbers[i]) * msolar / MassUnits / (dx*dx*dx);
      }
    } else if (mode == 2){
      metal_mass[0] = StellarYields_SNIaYieldsByNumber(0) * msolar / MassUnits / (dx*dx*dx);

      for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
        metal_mass[1 + i] = StellarYields_SNIaYieldsByNumber( StellarYieldsAtomicNumbers[i] ) * msolar / MassUnits / (dx * dx * dx);
      }
    }


    printf("Metal masses in array ");
    for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
      printf(" %"ESYM " %"ESYM " -- ",metal_mass[i] * dx*dx*dx *MassUnits / msolar , metal_mass[i]);
    }
    printf("\n");
  }


  m_eject    = m_eject    / (dx*dx*dx);   // now in code units (code mass / code volume)
  p_feedback = p_feedback / (dx*dx*dx);
  E_thermal  = E_thermal  / (dx*dx*dx);
  E_kin      = E_kin      / (dx*dx*dx);

  /* find coordinates of feedback center. This is nominally the particle position
     but is shifted if needed if particle is too close to grid boundary.
     This is taken care of below (xfc = x feedback center) */
  float xfc, yfc, zfc;
  CheckFeedbackCellCenter( xp, yp, zp, xstart, ystart, zstart, dx,
                           nx, ny, nz, ibuff, &xfc, &yfc, &zfc);

  printf("ISF: injecting feedback to grid\n");
  this->IndividualStarInjectFeedbackToGrid(xfc, yfc, zfc,
                                           up, vp, wp,
                                           m_eject, E_thermal, E_kin, p_feedback, metal_mass); // function call

  // subtract mass from particle (in solar units)
  *mp = (*mp) - m_eject * (dx*dx*dx)*MassUnits/msolar;

  if( *mp < 0 && mode != 2){ // ignore this check for SNIa yields
      ENZO_FAIL("IndividualStarFeedback: Ejected mass greater than current particle mass - negative particle mass!!!\n");
  }

  delete [] metal_mass;

  return SUCCESS;
}

void CheckFeedbackCellCenter(const FLOAT &xp, const FLOAT &yp, const FLOAT &zp,
                             const FLOAT &xstart, const FLOAT &ystart, const FLOAT &zstart,
                             const FLOAT &dx,
                             const int &nx, const int &ny, const int &nz,
                             const int &ibuff,
                             FLOAT *xfc, FLOAT *yfc, FLOAT *zfc){

  /* checks if cell center is O.K. and rescales if needed */

  int fbuff = ibuff + ((int) (IndividualStarFeedbackStencilSize+1)/2.0 - 1); // number of cells away from edge
  float xfcshift, yfcshift, zfcshift;

  if ( xp < xstart +         fbuff*dx ||
       xp > xstart + dx*nx - fbuff*dx ||
       yp < ystart +         fbuff*dx ||
       yp > ystart + dx*ny - fbuff*dx ||
       zp < zstart +         fbuff*dx ||
       zp > zstart + dx*nz - fbuff*dx   ){

    *xfc      = fmin( fmax(xp, xstart + fbuff*dx), xstart + dx*nx - fbuff*dx);
    *yfc      = fmin( fmax(yp, ystart + fbuff*dx), ystart + dx*ny - fbuff*dx);
    *zfc      = fmin( fmax(zp, zstart + fbuff*dx), zstart + dx*nz - fbuff*dx);

    xfcshift = xp - *xfc;
    yfcshift = yp - *yfc;
    zfcshift = zp - *zfc;

    printf("Warning: Shifting feedback zone %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n", *xfc, *yfc, *zfc, xfcshift, yfcshift, zfcshift);
  } else{
    *xfc = xp;
    *yfc = yp;
    *zfc = zp;
  }
}


int grid::IndividualStarInjectFeedbackToGrid(const FLOAT &xfc, const FLOAT &yfc, const FLOAT &zfc,
                               float up, float wp, float vp,
                               float m_eject, float E_thermal, float E_kin, float p_feedback, float *metal_mass){

  float dx = float(CellWidth[0][0]);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;

  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  if (MultiSpecies == TRUE){
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }


  FLOAT xstart = CellLeftEdge[0][0], ystart = CellLeftEdge[1][0], zstart = CellLeftEdge[2][0];
  int nx = *(GridDimension), ny = *(GridDimension+1), nz = *(GridDimension+2);
  int number_of_cells;
  int stencil = IndividualStarFeedbackStencilSize; //renaming in case variable size
                                                   // is added in future
  /* for now, check total gas mass in injection region and print warning if too low */

  number_of_cells = POW(stencil, 3);
  /* scale everything to be the injected values in each cell */
  m_eject    = m_eject   / ((float) number_of_cells);
  E_thermal  = E_thermal / ((float) number_of_cells);
  p_feedback = p_feedback / ((float) number_of_cells - 1); // no momentum in center cell
  E_kin      = E_kin; // E_kin is totoal?   //  / ((float) number_of_cells);

  /* */

  float face_shift;
  float xface, yface, zface, dxf, dyf, dzf;
  float xpos, ypos, zpos, dxc, dyc, dzc;
  int iface, jface, kface, ic, jc, kc;

  /* check hydro method and shift things around
     to face-centered */
  face_shift = 0.0;
  if (HydroMethod == 2){ face_shift = 0.5; }

  /* AJE: I Suspect -0.5 should not be there in xface/xpos - June 2016*/

  xface = (xfc - xstart)/dx  - face_shift;
  yface = (yfc - ystart)/dx  - face_shift;
  zface = (zfc - zstart)/dx  - face_shift;

  iface = ((int) floor(xface + 0.5));
  jface = ((int) floor(yface + 0.5));
  kface = ((int) floor(zface + 0.5));

  dxf = iface + 0.5 - xface;
  dyf = jface + 0.5 - yface;
  dzf = kface + 0.5 - zface;

  /* we now need the index of the cell to add mass */
  xpos = (xfc - xstart)/dx;
  ypos = (yfc - ystart)/dx;
  zpos = (zfc - zstart)/dx;

  ic   = ((int) floor(xpos + 0.5));
  jc   = ((int) floor(ypos + 0.5));
  kc   = ((int) floor(zpos + 0.5));

  dxc  = ic + 0.5 - xpos;
  dyc  = jc + 0.5 - ypos;
  dzc  = kc + 0.5 - zpos;

  /* allocate local field stencil - bigger than feedback stencil */
  float *u_local, *v_local, *w_local, *d_local, *ge_local, *te_local;
  float *ke_before;

  int local_number_of_cells = POW(stencil + 1,3);

  u_local  = new float[local_number_of_cells];
  v_local  = new float[local_number_of_cells];
  w_local  = new float[local_number_of_cells];
  d_local  = new float[local_number_of_cells];
  ge_local = new float[local_number_of_cells];
  te_local = new float[local_number_of_cells];


  if (HydroMethod != 2){
    ke_before = new float[local_number_of_cells];
  } else {ke_before = NULL;}

  /* should go up to stencil + 1 I think (stencil = 3, 0, 1, 2, 3) */
  for(int k = 0; k < stencil + 1; k++){
    for(int j = 0; j < stencil + 1; j++){
      for(int i = 0; i < stencil + 1; i++){
        int index = i + (j + k * (stencil + 1)) * (stencil + 1);
        u_local[index] = 0.0;
        v_local[index] = 0.0;
        w_local[index] = 0.0;
        d_local[index] = 0.0;
        ge_local[index] = 0.0;
        te_local[index] = 0.0;

        if( HydroMethod != 2){
          ke_before[index] = 0.0;
        }
      }
    }
  }
  /* done allocating the zeroed arrays */

  /* add up kinetic energy in the clocal region */
  if( HydroMethod != 2 ){
    int integer_sep = ((int) (stencil+1)/2.0 - 1); //floor((stencil + 1) / 2.0);

    int local_index = 0; // AJE 5 - 10 - 16
    for(int k = kc - integer_sep; k <= kc + integer_sep + 1; k ++){
      for(int j = jc - integer_sep; j <= jc + integer_sep + 1; j++){
        for(int i = ic - integer_sep; i <= ic + integer_sep + 1; i++){
          /* check this index, I think it should be nx ny not stencil+1 */
          int index = i + (j + k * (ny))*(nx);

          ke_before[local_index] = 0.5 * BaryonField[DensNum][index] *
                                 ( (BaryonField[Vel1Num][index] * BaryonField[Vel1Num][index]) +
                                   (BaryonField[Vel2Num][index] * BaryonField[Vel2Num][index]) +
                                   (BaryonField[Vel3Num][index] * BaryonField[Vel3Num][index]));
          local_index++;
        }
      }
    }
  } // end kinetic energy sum for Zeus


  /* function call to convert velocities to momentum in particle frame */
  Momentum(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[DensNum],
           up, vp, wp, nx, ny, nz, ic, jc, kc, iface, jface, kface, stencil, 1);
  /* end convert to momenta */


  /* If needed, convert metals */
  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
    int field_num;
    this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0);
    MetalConversion( BaryonField[field_num], BaryonField[DensNum],
                     dx, nx, ny, nz, ic, jc, kc, stencil, 1);

    for(int m = 0; m < StellarYieldsSNData.NumberOfYields; m++){
      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[m]);

      MetalConversion( BaryonField[field_num], BaryonField[DensNum],
                       dx, nx, ny, nz, ic, jc, kc, stencil, 1);
    }
  }
  /* -------------------------------------------------- */

  /* compute the total mass and energy in cells before the explosion */
  float mass_before, energy_before, kin_energy_before;

  SumMassEnergy(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[DensNum],
                  BaryonField[GENum],   BaryonField[TENum],
                  nx, ny, nz, iface, jface, kface, ic, jc, kc, stencil,
                  &mass_before, &energy_before, &kin_energy_before);


  /* Now add mass and momentum terms to the dummy fields constructed earlier */
  AddFeedbackToGridCells(u_local, v_local, w_local, d_local, ge_local, te_local,
                         stencil+1, stencil+1, stencil+1, 1, 1, 1, 1, 1, 1,
                         dxf, dyf, dzf, dxc, dyc, dzc, m_eject, 1.0, 0.0, stencil);

  /* quadratic equation in p - momentum energy injection */
  float mom_per_cell = 0.0;
  if (E_kin > 0) {
    float A, B, C;
    ComputeAbcCoefficients( BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
                            BaryonField[DensNum], BaryonField[GENum]  ,
                            u_local, v_local, w_local, d_local,
                            nx, ny, nz, iface, jface, kface, ic, jc, kc, stencil,
                            A, B, C);

    A = A - (kin_energy_before + E_kin);
    mom_per_cell  = (-B + sqrt(B*B - 4.0 * A * C)) / (2.0 * C);

    printf("ISF: Coeffs mom_per_cell %"ESYM"\n", mom_per_cell);
  } else { // no kinetic energy injection - add feedback will add mass directly
    mom_per_cell = p_feedback;
  }

  /* add metal feedback - mass in cells */
  printf("ISF: Starting metal injection feedback calls\n");
  if(TestProblemData.MultiMetals == 2 && IndividualStarFollowStellarYields){
    /* For the first call, add in general metallicity field */
    int field_num;
    this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0); // when atomic number is zero, gives metallicity field

    AddMetalSpeciesToGridCells(BaryonField[field_num], metal_mass[0] / ((float) number_of_cells),
                               nx, ny, nz, ic, jc, kc, dxc, dyc, dzc, stencil);

    for(int ii = 0; ii < StellarYieldsNumberOfSpecies; ii++){
      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[ii]);

       printf("ISF: Adding metal feedaback for field %"ISYM" and atomic number %"ISYM" %"ISYM"\n", field_num, StellarYieldsAtomicNumbers[ii], ii);
       AddMetalSpeciesToGridCells(BaryonField[field_num], metal_mass[1 + ii] / ((float) number_of_cells),
                                  nx, ny, nz, ic, jc, kc, dxc, dyc, dzc, stencil);

    }

  } // 


  /* Now call add feedback again to add the feedback into the grid cells */
  AddFeedbackToGridCells(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
                         BaryonField[DensNum], BaryonField[GENum],   BaryonField[TENum]  ,
                         nx, ny, nz, ic, jc, kc, iface, jface, kface,
                         dxf, dyf, dzf, dxc, dyc, dzc, m_eject, mom_per_cell, E_thermal,
                         stencil);



  /* compute total mass and energy after feedback has been added */
  float mass_after, energy_after, kin_energy_after;
  SumMassEnergy(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
                    BaryonField[DensNum], BaryonField[GENum],   BaryonField[TENum],
                    nx, ny, nz, iface, jface, kface, ic, jc, kc, stencil,
                    &mass_after, &energy_after, &kin_energy_after);

  /* error checking statments go here */
  printf("ISF energy cheks, before %"ESYM" after %"ESYM" eject %"ESYM"\n", energy_before, energy_after, E_thermal);
  printf("ISF Mass checks, before %"ESYM" after %"ESYM" eject %"ESYM"\n", mass_before, mass_after, m_eject);
  printf("ISF Kinetic energy checks, before, after, E_kin %"ESYM" %"ESYM" %"ESYM"\n", kin_energy_before, kin_energy_after, E_kin);
  /*           -------------           */


  /* Now, reset the velocity fields to simulation frame */
  Momentum(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
           BaryonField[DensNum], up, vp, wp, nx, ny, nz, ic, jc, kc,
           iface, jface, kface, stencil, -1);

  /* If needed, convert metals back to metal fractions */
  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
    int field_num;
    this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0);
    MetalConversion( BaryonField[field_num], BaryonField[DensNum],
                     dx, nx, ny, nz, ic, jc, kc, stencil, -1);


    for(int m = 0; m < StellarYieldsSNData.NumberOfYields; m++){
      int field_num;
      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[m]);

      MetalConversion( BaryonField[field_num], BaryonField[DensNum],
                       dx, nx, ny, nz, ic, jc, kc, stencil, -1);
    }
  }
  /* -------------------------------------------------- */


  /* For sanity checking purposes */
  int number_of_bad_cells = 0;
  for(int ii = 0; ii < nx*ny*nz; ii++){
    float vx, vy, vz, vmag;
    float VelocityUnits = 3128941.28086;

    vx = BaryonField[Vel1Num][ii] * VelocityUnits;
    vy = BaryonField[Vel2Num][ii] * VelocityUnits;
    vz = BaryonField[Vel3Num][ii] * VelocityUnits;
    vmag = sqrt(vx*vx + vy*vy + vz*vz);

    float v_threshold = 1.0E4 * 1.0E5; // 10,000 km /s

    if (vx > v_threshold || vy > v_threshold || vz > v_threshold || vmag > v_threshold){
//      printf("Velocities are too large %"ESYM " %"ESYM" %"ESYM" %"ESYM"\n", vx, vy, vz, vmag);
      number_of_bad_cells++;
    }
  }
  printf("Velocities are too large cells %"ISYM" %"ISYM"\n", nx*ny*nz, number_of_bad_cells);


  /* Adjust the total energy field if we are using PPM */
  if (HydroMethod != 2){
    float ke_injected = 0.0;
    float delta_ke = 0.0;
    float ke_after = 0.0;
    int integer_sep = ((int) (stencil+1)/2.0 - 1); // floor((stencil + 1) / 2.0);
    printf("ISF kinetic feedback: integer_separation = %"ISYM"\n",integer_sep);

    int local_index = 0; // AJE 5 - 10 - 16
    for(int k = kc - integer_sep; k <= kc + integer_sep + 1; k ++){
      for(int j = jc - integer_sep; j <= jc + integer_sep + 1 ; j++){
        for(int i = ic - integer_sep; i <= ic + integer_sep + 1; i++){


//    for(int k = -integer_sep; k < integer_sep; k++){
//      for(int j = -integer_sep; j < integer_sep; j++){
//        for(int i = -integer_sep; i < integer_sep; i++){

          int index   = i + ( j + k*ny)*nx;

          ke_after = 0.5 * BaryonField[DensNum][index] *
                    ( BaryonField[Vel1Num][index] * BaryonField[Vel1Num][index] +
                      BaryonField[Vel2Num][index] * BaryonField[Vel2Num][index] +
                      BaryonField[Vel3Num][index] * BaryonField[Vel3Num][index]);

          delta_ke = ke_after - ke_before[local_index];

          BaryonField[TENum][index] += delta_ke/BaryonField[DensNum][index];

          ke_injected += delta_ke;

          local_index++;
        }
      }
    }
    printf("IndividualStarFeedback: change in kinetic energy %"ESYM"\n", ke_injected);
  } // endif


  /* free up memory */
  delete[] u_local;
  delete[] v_local;
  delete[] w_local;
  delete[] d_local;
  delete[] ge_local;
  delete[] te_local;
  delete[] ke_before;

  return SUCCESS;
}


void Momentum(float *u, float *v, float *w, float *d,
              const float &up, const float &vp, const float &wp,
              const int &nx, const int &ny, const int &nz,
              const int &ic, const int &jc, const int &kc,
              const int &iface, const int &jface, const int &kface,
              const int &stencil, int idir){

  int xo, yo, zo;
  // unit one shift in each direction
  xo = 1;
  yo = nx;
  zo = (nx * ny);

  /* making this a bigger region.... for stencil = 3, do -2 -1 0 1 2  -- 4/20/16 */
//  printf("Momentum Conversion: Direction %"ISYM"\n",idir);
//  int integer_sep = floor((stencil+1)/2.0) + 1; - AJE 5/10/16
  int integer_sep = ( (int) (stencil+1)/2.0 -1);
  for(int k = -integer_sep; k <= integer_sep + 1; k++){
    for(int j = -integer_sep; j <= integer_sep + 1 ; j++){
      for(int i = -integer_sep; i <= integer_sep + 1; i++){

        int index=0, x_index=0, y_index = 0, z_index = 0;
        if ( idir == 1.0) {
          if ( HydroMethod == 2 ){

            x_index = (iface + i) + ( (jc    + j) + (kc    + k) * ny) * nx;
            y_index = (ic    + i) + ( (jface + j) + (kc    + k) * ny) * nx;
            z_index = (ic    + i) + ( (jc    + j) + (kface + k) * ny) * nx;

  //          printf("Input u = %"ESYM" v = %"ESYM" w = %"ESYM"\n", u[x_index], v[y_index], w[z_index]);

            u[x_index] = (u[x_index] - up)*( 0.5 *(d[x_index] + d[x_index + xo]));
            v[y_index] = (v[y_index] - vp)*( 0.5 *(d[y_index] + d[y_index + yo]));
            w[z_index] = (w[z_index] - wp)*( 0.5 *(d[z_index] + d[z_index + zo]));

    //        printf("Output u = %"ESYM" v = %"ESYM" w = %"ESYM"\n", u[x_index], v[y_index], w[z_index]);

          } else {
            index = (ic + i) + ( (jc + j) + (kc + k)*ny)*nx;
            u[index] = (u[index] - up)*d[index];
            v[index] = (v[index] - vp)*d[index];
            w[index] = (w[index] - wp)*d[index];
          } // hydro method check


        } else { // reverse
          if( HydroMethod == 2){

            x_index = (iface + i) + ( (jc    + j) + (kc    + k) * ny) * nx;
            y_index = (ic    + i) + ( (jface + j) + (kc    + k) * ny) * nx;
            z_index = (ic    + i) + ( (jc    + j) + (kface + k) * ny) * nx;


      //      printf("Input u = %"ESYM" v = %"ESYM" w = %"ESYM"\n", u[x_index], v[y_index], w[z_index]);

            u[x_index] = u[x_index] / (0.5 * (d[x_index] + d[x_index + xo])) + up;
            v[y_index] = v[y_index] / (0.5 * (d[y_index] + d[y_index + yo])) + vp;
            w[z_index] = w[z_index] / (0.5 * (d[z_index] + d[z_index + zo])) + wp;

        //    printf("Output u = %"ESYM" v = %"ESYM" w = %"ESYM"\n", u[x_index], v[y_index], w[z_index]);
          } else{
            index = (ic + i) + ( (jc + j) + (kc + k) *ny)*nx;

            u[index] = u[index] / d[index] + up;
            v[index] = v[index] / d[index] + vp;
            w[index] = w[index] / d[index] + wp;
          }
        } // end idir check

        //printf("momentum_indexes: x y z %"ISYM" %"ISYM" %"ISYM"\n",x_index, y_index, z_index);
      }
    }
  } // end loop


// done
}

void MetalConversion(float *m, float *d, const float &dx,
                     const int &nx, const int &ny, const int &nz,
                     const int &ic, const int &jc, const int &kc,
                     const int &stencil, int idir){
 /* -----------------------------------------------------------------
  * MetalConversion
  * -----------------------------------------------------------------
  * A. Emerick - 04/19/16
  * ----------------------
  * Converts an arbitrary metal field from metal density to metal
  * density (a proxy for mass since we are operating on a grid with
  * uniform dx) and vice versa depending on idir. Feedback injection
  * occurs in terms of metall density (mass), field stored as fraction
  * ------------------------------------------------------------------- */

  /* any add all changes to momentum above should (probably) be reflected here as well */
  /* Metal fields give metal density in a given cell - convert to mass */
//  int integer_sep = floor((stencil+1)/2.0) + 1;
  int integer_sep = ((int) (stencil+1)/2.0 - 1);
  for(int k = -integer_sep; k <= integer_sep + 1; k++){
    for(int j = -integer_sep; j <= integer_sep + 1; j++){
      for(int i = -integer_sep; i <= integer_sep + 1; i++){

        int index = (ic + i) + ( (jc + j) + (kc + k) * ny) * nx;

        if (idir == 1.0){
          m[index] = m[index] * d[index];
        } else {
          m[index] = m[index] / d[index];
        }

      } // end i
    } // end j
  } // end k
} // done with metal conversion

void SumMassEnergy(float *pu, float *pv, float *pw, float *d, float *ge, float *te,
                   const int &nx, const int &ny, const int &nz, 
                   const int &iface, const int &jface, const int &kface,
                   const int &ic, const int &jc, const int &kc, const int &stencil,
                   float *mass_sum, float *energy_sum, float *kin_energy_sum){

  *mass_sum = 0.0; *energy_sum = 0.0; *kin_energy_sum = 0.0;

  int xo, yo, zo;
  // unit one shift in each direction
  xo = 1;
  yo = nx;
  zo = (nx * ny);

//  int integer_sep = floor((stencil+1)/2.0);
  int integer_sep = ((int) (stencil+1)/2.0 - 1);
  for(int k = -integer_sep; k <= integer_sep+1; k++){
    for(int j = -integer_sep; j <= integer_sep+1; j++){
      for(int i = -integer_sep; i <= integer_sep+1; i++){
        float mass_term, mom_term, gas_energy = 0.0, kinetic_energy;

        int index   = (ic + i) + ( (jc + j) + (kc + k)*ny)*nx;
        int x_index = (iface + i) + ( (jc    + j) + (kc    + k)*ny)*nx;
        int y_index = (ic    + i) + ( (jface + j) + (kc    + k)*ny)*nx;
        int z_index = (ic    + i) + ( (jc    + j) + (kface + k)*ny)*nx;

        mass_term = d[index];
        mom_term  = pu[x_index]*pu[x_index] +
                    pv[y_index]*pv[y_index] +
                    pw[z_index]*pw[z_index] ;

        /* total mass and energy */
        kinetic_energy   = mom_term / (2.0 * mass_term);
        *mass_sum        = *mass_sum + mass_term;
        *kin_energy_sum  = *kin_energy_sum + kinetic_energy;

        /* account for thermal energy depending on Hydro solver */
        if (HydroMethod == 2){
          gas_energy = te[index] * d[index];
        }
        if (DualEnergyFormalism) {
          gas_energy = ge[index]*d[index];
        }
        if (HydroMethod != 2 && DualEnergyFormalism == 0) {
          gas_energy = te[index]*d[index] - kinetic_energy;
        }

        *energy_sum = *energy_sum + kinetic_energy + gas_energy;

      }
    }
  }// end loop


// done with routine
}


void ComputeAbcCoefficients(float *pu, float *pv, float *pw, float *d,
                            float *ge, float *pu_l, float *pv_l, float *pw_l,
                            float *d_l, const int &nx, const int &ny, const int &nz,
                            const int &iface, const int &jface, const int &kface,
                            const int &ic, const int &jc, const int &kc,
                            const int &stencil, float &A, float &B, float &C){

  float mass_sum, energy_sum;

  A = 0.0; B = 0.0; C = 0.0; mass_sum = 0.0; energy_sum = 0.0;

  int xo, yo, zo;
  // unit one shift in each direction
  xo = 1;
  yo = nx;
  zo = (nx * ny);

  float mass_term=0.0, mom_term=0.0, b_term=0.0, c_term = 0.0;

  int integer_sep = ((int) (stencil+1)/2.0 - 1);
  int l_index     = 0; int index = 0;

  for(int k = -integer_sep; k <= integer_sep + 1; k++){
    for(int j = -integer_sep; j <= integer_sep + 1; j++){
      for(int i = -integer_sep; i <= integer_sep + 1; i++){
        int x_index, y_index, z_index;

        /* this may be memory issue -- check this */

        index   = (ic + i) + ( (jc + j) + (kc + k)*ny)*nx;
        x_index = (iface + i) + ( (jc    + j) + (kc    + k)*ny)*nx;
        y_index = (ic    + i) + ( (jface + j) + (kc    + k)*ny)*nx;
        z_index = (ic    + i) + ( (jc    + j) + (kface + k)*ny)*nx;


        mass_term = d[index];
        mom_term  = pu[x_index]*pu[x_index] +
                    pv[y_index]*pv[y_index] +
                    pw[z_index]*pw[z_index] ;


        mass_term += d_l[l_index];

        b_term     = pu[x_index]*pu_l[l_index] +
                     pv[y_index]*pv_l[l_index] +
                     pw[z_index]*pw_l[l_index]  ;

        c_term     = pu_l[l_index] * pu_l[l_index] +
                     pv_l[l_index] * pv_l[l_index] +
                     pw_l[l_index] * pw_l[l_index] ;


        A         += mom_term/(2.0 * mass_term);
        B         += b_term/mass_term;
        C         += c_term/(2.0 * mass_term);

        l_index++;
      }
    }
  }
  printf("ComputeAbcCoefficients: local_index = %"ISYM" integer_sep = %"ISYM" A = %"ESYM" B = %"ESYM" C = %"ESYM"\n", l_index, index, A, B, C);


} // end comput coeff

void AddMetalSpeciesToGridCells(float *m, const float &mass_per_cell,
                                const int &nx, const int &ny, const int &nz,
                                const int &ic, const int &jc, const int &kc,
                                const float &dxc, const float &dyc, const float &dzc,
                                const int &stencil){
 /* -------------------------------------------------------------------------
  * AddMetalSpeciesToGridCells
  * --------------------------------------------------------------------------
  * Adds in metal species deposition for a given metal field and mass ejection.
  * This is a copy / adaptation of the below algorithm (AddFeedbackToGridCells).
  * Any modification to that function should be reflected here as well.
  * I know this is gross, but it is somewhat more efficient.
  * ------------------------------------------------------------------------- */

  int integer_sep = ((int) (stencil+1)/2.0 - 1);
  float delta_mass = 0.0, total_mass = 0.0;

  for (int k = -integer_sep; k <= integer_sep; k++){
    for (int j = -integer_sep; j <= integer_sep; j++){
      for (int i = -integer_sep; i <= integer_sep; i++){

        /*  */
        for(int i_loc = i; i_loc <= i + 1; i_loc++){
          float dxc_loc = ( (i_loc == i) ? dxc : 1.0 - dxc);

          for(int j_loc = j; j_loc <= j + 1; j_loc++){
            float dyc_loc = ( (j_loc == j) ? dyc : 1.0 - dyc);

            for( int k_loc = k; k_loc <= k + 1; k_loc++){
              float dzc_loc = ( (k_loc == k) ? dzc : 1.0 - dzc);


              int index = (ic + i_loc) + ( (jc + j_loc) + (kc + k_loc)*ny) * nx;

              delta_mass = mass_per_cell * dxc_loc * dyc_loc * dzc_loc;

              m[index] = m[index] + delta_mass;

              total_mass += delta_mass;
            } //

          } //

        } //

      }
    }
  } // end k loop 

  printf("MetalFeedback: Deposited total metal mass (density) %"ESYM"\n", total_mass);
}

void AddFeedbackToGridCells(float *pu, float *pv, float *pw,
                            float *d, float * ge, float *te,
                            const int &nx, const int &ny, const int &nz,
                            const int &ic, const int &jc, const int &kc,
                            const int &iface, const int &jface, const int &kface,
                            const float &dxf, const float &dyf, const float &dzf,
                            const float &dxc, const float &dyc, const float &dzc,
                            const float &mass_per_cell, const float &mom_per_cell,
                            const float &therm_per_cell, const int &stencil){


  int integer_sep = ((int) (stencil+1)/2.0 - 1);
  float total_mass = 0.0, delta_therm =0.0;

  // should go over a stencilxstencilxstencil region (if stencil is 3, -1, 0, 1)
  for(int k = -integer_sep; k <= integer_sep; k++){
    for(int j = -integer_sep; j <= integer_sep; j++){
      for(int i = -integer_sep; i <= integer_sep; i++){

        for(int i_loc = i; i_loc <= i + 1; i_loc++){
          float dxf_loc, dxc_loc;

          dxf_loc = dxf;  dxc_loc = dxc;
          if( i_loc == i + 1){
            dxf_loc = 1.0 - dxf;
            dxc_loc = 1.0 - dxc;
          }

          for(int j_loc = j; j_loc <= j + 1; j_loc++){
            float dyf_loc, dyc_loc;

            dyf_loc = dyf;  dyc_loc = dyc;
            if( j_loc == j + 1){
              dyf_loc = 1.0 - dyf;
              dyc_loc = 1.0 - dyc;
            }

            for(int k_loc = k; k_loc <= k + 1; k_loc++){
              float dzf_loc, dzc_loc;

              dzf_loc = dzf;    dzc_loc = dzc;
              if( k_loc == k + 1){
                dzf_loc = 1.0 - dzf;
                dzc_loc = 1.0 - dzc;
              }
              /* now we do the feedback */
              int index, x_index, y_index, z_index;
              float delta_mass, delta_pu, delta_pv, delta_pw, inv_dens;

              index = (ic + i_loc) + ( (jc + j_loc) + (kc + k_loc)*ny)*nx;

              x_index = (iface + i_loc) + ( (jc    + j_loc) + (kc    + k_loc)*ny)*nx;
              y_index = (ic    + i_loc) + ( (jface + j_loc) + (kc    + k_loc)*ny)*nx;
              z_index = (ic    + i_loc) + ( (jc    + j_loc) + (kface + k_loc)*ny)*nx;

              /* do things here finally */
              delta_mass  = mass_per_cell * dxc_loc * dyc_loc * dzc_loc;

              // add momentum, then do sign changing later
              delta_pu    = mom_per_cell * dxf_loc * dyc_loc * dzc_loc;
              delta_pv    = mom_per_cell * dxc_loc * dyf_loc * dzc_loc;
              delta_pw    = mom_per_cell * dxc_loc * dyc_loc * dzf_loc;

              // change sign to point away from central sell
              // +1 or -1 if index is pos or neg respectively
              // multiply by zero if we are the central cell
              delta_pu *=  ( (i > 0) ? 1.0 : ( i < 0 ? -1.0 : 0));
              delta_pv *=  ( (j > 0) ? 1.0 : ( j < 0 ? -1.0 : 0));
              delta_pw *=  ( (k > 0) ? 1.0 : ( k < 0 ? -1.0 : 0));

              delta_therm = therm_per_cell * dxc_loc * dyc_loc * dzc_loc;

              /* add mass momentum and thermal energy */
              inv_dens   = 1.0 / (d[index] + delta_mass);


              pu[x_index] = pu[x_index] + delta_pu;
              pv[y_index] = pv[y_index] + delta_pv;
              pw[z_index] = pw[z_index] + delta_pw;

              total_mass  += delta_mass;

              te[index] = (te[index]*d[index] + delta_therm) * inv_dens;

              if (DualEnergyFormalism) {
                ge[index] = (ge[index]*d[index] + delta_therm) * inv_dens;
              }
              d[index] = d[index] + delta_mass;

              if( d[index] < 0.0 || delta_mass < 0.0){
                printf("individual_star: Feedback producing negative densities %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n", d[index], delta_mass, te[index], delta_therm);
                printf("--------------------------------- %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",mass_per_cell,dxf_loc, dxc_loc, dyf_loc, dyc_loc, dzf_loc, dzc_loc);
              }
          //    printf("feedback_indexes: i x y z %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",index, x_index, y_index, z_index);

            } // k loc
          } // j loc
        } // i loc


      }
    }
  }//end loop
  printf("AddedFeedbackToGridCells: mom_per_cell = %"ESYM" therm_per_cell = %"ESYM"\n", mom_per_cell, therm_per_cell);
  printf("AddedFeedbackToGridCells: M_tot =  %"ESYM" therm = %"ESYM"\n", total_mass, delta_therm);

}
