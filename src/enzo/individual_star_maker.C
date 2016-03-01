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

/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);
float SampleIMF(void);
float GaussianRandomVariable(void);
float compute_lifetime(float *mp);
unsigned_long_int mt_random(void);


int grid::individual_star_maker(int *nx, int *ny, int *nz, int *size,
                                float *dm, float *temp, float *dt,
                                float *dx, FLOAT *t, int *procnum,
                                float *d1, float *x1, float *v1, float *t1,
                                int *nmax, FLOAT *xstart, FLOAT *ystart,
                                FLOAT *zstart, int *ibuff, int *imethod,
                                float *mu, float *metal, int *ctype,
                                int *np, float *ParticleMass,
                                int *ParticleType, FLOAT *ParticlePosition[],
                                float *ParticleVelocity[], float *ParticleAttribute[]){

  const double msolar = 1.989e33;
  const double sndspdC = 1.3095e8;
  double m1   = (*d1)*POW(*x1,3); // mass units

  int i, j, k, index, ii=0, istar=0, index_presf=0;
  int xo, yo, zo, rsign=1;
  float bmass, div, min_temp, star_mass=0.0, sum_mass=0.0;
  float pstar, mass_to_stars, mass_available, tdyn;
  float dtot, isosndsp2, jeansmass, star_fraction, odthreshold;
  float umean, vmean, wmean, px, py, pz, px_excess, py_excess, pz_excess;
  float rnum;
  const double m_h = 1.673e-24;

  const int max_random = (1<<16);

  int form_method = -1; // tracker for debugging purposes


  if ((*dt) == 0.0){
    printf("DT EQUAL TO ZERO\N");
    return FAIL;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;

  int CINum, NINum, OINum, MgINum, SiINum, FeINum, YINum, BaINum, LaINum, EuINum;

  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  if(TestProblemData.MultiMetals >= 2){
    if(IdentifyChemicalTracerSpeciesFields(CINum, NINum, OINum, MgINum, SiINum,
                                           FeINum, YINum, BaINum, LaINum, EuINum) == FAIL){
      ENZO_FAIL("Failure in identifying chemical tracer species fields.");
    }
  }

  // Particle attributes hard coded for chemical tagging numbers. This is NOT IDEAL
  // AJE : TO DO is to address this issue with a lookup method like the above for baryon fields
  //       this may necessitate adding a new ParticleAttribute like array to stars specifically
  //       for chemical tagging

    // 3D -> 1D index
    xo = 1;
    yo = *nx;
    zo = (*nx) * (*ny);

    min_temp = 1.0E5; // set conditional based on cooling and metals present

    // over density threshold in code units
    odthreshold = StarMakerOverDensityThreshold * m_h * (*mu) / (*d1);

    // loop over all cells, check condition, form stars stochastically
    ii = 0; index_presf = 0;
    for (k = *ibuff; k < *nz - *ibuff; k++){
      for (j = *ibuff; j < *ny - *ibuff; j++){
        for (i = *ibuff; i < *nx - *ibuff; i++){
          index = i + (j + k * (*ny)) * (*nx);

          bmass = (BaryonField[DensNum][index]*(*dx)*(*dx)*(*dx)) * m1 / msolar; // in solar masses

          // perform the following easy checks for SF before proceeding
          // 1) Is density greater than the density threshold?
          // 2) Is temperature < the minimum temperature?
          // 3) Do not allow star formation if minimum star particle mass is
          //    above some fraction of cell mass. This is very unlikely to occur
          //    in intended use case:
          //    (i.e. baryon mass likely always going to be > ~10 solar masses)

          sum_mass = 0.0; index_presf = ii;
          if (   BaryonField[DensNum][index]      > odthreshold
              && temp[index] <= min_temp
              && IndividualStarMassFraction*bmass > IndividualStarIMFLowerMassCutoff
              && 0.5*bmass > IndividualStarIMFUpperMassCutoff){


            // star formation may be possible
            // compute values and check jeans mass unstable
            dtot = ( BaryonField[DensNum][index] + dm[index] ) * (*d1);         // total density
            tdyn = sqrt(3.0 * pi / 32.0 / GravConst / dtot) / (*t1);            // in code units
            isosndsp2 = sndspdC * temp[index] ;
            jeansmass = pi / (6.0 * sqrt(BaryonField[DensNum][index]*(*d1)) *
                            POW(pi * isosndsp2 / GravConst ,1.5)) / msolar; // in solar masses

            if (jeansmass <= bmass){

              // calculate mass in cell that can be converted to stars in timestep
              // generally this should be small (comparable to or less than the lower mass
              // cutoff of the IMF)

              star_fraction  = min(StarMakerMassEfficiency*(*dt)/tdyn, 1.0);
              mass_to_stars  = star_fraction * bmass;
              mass_available = StarMakerMassEfficiency * bmass;
              mass_to_stars  = min(mass_to_stars, mass_available);

              // If mass_to_stars greater than available mass, convert
              // all of available mass into stars
              // Frankly this is very unlikely to occur...
              // Tests as of 2/22/16 show NO SF here for at least 10^5 stars in a LMC dwarf galaxy
              if(mass_to_stars >= mass_available){
                mass_to_stars = mass_available;
                while( ii < *nmax && mass_to_stars > IndividualStarIMFUpperMassCutoff){
                  ParticleMass[ii] = SampleIMF();
                  sum_mass        += ParticleMass[ii]; // counter for mass formed in this cell
                  mass_to_stars   -= ParticleMass[ii]; // reduce available mass
                  ii++;
                }
              }

              // Tests (as of 2/22/16) show NO SF here for at least the first 10^5 stars
              if (mass_to_stars > IndividualStarIMFUpperMassCutoff){
                while (ii < *nmax && mass_to_stars > IndividualStarIMFUpperMassCutoff){
                  ParticleMass[ii]  = SampleIMF();
                  sum_mass         += ParticleMass[ii];
                  mass_to_stars    -= ParticleMass[ii];
                  ii++;
                }
              }

              // If mass is above IMF lower limit, star formation will happen.
              // Just form stars randomly over IMF until mass dips below lower cutoff
              if(mass_to_stars > IndividualStarIMFLowerMassCutoff){
                while( ii < *nmax && mass_to_stars > IndividualStarIMFLowerMassCutoff){
                  ParticleMass[ii]  = SampleIMF();
                  sum_mass         += ParticleMass[ii];
                  mass_to_stars    -= ParticleMass[ii];
                  ii++;

                  if (mass_to_stars < 0.0){
                    mass_to_stars = 0.0;
                  }
                }
              }

              // now we are in the Goldbaum et. al. 2015 regime (star_maker_ssn.F)
              // Calculate probability of star forming and form stars stochastically
              if (mass_to_stars < IndividualStarIMFLowerMassCutoff && mass_to_stars > tiny_number){
                star_mass = SampleIMF();
                pstar     = mass_to_stars / star_mass;
                rnum =  (float) (random() % max_random) / (float) (max_random);
                if (rnum < pstar){
                  ParticleMass[ii]  = star_mass;
                  sum_mass         += ParticleMass[ii];
                  ii++;
                }
              }

              // prepare for assigning star properties by computing the local
              // gas velocity properties (this is for velocity assignment)
              // 2 = Zeus .. otherwise PPM
              // copied from pop3_maker.F
              if (*imethod == 2){
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

                ParticleType[istar]            = -(*ctype);   // negative is a "new" star
                ParticleAttribute[0][istar]    = *t;          // formation tim
                ParticleAttribute[1][istar]    = compute_lifetime( &ParticleMass[istar] ) / (*t1); // lifetime
                ParticleMass[istar]            = ParticleMass[istar] * msolar / m1;                // mass

                // give the star particle a position chosen at random
                // within the cell ( so they are not all at cell center )

                rnum =  (float) (random() % max_random) / (float) (max_random);
                ParticlePosition[0][istar] = (*dx) * rnum + *xstart + ((float) i + 0.5)*(*dx);
                rnum =  (float) (random() % max_random) / (float) (max_random);
                ParticlePosition[1][istar] = (*dx) * rnum + *ystart + ((float) j + 0.5)*(*dx);
                rnum =  (float) (random() % max_random) / (float) (max_random);
                ParticlePosition[2][istar] = (*dx) * rnum + *zstart + ((float) k + 0.5)*(*dx);

                // assume velocity dispersion is isotropic in each velocity component. Multiply disp by
                // sqrt(1/3) to get disp in each component... taking above velocities as the mean
                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[0][istar] = umean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0E5*(*t1)/(*x1);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[1][istar] = vmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0e5*(*t1)/(*x1);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[2][istar] = wmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0E5*(*t1)/(*x1);

                // ENSURE MOMENTUM CONSERVATION!!!!!
                // make running total of momentum in each direction
                px += ParticleVelocity[0][istar]*ParticleMass[istar];
                py += ParticleVelocity[1][istar]*ParticleMass[istar];
                pz += ParticleVelocity[2][istar]*ParticleMass[istar];

                // this is where code would go to assign
                // chemical tags to all of the particles
                // depending on whether or not multimetals is ON
                if(TestProblemData.MultiMetals >= 2){
                  if(MULTIMETALS_METHOD(MULTIMETALS_ALPHA)){
                    ParticleAttribute[ 4][istar] = BaryonField[ CINum][index];
                    ParticleAttribute[ 5][istar] = BaryonField[ NINum][index];
                    ParticleAttribute[ 6][istar] = BaryonField[ OINum][index];
                    ParticleAttribute[ 7][istar] = BaryonField[MgINum][index];
                    ParticleAttribute[ 8][istar] = BaryonField[SiINum][index];
                    ParticleAttribute[ 9][istar] = BaryonField[FeINum][index];
                  }
                  if(MULTIMETALS_METHOD(MULTIMETALS_SPROCESS)){
                    ParticleAttribute[10][istar] = BaryonField[ YINum][index];
                    ParticleAttribute[11][istar] = BaryonField[BaINum][index];
                    ParticleAttribute[12][istar] = BaryonField[LaINum][index];
                  }
                  if(MULTIMETALS_METHOD(MULTIMETALS_RPROCESS)){
                    ParticleAttribute[13][istar] = BaryonField[EuINum][index];
                  }
                } // end multi metals

              } // end while loop for assigning particle properties
              // ---------------------------------------------------

              // ensure zero net momentum from mean velocity
              // momentum of gas converted into stars (sum_mass * umean)
              // should be equal to the total momentum of the stars
              // compute excess momentum and modify star velocity evenly (mass weighted)
              // this is not completely physical, as pre-SF and post-SF gas vel is the same
              sum_mass = sum_mass * msolar / m1; // in code units


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

              // now remove mass from grid
              BaryonField[DensNum][index] = (bmass*msolar/m1 - sum_mass) / ((*dx)*(*dx)*(*dx)) ;

            } // if jeans mass unstable
          } // resolution and density


        } // enx x loop
      } // end y loop
    } // end z loop

    // Done forming stars!!! Output and exit
    if (ii > 0){
      printf("P(%"ISYM"): individual_star_maker[add]: %"ISYM" new star particles\n", *procnum, ii);
    }
    if (ii >= *nmax){
      fprintf(stdout, "individual_star_maker: reached max new particle count!! Available: %"ISYM". Made: %"ISYM"\n", *nmax, ii);

    }

  // star masses are recorded as densities (mass / cell volume)
  for (int counter = 0; counter < ii; counter++){
    ParticleMass[counter] = ParticleMass[counter] / ((*dx)*(*dx)*(*dx)); // code units / cell volume
  }

  *np = ii - 1; // number of stars formed : AJE 2/29 check if this is a bug with the -1

  return SUCCESS;
}


float SampleIMF(void){
  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  float x = (float) (random_int%max_random) / (float) (max_random);
  float dm = log10(IndividualStarIMFUpperMassCutoff / IndividualStarIMFLowerMassCutoff)/ (float) (IMF_TABLE_ENTRIES-1);
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

  if (IndividualStarIMFUpperMassCutoff != 100.0 || IndividualStarIMFLowerMassCutoff != 1.0){
    printf("individual_star_maker: Bounds are wrong!!!!\n");
  }

  return m;
}

float GaussianRandomVariable(void){
 // returns gaussian random variable y1

    const int max_random = (1<<16);

    float y1, y2, w;
    float x1 = (float) (random() % max_random) / ((float) (max_random));
    float x2 = (float) (random() % max_random) / ((float) (max_random));

    do {
       x1 = 2.0 * (float) (random() % max_random) / ((float) (max_random)) - 1.0;
       x2 = 2.0 * (float) (random() % max_random) / ((float) (max_random)) - 1.0;
       w  = x1 * x1 + x2 * x2;
    } while ( w >= 1.0);


  w = sqrt( ( -2.0 * log( w ) ) / w );
  y1 = x1 * w;
  y2 = x2 * w;

  return y1;

}

int grid::individual_star_feedback(int *nx, int *ny, int *nz,
                                   float *dx,
                                   FLOAT *current_time, float *d1, float *x1,
                                   float *v1, float *t1, FLOAT *xstart, FLOAT *ystart,
                                   FLOAT *zstart, int *ibuff, int *np,
                                   float *ParticleMass, FLOAT *ParticlePosition[],
                                   float *ParticleVelocity[], float *ParticleAttribute[]){

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  int CINum, NINum, OINum, MgINum, SiINum, FeINum, YINum, BaINum, LaINum, EuINum;


  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  if (CRModel){
    if ((CRNum = FindField(CRDensity, FieldType, NumberOfBaryonFields)) < 0)
      ENZO_FAIL("Cannot Find Cosmic Rays");
  }

  if(TestProblemData.MultiMetals >=2){

    if(IdentifyChemicalTracerSpeciesFields(CINum, NINum, OINum, MgINum, SiINum,
                                           FeINum, YINum, BaINum, LaINum, EuINum) == FAIL){
      ENZO_FAIL("Error in IdentifyChemicalTracerSpeciesFields.");
    }

  }


  float mp, particle_age, lifetime, distmass, energy, dratio;
  const double msolar = 1.989e33;                 // solar mass in cgs
  const double speed_of_light = 2.99792458e10 ;
  double m1 = (*d1)*(*x1)*(*x1)*(*x1);            // code mass units
  double e1 = m1 * (*x1) * (*x1) / (*t1) * (*t1); // code energy units

  int ip, jp, kp, index; // particle location

  printf("IF Feedback about to loop over particles %"ISYM"\n",(*np));

  // loop over all star particles
  for(int i = 0; i < (*np); i++){

    // where is the particle?
    ip = int ( (ParticlePosition[0][i] - (*xstart)) / (*dx));
    jp = int ( (ParticlePosition[1][i] - (*ystart)) / (*dx));
    kp = int ( (ParticlePosition[2][i] - (*zstart)) / (*dx));

    mp = ParticleMass[i] * (*dx) * (*dx) * (*dx); // mass in code units
    mp = mp * (m1) / msolar ; // Msun

    // warning if outside current grid
    if( ip < 0 || ip > (*nx) || jp < 0 || jp > (*ny) || kp < 0 || kp > (*nz)){
      printf("Warning: star particle is outside of grid\n");
      return FAIL;
    }

    // Modify position if too close to grid boundary
    // If too close, pretends particles are at a new position, which
    // is the to closest they can be to grid without causing
    // feedback require MPI calls. Does not actually move particle
    if (StarFeedbackDistRadius > 0){
      ip = fmax( (*ibuff) + StarFeedbackDistRadius,

               fmin( (*nx) - (*ibuff) - StarFeedbackDistRadius, ip) );
      jp = fmax( (*ibuff) + StarFeedbackDistRadius,
               fmin( (*ny) - (*ibuff) - StarFeedbackDistRadius, jp) );
      kp = fmax( (*ibuff) + StarFeedbackDistRadius,
               fmin( (*nz) - (*ibuff) - StarFeedbackDistRadius, kp) );
    }

    particle_age = abs( (*current_time) - ParticleAttribute[0][i] );
    lifetime     = compute_lifetime( &mp ) / (*t1);

    // if true, do stellar winds
//    if(abs(*current_time / lifetime) < 1.0){
      // do stellar winds and radiation here
//    }

    if( particle_age > lifetime ){

      if(mp < IndividualStarTypeIIMassCutoff){

        ParticleMass[i] = 0.0; // delete star and do nothing

      } else if (IndividualStarTypeIIMassCutoff <= mp &&
                 mp < IndividualStarPSNMassCutoff){

        // mass and energy to distribute over local cells
        distmass = StarMassEjectionFraction * ParticleMass[i] / ((float) StarFeedbackDistTotalCells);
        energy   = ParticleMass[i] * StarEnergyToThermalFeedback * 
                   (speed_of_light*speed_of_light/((*v1)*(*v1))) / ((float) StarFeedbackDistTotalCells);

        // add energy to surroundings
        for(int kc = kp - StarFeedbackDistRadius; kc <= kp + StarFeedbackDistRadius; kc++){
          int stepk = abs(kc - kp);
          for(int jc = jp - StarFeedbackDistRadius; jc <= jp + StarFeedbackDistRadius; jc++){
            int stepj = stepk + abs(jc - jp);
            for(int ic = ip - StarFeedbackDistRadius; ic <= ip + StarFeedbackDistRadius; ic++){
              int cellstep = stepj + abs(ic - ip);
              index = ic + (jc + kc * (*ny)) * (*nx);


              if (cellstep <= StarFeedbackDistCellStep){ // not sure
                dratio    = 1.0 / (BaryonField[DensNum][index] + distmass);
                
// add CR here:
                if( CRModel ){
                  BaryonField[TENum][index] =
                     (BaryonField[TENum][index]*BaryonField[DensNum][index] +
                                                     energy*(1.0-CRFeedback)) * dratio;

                  if(DualEnergyFormalism){
                    BaryonField[GENum][index] =
                      ((BaryonField[GENum][index] * BaryonField[DensNum][index]) +
                                                     energy*(1.0-CRFeedback)) * dratio;
                  }


                  BaryonField[CRNum][index] += energy*CRFeedback; // add CR
                } else{
                  BaryonField[TENum][index] =
                     (BaryonField[TENum][index]*BaryonField[DensNum][index] + energy) * dratio;

                  if(DualEnergyFormalism){
                    BaryonField[GENum][index] =
                        ((BaryonField[GENum][index] * BaryonField[DensNum][index]) +
                                                                               energy) * dratio;
                  }
                } // end CR

                // now do metal feedback and chemical yields
                // here
                // if(TestProblemData.UseMetallicityField){
                //
                //  if(MULTIMETALS_METHOD(MULTIMETALS_ALPHA){
                //
                //    }
                //  BaryonField[CINum][index] = BaryonField[CINum][index] * 1.01;
                //  BaryonField[NINum][index] = BaryonField[NINum][index] * 1.05;
                //  BaryonField[OINum][index] = BaryonField[OINum][index] * 0.9;
                //}

                // do mass and momentum feedback
                BaryonField[Vel1Num][index] = BaryonField[Vel1Num][index] * BaryonField[DensNum][index] +
                           distmass * ParticleVelocity[0][i];
                BaryonField[Vel2Num][index] = BaryonField[Vel2Num][index] * BaryonField[DensNum][index] +
                           distmass * ParticleVelocity[1][i];
                BaryonField[Vel3Num][index] = BaryonField[Vel3Num][index] * BaryonField[DensNum][index] +
                           distmass * ParticleVelocity[2][i];

                // add mass to cells and convert vels to des again
                BaryonField[DensNum][index] += distmass;
                BaryonField[Vel1Num][index] /= BaryonField[DensNum][index];
                BaryonField[Vel2Num][index] /= BaryonField[DensNum][index];
                BaryonField[Vel3Num][index] /= BaryonField[DensNum][index];

                if (HydroMethod != 2 && DualEnergyFormalism){
                  BaryonField[TENum][index] = 0.5 * (BaryonField[Vel1Num][index]*BaryonField[Vel1Num][index] + 
                                                     BaryonField[Vel2Num][index]*BaryonField[Vel2Num][index] +
                                                     BaryonField[Vel3Num][index]*BaryonField[Vel3Num][index])
                                                   + BaryonField[GENum][index];
                }

              } // end if distribute


            } // i dist
          } // j dist
        }// k dist


        // particle is now dead
      } // end if Type II SN

      // kill particle
      ParticleMass[i] = 0.0;
    } // end end of life feedback


  } // loop over particles

  return SUCCESS;
}

//
// Function to compute star lifetime based on age
// and metallicity. Input mass in Msun.
// Age in cgs
float compute_lifetime(float *mp){

  const double myr       = 3.1536e13;
  const double tau_solar = myr * 1.0E4   ; // lifetime of 1 solar mass star, seconds

  float tau;

  // For now, just do an approximate scaling:
  // until have full stellar evolution models
  tau = tau_solar * POW( *mp, -2.5);

  return tau;
}
