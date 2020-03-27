/*--------------------------------------------------------------------------------
/
/ INDIVIDUAL STAR PROPERTIES
/
/ Author: Andrew Emerick
/ Date  : March, 2016
/ modified 1:
/
/ PURPOSE:
/   Contains helper routines to compute star properties used to calculate
/   radiation properties. Includes interpolation method for trillinear
/   interpolation of radiation data as well as a gaussian random variable
/   generator used here and in some random draws in individual_star_maker
/
/   WARNING: In every case where interpolation is over metallicity, we institute
/            a 'metallicity floor' whereby the star's metallicity is set to Z_min
/            during interpolation for that given table, whenever Z_star < Z_min
/            for that given table. This *may* lead to potential inconsistencies if
/            Z_min is different in each table, but really that means certain
/            properties will be *better* than others in this regime. This is also
/            applied to the nucleosynethsis yields tables.
/
/
--------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "flowdefs.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "communication.h"
#include "CommunicationUtilities.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"


#include "StarParticleData.h"
#include "IndividualStarProperties.h"
#include "StellarYieldsRoutines.h"



int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
             
extern "C" void FORTRAN_NAME(pop3_properties)(FLOAT *mass, FLOAT* luminosity,
                                              FLOAT *lifetime);

/* internal function prototypes */
int IndividualStarRadDataEvaluateInterpolation(float &y, float **ya[],
                                               const float &t, const float &u, const float &v,
                                               const int &i, const int &j, const int &k);

float LinearInterpolationCoefficient(const int &i, const float &x1, const float *x1a);

int LinearInterpolationCoefficients(float &t, int &i,
                                    const float &x1, const float *x1a, const int &x1a_size);

int LinearInterpolationCoefficients(float &t, float &u, int &i, int &j,
                                    const float &x1, const float &x2,
                                    const float *x1a, const float *x2a,
                                    const int &x1a_size, const int &x2a_size);

int LinearInterpolationCoefficients(float &t, float &u, float &v, int &i, int &j, int &k,
                                    const float &x1, const float &x2, const float &x3,
                                    const float *x1a, const float *x2a, const float *x3a,
                                    const int &x1a_size, const int &x2a_size, const int &x3a_size);

int search_lower_bound(float *arr, float value, int low, int high, int total);

unsigned_long_int mt_random(void);

float SNIaProbability(const float &current_time, const float &formation_time,
                      const float &lifetime, const float &TimeUnits){
/* --------------------------------
 * SNIaProbability
 * --------------------------------
 * Using DTD model, computes probability that WD will explode as SNIa
 * at a specified time after its formation
 *
 * CURRENTLY NOT USED IN SNIa MODEL - OCT 2016 - SLATED FOR REMOVAL
 * ---------------------------------
 */

  float  dPdt = IndividualStarSNIaFraction;
  float hubble_time;
  if (ComovingCoordinates){
    hubble_time = (1.0/HubbleConstantNow)*1.0E4 * Myr_s;
  } else{
    hubble_time = (1.0/0.701)*1.0E4*Myr_s;
  }



  if (IndividualStarDTDSlope == 1.0){
    dPdt /= log( ( (hubble_time / TimeUnits) + formation_time ) / (formation_time + lifetime));
  } else {
    dPdt *= (- IndividualStarDTDSlope + 1.0);
    dPdt /= ( POW( (hubble_time / TimeUnits) + formation_time, -IndividualStarDTDSlope + 1) -
              POW( formation_time + lifetime, -IndividualStarDTDSlope + 1));
  }

  dPdt = dPdt * POW(current_time, -IndividualStarDTDSlope);

  return dPdt;
}

int grid::IndividualStarSetWDLifetime(void){
/* --------------------------------------------------
 * IndividualStarSetWDLifetime
 * --------------------------------------------------
 * Updates WD lifetimes if not yet initialized using
 * DTD SNIa model
 * --------------------------------------------------
 */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
              &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }


  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0)
    return SUCCESS;

  for (int i = 0; i < NumberOfParticles; i++){

    if( ParticleType[i] != -PARTICLE_TYPE_INDIVIDUAL_STAR_WD ){
      continue;
    }

    //
    // lifetime is set relative to WD formation time (now)
    //
    float new_lifetime = -1;

    int result = SetWDLifetime(new_lifetime, this->Time, ParticleAttribute[0][i],
                                            ParticleAttribute[1][i], TimeUnits);

    ParticleAttribute[1][i] = new_lifetime;
    ParticleType[i]         = ABS(ParticleType[i]);

    if (ParticleAttribute[1][i] < 0){
      return FAIL;
    }
    //
    // feedback operates computing death time = lifetime + birth time
    // renormalize so as to keep birth time the original star particle birth time
    //  - original lifetime of progenitor star to WD can be backed out via postprocessing, but not birth time
    //
    if (result > 0){ // negative result means WD never exploding
      // set lifetime correctly as new_lifetime + main_sequence_lifetime
      // fmax forces any explosion NOW to happen next 1-2 timesteps
      // since machinery may not catch appropriately
      ParticleAttribute[1][i] = fmax(new_lifetime,1.5*this->dtFixed) + (this->Time - ParticleAttribute[0][i]);
    }
  }


  return SUCCESS;
}

int SetWDLifetime(float &WD_lifetime,
                    const float &current_time, const float &formation_time,
                    const float &lifetime, const float &TimeUnits){
/* --------------------------------------------------------------------------
 * SetWDLifetime
 * --------------------------------------------------------------------------
 * For a given WD, assigns its lifetime based on the formation of the MS star
 * and the lifetime of that star. Computes this using a DTD model to compute
 * when over a hubble time the star will explode as a SNIa, using a random
 * number draw with total normalized probability set as input (on order of
 * a few percent) [see SNIaProbabiltiy function]. If the star does not explode
 * in a Hubble time, then lifetime is set to some arbitrarily large value.
 * --------------------------------------------------------------------------
 */

  if (IndividualStarWDFixedLifetime > 0){ // testing parameter to force SNIA immediately
    WD_lifetime = IndividualStarWDFixedLifetime * Myr_s / TimeUnits;
    return 1;
  }

  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  float x  = (float) (random_int % max_random) / ((float) (max_random));

  const int size = 1000; // number of table entries
  float hubble_time;

  if (ComovingCoordinates){
    hubble_time = (1.0/HubbleConstantNow)*1.0E4 * Myr_s;
  } else{
    hubble_time = (1.0/0.701)*1.0E4*Myr_s;
  }

  float tabulated_probability [size];
  float time [size];

  /* logspace in time */
  float min_time = log10(lifetime/10.0);
  float max_time = log10(hubble_time / TimeUnits);
  float dt       = (max_time - min_time) / (1.0*(size - 1)) ;

  time[0] = POW(10.0, min_time);
  tabulated_probability[0] = SNIaProbability(time[0] + formation_time + lifetime,
                                             formation_time, lifetime, TimeUnits);

  const float inv_six = 1.0 / 6.0;

  float f_a = tabulated_probability[0];

  /* tabulate dPdt and do Simpson's int to obtain probability distribution */
  for( int i = 1; i < size; i ++){
    time[i] = POW(10.0, min_time + dt*i);

    float f_b  = SNIaProbability(time[i] + formation_time + lifetime,
                                 formation_time, lifetime, TimeUnits);

    float f_ab = SNIaProbability( 0.5*(time[i] + time[i-1]) + formation_time + lifetime,
                                 formation_time, lifetime, TimeUnits);

    tabulated_probability[i]  = inv_six*(time[i] - time[i-1])*(f_a + 4.0*f_ab + f_b);

    f_a = f_b;
  }

  /* this is strictly wrong, but should never really be an issue */
  tabulated_probability[0] = tabulated_probability[1] * 0.5;

  // change to cumulative probability distribution
  // This is normalized to IndividualStarSNIa fraction, or
  // the probability that the star will explode as SNIa in a Hubble time
  // -- value is around a few percent
  for( int i = 1; i < size; i ++){
    tabulated_probability[i] += tabulated_probability[i-1];
  }

  /* now use the random number to set lifetime */
  if (x < tabulated_probability[0] ){

    // explosion is happening immenently - This is unlikely

    WD_lifetime = time[0];

    return 1;

  } else if (x > tabulated_probability[size-1]){

    // explosion never happens
    // most likely outcome assuming prob is around a few percent

    WD_lifetime = 1000.0 * hubble_time/TimeUnits;

    return -1;
  } else {
    // explosion is happening at some time, find that time

    int bin_number = search_lower_bound(tabulated_probability, x, 0,
                                        size, size);
    WD_lifetime = time[bin_number];

    return 1;
  }

  return FAIL;
}

void ComputeStellarWindVelocity(Star *cstar, float *v_wind){
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

  int* se_table_position = cstar->ReturnSETablePosition();

  /* get properties */
  IndividualStarInterpolateLuminosity(L, se_table_position[0], se_table_position[1],
                                          cstar->ReturnBirthMass(), cstar->ReturnMetallicity());
  IndividualStarInterpolateProperties(Teff, R, se_table_position[0], se_table_position[1],
                                          cstar->ReturnBirthMass(), cstar->ReturnMetallicity());

  // wind is in units of km / s
  // L - solar units
  // T - Kelvin
  // M - solar units
  // Z - solar units
  *v_wind =   1.23 - 0.30*log10(L) + 0.55 * log10(cstar->ReturnBirthMass())
            + 0.64 * log10(Teff) + 0.13*log10(cstar->ReturnMetallicity()/solar_z);
  *v_wind = POW(10.0, *v_wind);

  return;
}

void ComputeStellarWindMassLossRate(const float &mproj, const float &metallicity,
                                    float *dMdt){
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

  /* get properties */
  if(IndividualStarInterpolateLuminosity(L, mproj, metallicity) == FAIL){
    ENZO_FAIL("ComputeStellarWindMassLossRate: Failed to interpolate luminosity");
  }

  if(IndividualStarInterpolateProperties(Teff, R, mproj, metallicity) == FAIL){
    ENZO_FAIL("ComputeStellarWindMassLossRate: Failed to interpolate stellar properties");
  }

  /* compute logged mass loss rate */
  *dMdt = -24.06 + 2.45 * log10(L) - 1.10 * log10(mproj) + 1.31 * log10(Teff)
                                   + 0.80 * log10(metallicity / solar_z);

  *dMdt = POW(10.0, *dMdt) / yr_s ; // SolarMass / yr -> SolarMass / s
  return;

}

float IndividualStarSurfaceGravity(const float &mp, const float &R){
  /* ----------------------------------------------------
   * IndividualStarSurfaceGravity
   * ----------------------------------------------------
   * A. Emerick - March 2016
   *
   * Given stellar mass (solar) and radius (cm), compute
   * the star's surface gravity (cgs).
   * ---------------------------------------------------*/

  return GravConst*(mp * SOLAR_MASS) / ( R*R );
}

int IndividualStarInterpolateLWFlux(float &LW_flux,
                                    const int &i, const int &j, const int &k,
                                    const float &Teff, const float &g, const float &metallicity){

  float t, u, v;

  // convert metallicity to solar
  float Z = (metallicity) / IndividualStarRadData.Zsolar;

  // WARNING: see statement of metallicity floor at top of file
  if (Z < IndividualStarRadData.Z[0]){
    Z = IndividualStarRadData.Z[0];
  }

  /* not on grid - use black body instead */
  if( i == -9 || j == -9 || k == -9){
    return FAIL;
  }

  // compute coefficients
  t = LinearInterpolationCoefficient(i, Teff, IndividualStarRadData.T);
  u = LinearInterpolationCoefficient(j, g   , IndividualStarRadData.g);
  v = LinearInterpolationCoefficient(k, Z   , IndividualStarRadData.Z);

  // do the interpolation
  if(IndividualStarRadDataEvaluateInterpolation(LW_flux, IndividualStarRadData.LW_flux,
                                                t, u, v, i, j, k) == FAIL){
    printf("IndividualStarLWHeating: outside sample gird points, using black body instead\n");
    return FAIL;
  }

  return SUCCESS;
}


int IndividualStarInterpolateLWFlux(float & LW_flux, const float &Teff, const float &g,
                                    const float &metallicity){

  int i, j, k;    // i,j,k are Teff, g, and Z respectively
  float t, u ,v;

  // convert metallicity to solar
  float Z = (metallicity) / IndividualStarRadData.Zsolar;

  // WARNING: see statement of metallicity floor at top of file
  if (Z < IndividualStarRadData.Z[0]){
    Z = IndividualStarRadData.Z[0];
  }

  if( LinearInterpolationCoefficients(t, u ,v, i ,j, k, Teff, g, Z,
                                      IndividualStarRadData.T, IndividualStarRadData.g,
                                      IndividualStarRadData.Z,
                                      IndividualStarRadData.NumberOfTemperatureBins,
                                      IndividualStarRadData.NumberOfSGBins,
                                      IndividualStarRadData.NumberOfMetallicityBins) ==FAIL){

    float value, value_min, value_max;

    // if interpolation fails due to temperature, do BB interpolation
    // otherwise FAIL
    if (t < 0 || t > 1){
      return FAIL;
    }

    printf("IndividualStarInterpolateLWFlux: Failure in interpolation");

    if( t < 0 || t > 1) {
      printf("Temperature out of bounds ");
      value = Teff; value_min = IndividualStarRadData.T[0]; value_max = IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins-1];
    } else if (u < 0 || u > 1) {
      printf("Surface Gravity out of bounds ");
      value = g; value_min = IndividualStarRadData.g[0]; value_max = IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins-1];
    } else if (v < 0 || v > 1) {
      printf("Metallicity out of bounds ");
      value = Z; value_min = IndividualStarRadData.Z[0]; value_max = IndividualStarRadData.Z[IndividualStarRadData.NumberOfMetallicityBins-1];
    }

    printf(" with value = %"ESYM" for minimum = %"ESYM" and maximum %"ESYM"\n", value, value_min, value_max);
    return FAIL;
  }

  // otherwise, do the interpolation
  if(IndividualStarRadDataEvaluateInterpolation(LW_flux, IndividualStarRadData.LW_flux,
                                                t, u, v, i, j, k) == FAIL){
    printf("IndividualStarLWHeating: outside sample gird points, using black body instead\n");
    return FAIL;
  }

  return SUCCESS;
}

int IndividualStarInterpolateFUVFlux(float &FUV_flux,
                                    const int &i, const int &j, const int &k,
                                    const float &Teff, const float &g, const float &metallicity){

  float t, u, v;

  // convert metallicity to solar
  float Z = (metallicity) / IndividualStarRadData.Zsolar;

  // WARNING: see statement of metallicity floor at top of file
  if (Z < IndividualStarRadData.Z[0]){
    Z = IndividualStarRadData.Z[0];
  }

 /* not on grid - use black body instead */
  if( i == -9 || j == -9 || k == -9){
    return FAIL;
  }

  // compute coefficients
  t = LinearInterpolationCoefficient(i, Teff, IndividualStarRadData.T);
  u = LinearInterpolationCoefficient(j, g   , IndividualStarRadData.g);
  v = LinearInterpolationCoefficient(k, Z   , IndividualStarRadData.Z);

  // do the interpolation
  if(IndividualStarRadDataEvaluateInterpolation(FUV_flux, IndividualStarRadData.Fuv,
                                                t, u, v, i, j, k) == FAIL){
    printf("IndividualStarLWHeating: outside sample gird points, using black body instead\n");
    return FAIL;
  }

  return SUCCESS;
}

int IndividualStarInterpolateFUVFlux(float & Fuv, const float &Teff, const float &g,
                                     const float &metallicity){

  int i, j, k; // i -> Teff , j -> g, k -> Z
  float t, u, v;

  // convert metallicity to solar
  float Z; // Z is in units of solar
  Z = (metallicity) / IndividualStarRadData.Zsolar;

  // WARNING: See statement at beginning of file
  if (Z < IndividualStarRadData.Z[0]){
    Z = IndividualStarRadData.Z[0];
  }

  if( LinearInterpolationCoefficients(t, u, v, i, j, k, Teff, g, Z,
                                      IndividualStarRadData.T, IndividualStarRadData.g, IndividualStarRadData.Z,
                                      IndividualStarRadData.NumberOfTemperatureBins, IndividualStarRadData.NumberOfSGBins, IndividualStarRadData.NumberOfMetallicityBins) == FAIL){
    /* if interpolation fails, fail here */
    float value, value_min, value_max;

    if ( t < 0 || t > 1){
      return FAIL;
    } // no print statements if temperature is out of bounds -
      // that is O.K. and covered with black body integral

    printf("IndividualStarInterpolateFUVFlux: Failure in interpolation ");


    if( t < 0 || t > 1) {
      printf("Temperature out of bounds ");
      value = Teff; value_min = IndividualStarRadData.T[0]; value_max = IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins-1];
    } else if (u < 0 || u > 1) {
      printf("Surface Gravity out of bounds ");
      value = g; value_min = IndividualStarRadData.g[0]; value_max = IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins-1];
    } else if (v < 0 || v > 1) {
      printf("Metallicity out of bounds ");
      value = Z; value_min = IndividualStarRadData.Z[0]; value_max = IndividualStarRadData.Z[IndividualStarRadData.NumberOfMetallicityBins-1];
    }

    printf(" with value = %"ESYM" for minimum = %"ESYM" and maximum %"ESYM"\n", value, value_min, value_max);
    return FAIL;
  }


  if( IndividualStarRadDataEvaluateInterpolation(Fuv, IndividualStarRadData.Fuv,
                                                 t, u, v, i, j, k) == FAIL){
    printf("IndividualStarFUVHeating: outside sampled grid points, using black body instead\n");
    return FAIL;
  }

  return SUCCESS;
}

int IndividualStarComputeLWLuminosity(float &L_Lw, //Star *cstar){
                                      const int &i, const int &j, const int &k,
                                      const float &Teff, const float &R, const float &g,
                                      const float &Z){


  float LW_flux;
/*
  int * se_table_pos = cstar->ReturnSETablePosition();
  int * rad_table_pos = cstar->ReturnRadTablePosition();

  IndividualStarInterpolateProperties(Teff, R,
                                      se_table_pos[0], se_table_pos[1],
                                      cstar->ReturnBirthMass(), cstar->ReturnMetallicity());

  g = IndividualStarSurfaceGravity(cstar->ReturnBirthMass(), R);
*/
  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateLWFlux(LW_flux,
                                       i, j, k,
                                       Teff, g, Z) == SUCCESS){

      L_Lw = 4.0 * pi * R * R * LW_flux;

      return SUCCESS; // rates computed from table
    }
  }

  /* if we are here it is because we need to do the black body integration */
  const double lw_emin = LW_threshold_energy * erg_eV ; // eV -> cgs
  const double lw_emax = HI_ionizing_energy * erg_eV ; // eV -> cgs

  ComputeBlackBodyFlux(LW_flux, Teff, lw_emin, lw_emax);

  if (IndividualStarBlackBodyOnly == FALSE){ // adjust to match OSTAR
      int index;
      // If we are here because star is below temperature limit, we are a
      // low mass star and black body vastly overestimates
      if ( (Teff) < IndividualStarRadData.T[0] ){
          index = 0;
      } else { // else we are too hot and black body is even worse
          index = 1;
      }

      LW_flux *= IndividualStarBlackBodyLWFactors[index];
  }


  L_Lw = 4.0 * pi * R * R * LW_flux;

  return SUCCESS;

}

int IndividualStarComputeFUVLuminosity(float &L_fuv,
                                       const int &i, const int &j, const int &k,
                                       const float & Teff, const float &R, const float &g,
                                       const float & Z){

  float Fuv;
/*
  int * se_table_pos = cstar->ReturnSETablePosition();
  int * rad_table_pos = cstar->ReturnRadTablePosition();

  IndividualStarInterpolateProperties(Teff, R,
                                      se_table_pos[0], se_table_pos[1],
                                      cstar->ReturnBirthMass(), cstar->ReturnMetallicity());

  g = IndividualStarSurfaceGravity(cstar->ReturnBirthMass(), R);
*/

  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateFUVFlux(Fuv,
                                        i, j, k,
                                        Teff, g, Z) == SUCCESS){

      L_fuv = 4.0 * pi * R * R * Fuv;

      return SUCCESS; // rates computed from table
    }
  }

  /* if we are here it is because we need to do the black body integration */
  const double fuv_emin = FUV_threshold_energy * erg_eV ; // eV -> cgs
  const double fuv_emax = HI_ionizing_energy * erg_eV ; // eV -> cgs

  ComputeBlackBodyFlux(Fuv, Teff, fuv_emin, fuv_emax);

  if (IndividualStarBlackBodyOnly == FALSE){ // adjust to match OSTAR
      int index;
      // If we are here because star is below temperature limit, we are a
      // low mass star and black body vastly overestimates
      if ( (Teff) < IndividualStarRadData.T[0] ){
          index = 0;
      } else { // else we are too hot and black body is even worse
          index = 1;
      }

      Fuv *= IndividualStarBlackBodyFUVFactors[index];
  }


  L_fuv = 4.0 * pi * R * R * Fuv;

  return SUCCESS;


}

int IndividualStarComputeFUVLuminosity(float &L_fuv, const float &mp, const float &metallicity){

  float Teff, R, g, Fuv;

  IndividualStarInterpolateProperties(Teff, R, mp, metallicity);
  g = IndividualStarSurfaceGravity(mp, R);

  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateFUVFlux(Fuv, Teff, g, metallicity) == SUCCESS){

      L_fuv = 4.0 * pi * R * R * Fuv;

      return SUCCESS; // rates computed from table
    }
  }


  /* if we are here it is because we need to do the black body integration */
  const double fuv_emin = FUV_threshold_energy * erg_eV ; // eV -> cgs
  const double fuv_emax = HI_ionizing_energy * erg_eV ; // eV -> cgs

  ComputeBlackBodyFlux(Fuv, Teff, fuv_emin, fuv_emax);

  if (IndividualStarBlackBodyOnly == FALSE){ // adjust to match OSTAR
      int index;
      // If we are here because star is below temperature limit, we are a
      // low mass star and black body vastly overestimates
      if ( (Teff) < IndividualStarRadData.T[0] ){
          index = 0;
      } else { // else we are too hot and black body is even worse
          index = 1;
      }

      Fuv *= IndividualStarBlackBodyFUVFactors[index];
  }


  L_fuv = 4.0 * pi * R * R * Fuv;

  return SUCCESS;

}

int IndividualStarComputeIonizingRates(float &q0, float &q1,
                                       const int &i, const int &j, const int &k,
                                       const float &Teff, const float &g, const float &metallicity){

  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateRadData(q0, q1,
                                        i, j, k,
                                        Teff, g, metallicity) == SUCCESS){
      return SUCCESS;
    }

  }

  // if we are here, it means star was outside tabulated data. Use a black
  // body to compute radiation:
  float x;

  x = (HI_ionizing_energy / eV_erg) / (kboltz * (Teff));
  // compute the black body radiance in unitless numbers
  if( (PhotonRadianceBlackBody(q0, x)) == FAIL){
    ENZO_FAIL("IndividualStarComputeIonizingRates: Summation of black body integral failed to converge\n");
  }

  x = (HeI_ionizing_energy / eV_erg) / (kboltz * (Teff));
  if ( PhotonRadianceBlackBody(q1, x) == FAIL ){
    ENZO_FAIL("IndividualStarComputeIonizingRates: Summation of black body integral failed to converge\n");
  }

  //
  // now adjust the black body curve to be roughly continious
  // with OSTAR, but only do this if OSTAR is actually being used
  //
  if(IndividualStarBlackBodyOnly == FALSE){
    int index;
    // If we are here because star is below temperature limit, we are a
    // low mass star and black body overestimates radiation
    if ( (Teff) < IndividualStarRadData.T[0] ){
        index = 0;
    } else { // else we are too hot
        index = 1;
    }

    q0 = (q0) * IndividualStarBlackBodyq0Factors[index];
    q1 = (q1) * IndividualStarBlackBodyq1Factors[index];
  }


  float A = 2.0 * kboltz * kboltz * kboltz * (Teff*Teff*Teff) /
                               (h_planck*h_planck*h_planck*clight*clight);

  // convert to units of # / s / m-2 / sr-1
  q0 = A * (q0);
  q1 = A * (q1);

  return SUCCESS;



}

int IndividualStarInterpolateRadData(float &q0, float &q1,
                                     const int &i, const int &j, const int &k,
                                     const float &Teff, const float &g, const float &metallicity){

  float t, u, v;

  // convert metallicity to solar
  float Z = (metallicity) / IndividualStarRadData.Zsolar;

  // WARNING: see statement of metallicity floor at top of file
  if (Z < IndividualStarRadData.Z[0]){
    Z = IndividualStarRadData.Z[0];
  }

  /* not on grid - use black body instead */
  if( i == -9 || j == -9 || k == -9){
    return FAIL;
  }

  // compute coefficients
  t = LinearInterpolationCoefficient(i, Teff, IndividualStarRadData.T);
  u = LinearInterpolationCoefficient(j, g   , IndividualStarRadData.g);
  v = LinearInterpolationCoefficient(k, Z   , IndividualStarRadData.Z);

  // do the interpolation
  if(((IndividualStarRadDataEvaluateInterpolation(q0, IndividualStarRadData.q0,
                                                t, u, v, i, j, k) == FAIL)) ||
     ((IndividualStarRadDataEvaluateInterpolation(q1, IndividualStarRadData.q1,
                                                t, u, v, i, j, k) == FAIL))){
    printf("IndividualStarRadData: outside sample gird points, using black body instead\n");
    return FAIL;
  }

  return SUCCESS;

}

int IndividualStarComputeIonizingRates(float &q0, float &q1, const float &Teff,
                                       const float &g, const float &metallicity){
 /*============================================================
  * IndividualStarComputeIonizingRates
  * ===========================================================
  * A. Emerick - March 2016
  *
  * Computes the ionizing photon rates for an individual star.
  * Attempts to use OSTAR2002 gridded ionizing rate data first
  * which is an interpolation over a grid of T, g, and Z. If
  * the desired star properties are outside of grid, then
  * use a black body instead
  * ============================================================
  */

  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateRadData(q0, q1, Teff, g, metallicity) == SUCCESS){
      return SUCCESS; // rates computed from table
    }
  }

  // if we are here, it means star was outside tabulated data. Use a black
  // body to compute radiation:
  float x;

  x = (HI_ionizing_energy / eV_erg) / (kboltz * (Teff));
  // compute the black body radiance in unitless numbers
  if( (PhotonRadianceBlackBody(q0, x)) == FAIL){
    ENZO_FAIL("IndividualStarComputeIonizingRates: Summation of black body integral failed to converge\n");
  }

  x = (HeI_ionizing_energy / eV_erg) / (kboltz * (Teff));
  if ( PhotonRadianceBlackBody(q1, x) == FAIL ){
    ENZO_FAIL("IndividualStarComputeIonizingRates: Summation of black body integral failed to converge\n");
  }

  //
  // now adjust the black body curve to be roughly continious
  // with OSTAR, but only do this if OSTAR is actually being used
  //
  if(IndividualStarBlackBodyOnly == FALSE){
    int index;
    // If we are here because star is below temperature limit, we are a
    // low mass star and black body overestimates radiation
    if ( (Teff) < IndividualStarRadData.T[0] ){
        index = 0;
    } else { // else we are too hot
        index = 1;
    }

    q0 = (q0) * IndividualStarBlackBodyq0Factors[index];
    q1 = (q1) * IndividualStarBlackBodyq1Factors[index];
  }


  float A = 2.0 * kboltz * kboltz * kboltz * (Teff*Teff*Teff) /
                               (h_planck*h_planck*h_planck*clight*clight);

   // convert to units of # / s / m-2 / sr-1
  q0 = A * (q0);
  q1 = A * (q1);

  return SUCCESS;
}

int IndividualStarInterpolateLifetime(float &tau, const float &M,
                                      const float &metallicity, const int &mode){
  /* ================================================================
   * IndividualStarInterpolateLifetime
   * ================================================================
   * A. Emerick - May 2016
   *
   * Performs billinear interpolation over star mass and metallicity
   * to compute stellar lifetime using the PARSEC stellar evolution
   * tracks. Returns either the total stellar lifetime (mode = 1)
   * or the main sequence lifetime (mode = 2) of the star in seconds.
   * ================================================================
   */

  int   i, j; // indexes
  float t, u; // interpolation coefficients

  float Z = metallicity;

  if (mode == 3){ /// PopIII
    float temp_mass=0.0, temp_luminosity = 0.0;
    temp_mass = M;

    FORTRAN_NAME(pop3_properties)(&temp_mass, &temp_luminosity, &tau);

    tau *= yr_s;
    return SUCCESS;
  }

  // WARNING: see statement at top of file
  if (Z < IndividualStarPropertiesData.Z[0]){
    Z = IndividualStarPropertiesData.Z[0];
  }

  if( LinearInterpolationCoefficients(t, u, i, j, M, Z,
                                      IndividualStarPropertiesData.M, IndividualStarPropertiesData.Z,
                                      IndividualStarPropertiesData.NumberOfMassBins, IndividualStarPropertiesData.NumberOfMetallicityBins) == FAIL){
    /* if interpolation fails, fail here */
    float value, value_min, value_max;
    printf("IndividualStarInterpolateProperties1: Failure in interpolation ");

    if( t < 0 || t > 1){
      printf("Mass out of bounds ");
      value = M; value_min = IndividualStarPropertiesData.M[0]; value_max = IndividualStarPropertiesData.M[IndividualStarPropertiesData.NumberOfMassBins-1];
    } else if (u < 0 || u > 1){
      printf("Metallicity out of bounds ");
      value = Z; value_min = IndividualStarPropertiesData.Z[0]; value_max = IndividualStarPropertiesData.Z[IndividualStarPropertiesData.NumberOfMetallicityBins-1];
    }
    printf(" with value = %"ESYM" for minimum = %"ESYM" and maximum %"ESYM"\n", value, value_min, value_max);
    return FAIL;
  }

  if (mode == 1){
    IndividualStarEvaluateInterpolation(tau,
                                         IndividualStarPropertiesData.lifetime, i,j,t,u);
  } else if (mode == 2){
    IndividualStarEvaluateInterpolation(tau,
                                        IndividualStarPropertiesData.agb_start, i,j,t,u);
  }

  return SUCCESS;
}


void IndividualStarSetCoreCollapseSupernovaProperties(Star *cstar,
                                                      float &m_eject, float &E_thermal, float *metal_mass){
/* -------------------------------------------------------
 * IndividualStarCoreCollapseSupernovaProperties
 * -------------------------------------------------------
 * A. Emerick - Sep 2016
 *
 * Set the ejected mass, energy, and metal masses for a
 * core collapse supernova, given star's birth mass and
 * metallicity.
 * -------------------------------------------------------
 */

  int *yield_table_position = cstar->ReturnYieldTablePosition();

  /* compute total ejected yield */
  if ( IndividualStarFollowStellarYields && MultiMetals == 2){
    // 0 in first argument signifies use CC supernova yield table
    m_eject   = StellarYieldsInterpolateYield(0, yield_table_position[0], yield_table_position[1],
                                              cstar->ReturnBirthMass(), cstar->ReturnMetallicity(), -1);
  } else{
    m_eject   = StarMassEjectionFraction * cstar->ReturnMass();
  }

  /* Fail if we are injecting a second time */
  if (cstar->ReturnSNMassEjected() > 0.0){
    ENZO_FAIL("Somehow ejected SN mass twice for this particle\n");
  }

  /* set thermal energy of explosion */
  if( IndividualStarSupernovaEnergy < 0){
    E_thermal = m_eject * StarEnergyToThermalFeedback * (clight * clight);
  } else{
    E_thermal = IndividualStarSupernovaEnergy * 1.0E51;
  }

  /* metal masses for tracer species */
  if(IndividualStarFollowStellarYields && MultiMetals == 2){
    metal_mass[0] = StellarYieldsInterpolateYield(0, yield_table_position[0], yield_table_position[1],
                                                  cstar->ReturnBirthMass(), cstar->ReturnMetallicity(), 0);

    for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
      metal_mass[1+i] = StellarYieldsInterpolateYield(0, yield_table_position[0], yield_table_position[1],
                                                      cstar->ReturnBirthMass(), cstar->ReturnMetallicity(),
                                                      StellarYieldsAtomicNumbers[i]);
    }
  }

  return;
}

void IndividualStarSetStellarWindProperties(Star *cstar, const float &Time,
                                            const float &dtFixed, const float &TimeUnits,
                                            float &m_eject,
                                            float &E_thermal, float *metal_mass){

  float wind_lifetime, agb_start_time, wind_dt;

  /* New variables to make code slightly cleaner + handle units */
  float mproj        = cstar->ReturnBirthMass();
  float lifetime     = cstar->ReturnLifetime() * TimeUnits;
  float metallicity  = cstar->ReturnMetallicity();
  float particle_age = (Time - cstar->ReturnBirthTime())*TimeUnits;
  float dt           = dtFixed * TimeUnits;

  int *yield_table_position = cstar->ReturnYieldTablePosition();
  int *se_table_position    = cstar->ReturnSETablePosition();

  float m_eject_total = 0.0;

  if( IndividualStarFollowStellarYields && MultiMetals == 2){

    // 1 = wind, -1 = return total mass
    m_eject = StellarYieldsInterpolateYield(1, yield_table_position[0], yield_table_position[1],
                                            mproj, metallicity, -1); // total ejecta mass in SolarMass
    m_eject_total = m_eject;

    wind_lifetime = lifetime;   // CGS units

    if( mproj < IndividualStarAGBThreshold) {

      // 2 at end of argument implies compute start time of AGB phase
      IndividualStarInterpolateLifetime(agb_start_time, se_table_position[0], se_table_position[1],
                                           mproj, metallicity, 2);
      /* sanity check */
      float temp_lifetime;
      IndividualStarInterpolateLifetime(temp_lifetime, se_table_position[0], se_table_position[1], mproj, metallicity, 1);

      wind_lifetime = lifetime - agb_start_time; // CGS Units
      if (wind_lifetime < 0.0){
        printf("WARNING LIFETIME ISSUE --- lifetime = %"ESYM" agb_start = %"ESYM" temp_lifetime = %"ESYM"\n", lifetime, agb_start_time, temp_lifetime);
      }

        //
        // To ensure total (integrated) mass ejected is accurate, make sure we don't overinject
        // mass when winds should only be "ON" for part of a timestep, either at beginning or end
        // of AGB phase, or when AGB phase is unresolved (i.e. AGB time < dt)
        //
/*
        if (particle_age > lifetime && particle_age - dt < lifetime){

          // wind_dt = fmin( fmax(0.0, lifetime - (particle_age - dt)) , lifetime - agb_start_time);
          wind_dt = fmax(0.0, lifetime - (particle_age - dt));
          wind_dt = fmin( wind_dt, lifetime - agb_start_time);

          // printf("wind lifetime mode 1\n");
        } else if (particle_age > agb_start_time && particle_age < lifetime ) {
          wind_dt = fmin(particle_age - agb_start_time,dt); // wind only occurs for part of timestep + star dies

          // printf("wind lifetime mode 2\n");
        } else if (particle_age < agb_start_time && particle_age + dt > lifetime) {
          //
          // AGB phase is unresolved. Set wind timestep to lifetime to do all ejecta this timestep
          //
          wind_dt = wind_lifetime;
          // printf("wind lifetime mode 3\n");
        } else if (particle_age < agb_start_time && particle_age + dt > agb_start_time){
          wind_dt = particle_age + dt - agb_start_time; // wind only occurs for part of timestep
          // printf("wind lifeitme mode 4\n");
        } else{
          wind_dt = fmin( lifetime - agb_start_time, dt);
          // printf("PROBLEM IN AGB WIND PHASE\n");
        }
*/

      // printf("Wind lifetime = %"ESYM" - wind_dt = %"ESYM"  %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",wind_lifetime, wind_dt, lifetime, agb_start_time, particle_age, dt);


      // end AGB check
    } else { // massive stars (constant wind)

        //
        // Check timestep to make sure we don't overinject yields at end of life
        //
        wind_dt = fmin(fmax(lifetime - particle_age, 0.0) , dt);
        if (wind_dt < 0.0){
           wind_dt = dt - (particle_age - lifetime);

           if(abs(wind_dt) > dt){
             printf("DEBUG WARNING: Something very wrong is happending at stellar wind end of life\n");
             wind_dt = 0.001*dt;
           }
        }
    }

    /* Gaurd against cases where agb phase is zero */
    wind_lifetime = (wind_lifetime < tiny_number) ? dt : wind_lifetime;
//    wind_dt       = (wind_dt       < tiny_number) ? dt : wind_dt;
    wind_dt = dt;
    //printf("corrected Wind lifetime = %"ESYM" - wind_dt = %"ESYM"  %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",wind_lifetime, wind_dt, lifetime, agb_start_time, particle_age, dt);

    if (dt == 0){
        m_eject = 0.0;
        // printf("WARNING: ZERO TIME STEP SIZE IN WIND LAUNCHING");
    } else{
        m_eject  /= wind_lifetime ; // average mass loss rate over entire wind lifetime
    }

    // end yields methods
  } else {

    // use model to compute mass loss rate instead

    ComputeStellarWindMassLossRate(mproj, metallicity, &m_eject);

    wind_dt = fmin( fmax(lifetime - particle_age, 0.0), dt);
    if (wind_dt < 0.0){
       wind_dt = dt - (particle_age - lifetime);

       if(abs(wind_dt) > dt){
         printf("DEBUG WARNING: Something very wrong is happending at stellar wind end of life\n");
         wind_dt = 0.001*dt;
       }
    }

    wind_dt = dt;
    m_eject_total = m_eject;

  } // end  checking for yields


  float v_wind;

  if(mproj < IndividualStarAGBThreshold){
    /* no good model for AGB wind - use constant user velocity */
    v_wind = IndividualStarAGBWindVelocity;

  } else  if (IndividualStarStellarWindVelocity < 0){
    ComputeStellarWindVelocity(cstar, &v_wind); // v in km /s

  } else{
    v_wind = IndividualStarStellarWindVelocity; // user chosen, in km/s
  }

  if (v_wind > IndividualStarMaximumStellarWindVelocity &&
      IndividualStarMaximumStellarWindVelocity > 0)
            v_wind = IndividualStarMaximumStellarWindVelocity;

  v_wind *= km_cm; // now in cgs

  /* Now that we have wind lifetime and ejected mass, compute properties of wind*/

  float correction_factor = 1.0;
  m_eject = m_eject * wind_dt; // convert Mdot to M_ej

  /* If mass injection exceeds total, correct this fractionally -
     this should rarely happen, but is here to ensure correct chemical evolution.
     Also, if this is likely the last time step the particle is alive, dump
     everything. Again, this is not physically the best solution, but
     ensures that all of the correct yields get deposited, but sacrifices
     correct temporal injection. Given that dt is generally small, this means
     the time at which yields get injected may only be off by 10^3 - 10^4 years...
     this should be irrelevant for galaxy scale (100-1000 Myr) simulations. */
  if( (cstar->ReturnWindMassEjected() + m_eject > m_eject_total) ||
      (particle_age + 2.5*dt > lifetime) ){

    float old_ejection = m_eject;
    m_eject = fmax(m_eject_total - cstar->ReturnWindMassEjected(), 0.0);
    correction_factor = m_eject / old_ejection;
  }

  float Teff, R; // Need Teff for computing thermal energy of ejecta
  IndividualStarInterpolateProperties(Teff, R,
                                      se_table_position[0], se_table_position[1],
                                      cstar->ReturnBirthMass(), cstar->ReturnMetallicity());

  E_thermal = 1.5 * Teff * (m_eject*SolarMass / (mh)) * kboltz; // current T of wind

  if( v_wind > IndividualStarMaximumStellarWindVelocity * km_cm){ // so we don't waste CPU
    v_wind = IndividualStarMaximumStellarWindVelocity * km_cm;
  }

  E_thermal = E_thermal + 0.5 * (m_eject * SolarMass) * v_wind * v_wind; // assume 100% KE thermalization

  /* finally, compute metal masses if needed */
  float wind_scaling = wind_dt / wind_lifetime  * correction_factor;

//  if (wind_lifetime <= tiny_number){
//    wind_scaling = 0.0;
//  }

  if(IndividualStarFollowStellarYields && MultiMetals==2){

    metal_mass[0] = StellarYieldsInterpolateYield(1, yield_table_position[0], yield_table_position[1],
                                                     mproj, metallicity, 0); // total metal in SolarMass

    for (int i = 0; i < StellarYieldsNumberOfSpecies; i++){
      metal_mass[1 + i] = StellarYieldsInterpolateYield(1, yield_table_position[0], yield_table_position[1],
                                                        mproj, metallicity,
                                                        StellarYieldsAtomicNumbers[i]);
    }

    /* scale for wind dt and lifetime */
    for(int i = 0; i < StellarYieldsNumberOfSpecies+1; i ++){
      metal_mass[i] *= wind_scaling;
    }

  }


  for(int i = 0; i < StellarYieldsNumberOfSpecies+1; i++){
      if(metal_mass[i] < 0.0){
        printf("particle age = %"ESYM" lifetim - age = %"ESYM" dt %"ESYM" %"ESYM"\n", particle_age, lifetime-particle_age, dt, wind_scaling);
        printf("metal mass = %"ESYM" wind_dt = %"ESYM" wind_lifetime = %"ESYM" eject = %"ESYM"\n",metal_mass[i], wind_dt, wind_lifetime, m_eject);
        if(i>0){
            printf("i = %"ISYM" anum = %"ISYM"\n", i, StellarYieldsAtomicNumbers[i-1]);
        } else{
            printf("i = %"ISYM"\n",i);
        }
        cstar->PrintInfo();

        ENZO_FAIL("Negative metal mass in wind setup");
      }
  }


  // done computing stellar wind properties
  return;
}


void IndividualStarSetPopIIISupernovaProperties(Star *cstar, float &m_eject, float &E_thermal, float *metal_mass){
/* -------------------------------------------------------
 * IndividualStarSetPopIIISupernovaProperties
 * -------------------------------------------------------
 * A. Emerick - Oct 2018
 *
 * Set the ejected mass, energy, and metal masses for a
 * PopIII supernova explosion.
 * -------------------------------------------------------
 */

  int *yield_table_position = cstar->ReturnYieldTablePosition();
  float birth_mass = cstar->ReturnBirthMass();
  /* compute total ejected yield */

  if ( ((birth_mass >= TypeIILowerMass) && (birth_mass <= TypeIIUpperMass)) ||
       ((birth_mass >= PISNLowerMass)   && (birth_mass <= PISNUpperMass)) ){

    if ( IndividualStarFollowStellarYields && MultiMetals == 2){
      m_eject   = StellarYieldsInterpolatePopIIIYield(yield_table_position[0],
                                                      birth_mass, -1);
    } else{
      m_eject   = StarMassEjectionFraction * cstar->ReturnMass();
    }

    /* metal masses for tracer species */
    if(IndividualStarFollowStellarYields && MultiMetals == 2){
      metal_mass[0] = StellarYieldsInterpolatePopIIIYield(yield_table_position[0],
                                                          cstar->ReturnBirthMass(),
                                                          0);
      for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
        metal_mass[1+i] = StellarYieldsInterpolatePopIIIYield(yield_table_position[0],
                                                              cstar->ReturnBirthMass(),
                                                              StellarYieldsAtomicNumbers[i]);
      }
    }

    /* Set energy for normal SN */
    if (birth_mass <= TypeIIUpperMass){
      E_thermal = IndividualStarSupernovaEnergy * 1.0E51;
    } else {

      if (PopIIIPISNEnergy < 0.0){
        // taken from pop3_maker.F (heger and woosley??)
        float he_core    = (13.0 / 24.0) * (birth_mass - 20.0);
        float sne_factor = 5.0 + 1.304 * (he_core - 64.0);
        E_thermal = sne_factor * 1.0E51;
      } else {
        E_thermal = PopIIIPISNEnergy * 1.0E51;
      }

    }

  } else {
    m_eject = 0.0; // no yields -- but we should NEVER have to specify this if
                   //   feedback routines are operating correctly
    for (int i = 0; i <= StellarYieldsNumberOfSpecies; i++){
      metal_mass[i] = 0.0;
    }
    E_thermal = 0.0;
  }

  return;
}

void IndividualStarSetTypeIaSupernovaProperties(float &m_eject, float &E_thermal, float *metal_mass){
/* -------------------------------------------------------
 * IndividualStarSetTypeIaSupernovaProperties
 * -------------------------------------------------------
 * A. Emerick - Oct 2016
 *
 * Set the ejected mass, energy, and metal masses for a
 * Type Ia supernova explosion.
 *
 * Current model treats all Type Ia uniformly
 * -------------------------------------------------------
 */



  m_eject = StellarYields_SNIaYieldsByNumber(-1); // total ejecta in solar masses

  /* set energy given user input */
  if( IndividualStarSupernovaEnergy < 0){
    E_thermal = m_eject * StarEnergyToThermalFeedback * (clight * clight);
  } else {
    E_thermal = IndividualStarSupernovaEnergy * 1.0E51;
  }

  /* populate metal species array if needed */
  if (IndividualStarFollowStellarYields && MultiMetals == 2){
    metal_mass[0] = StellarYields_SNIaYieldsByNumber(0); // total metal mass

    for( int i = 0; i < StellarYieldsNumberOfSpecies; i++){
      metal_mass[i+1] = StellarYields_SNIaYieldsByNumber(StellarYieldsAtomicNumbers[i]);
    }
  }

  return;
}

int IndividualStarEvaluateInterpolation(float &y, float *ya[],
                                        const int &i, const int &j,
                                        const float &t, const float &u){
/* ------------------------------------------------------------------------
 * IndividualStarEvaluateInterpolation
 * ------------------------------------------------------------------------
 * Simple billinear interpolation evaluation when all factors are known
 * ------------------------------------------------------------------------
 */

  y = (1.0 - t)*(1.0 - u) * ya[i  ][j  ] +
      (1.0 - t)*(      u) * ya[i  ][j+1] +
      (      t)*(      u) * ya[i+1][j+1] +
      (      t)*(1.0 - u) * ya[i+1][j  ];

  return SUCCESS;
}

int IndividualStarInterpolateLifetime(float   &tau,
                                      const int &i, const int &j,
                                      const float &M, const float &metallicity, const int &mode){
  /* new function - oct 2016 - for new interplation methods */
  // convert metallicity to solar
  float Z = metallicity;

  // WARNING: See warning at beginning of file
  if( Z < IndividualStarPropertiesData.Z[0]){
    Z = IndividualStarPropertiesData.Z[0];
  }

  float t,u;

  t = LinearInterpolationCoefficient(i, M, IndividualStarPropertiesData.M);
  u = LinearInterpolationCoefficient(j, Z, IndividualStarPropertiesData.Z);

  if (mode == 1){ // total star lifetime
    IndividualStarEvaluateInterpolation(tau,
                                        IndividualStarPropertiesData.lifetime, i,j,t,u);

  } else if (mode == 2){ // main sequence lifetime
    IndividualStarEvaluateInterpolation(tau,
                                        IndividualStarPropertiesData.agb_start, i,j,t,u);

  } else{
    ENZO_FAIL("IndividualStarInterpolateProperties2: Failure in lifetime interpolation, mode must be set to either 1 or 2");
    return FAIL;
  }

  return SUCCESS;
}

int IndividualStarGetRadTablePosition(int &i, int &j, int &k,
                                      const float &Teff, const float &g, const float &metallicity){
  // convert metallicity to solar
  float Z; // Z is in units of solar
  Z = (metallicity) / IndividualStarRadData.Zsolar;

  // WARNING: See warning at beginning of file
  if( Z < IndividualStarRadData.Z[0]){
    Z = IndividualStarRadData.Z[0];
  }

  float t, u, v;
  if( LinearInterpolationCoefficients(t, u, v, i, j, k, Teff, g, Z,
                                      IndividualStarRadData.T, IndividualStarRadData.g, IndividualStarRadData.Z,
                                      IndividualStarRadData.NumberOfTemperatureBins, IndividualStarRadData.NumberOfSGBins, IndividualStarRadData.NumberOfMetallicityBins) == FAIL){
    /* if interpolation fails, fail here */
    float value, value_min, value_max;
    if ( t < 0  || t > 1){
      i = -9; j = -9; k = -9;
      return SUCCESS;
    }

    printf("IndividualStarInterpolateRadData: Failure in interpolation ");

    if( t < 0 || t > 1) {
      printf("Temperature out of bounds ");
      value = Teff; value_min = IndividualStarRadData.T[0]; value_max = IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins-1];
    } else if (u < 0 || u > 1) {
      printf("Surface Gravity out of bounds ");
      value = g; value_min = IndividualStarRadData.g[0]; value_max = IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins-1];
    } else if (v < 0 || v > 1) {
      printf("Metallicity out of bounds ");
      value = Z; value_min = IndividualStarRadData.Z[0]; value_max = IndividualStarRadData.Z[IndividualStarRadData.NumberOfMetallicityBins-1];
    }

    printf(" with value = %"ESYM" for minimum = %"ESYM" and maximum %"ESYM"\n", value, value_min, value_max);
    return FAIL;
  }

  return SUCCESS;
}


int IndividualStarGetSETablePosition(int &i, int &j, const float &M, const float &metallicity){

  float t, u; // interpolation coefficients

  float Z = metallicity;

  // WARNING: see statement at beginning of file
  if (Z < IndividualStarPropertiesData.Z[0]){
    Z = IndividualStarPropertiesData.Z[0];
  }

  if( LinearInterpolationCoefficients(t, u, i, j, M, Z,
                                      IndividualStarPropertiesData.M,
                                      IndividualStarPropertiesData.Z,
                                      IndividualStarPropertiesData.NumberOfMassBins,
                                      IndividualStarPropertiesData.NumberOfMetallicityBins) == FAIL){
    /* if interpolation fails, fail here */
    float value, value_min, value_max;
    printf("IndividualStarInterpolateProperties3: Failure in interpolation ");

    if( t < 0 || t > 1){
      printf("Mass out of bounds ");
      value = M;
      value_min = IndividualStarPropertiesData.M[0];
      value_max = IndividualStarPropertiesData.M[IndividualStarPropertiesData.NumberOfMassBins-1];
    } else if (u < 0 || u > 1){
      printf("Metallicity out of bounds ");
      value = Z;
      value_min = IndividualStarPropertiesData.Z[0];
      value_max = IndividualStarPropertiesData.Z[IndividualStarPropertiesData.NumberOfMetallicityBins-1];
    }
    printf(" with value = %"ESYM" for minimum = %"ESYM" and maximum %"ESYM"\n", value, value_min, value_max);
    return FAIL;
  }

  return SUCCESS;
}

void IndividualStarInterpolateLuminosity(float &L, const int &i, const int &j,
                                        const float &M, const float &metallicity){
 /* -----------------------------------------------------
  * IndividualStarInterpolateLuminosity
  * -----------------------------------------------------
  * A. Emerick - Oct 2016
  * -----------------------------------------------------
  */

  float t, u; // interpolation coefficients

  float Z = metallicity;

  // WARNING: see statement at beginning of file
  if (Z < IndividualStarPropertiesData.Z[0]){
    Z = IndividualStarPropertiesData.Z[0];
  }

  t = LinearInterpolationCoefficient(i, M, IndividualStarPropertiesData.M);
  u = LinearInterpolationCoefficient(j, Z, IndividualStarPropertiesData.Z);

  IndividualStarEvaluateInterpolation(L, IndividualStarPropertiesData.L,
                                      i, j, t, u);

  return;
}


int IndividualStarInterpolateLuminosity(float &L, const float &M, const float &metallicity){
  /* =================================================================
   * IndividualStarInterpolateLuminosity
   * =================================================================
   * A. Emerick - March 2016
   *
   * Performs billinear interpolation over star mass and metallicity
   * to compute star luminosity using the PARSEC stellar evolution
   * tracks. Luminosity is used to set the star particle lifetime.
   * Luminosity is returned in SOLAR UNITS
   * =================================================================
   */

  int   i, j; // indexes
  float t, u; // interpolation coefficients

  float Z = metallicity;

  // WARNING: see statement at beginning of file
  if (Z < IndividualStarPropertiesData.Z[0]){
    Z = IndividualStarPropertiesData.Z[0];
  }

  if( LinearInterpolationCoefficients(t, u, i, j, M, Z,
                                      IndividualStarPropertiesData.M, IndividualStarPropertiesData.Z,
                                      IndividualStarPropertiesData.NumberOfMassBins, IndividualStarPropertiesData.NumberOfMetallicityBins) == FAIL){
    /* if interpolation fails, fail here */
    float value, value_min, value_max;
    printf("IndividualStarInterpolateProperties4: Failure in interpolation ");

    if( t < 0 || t > 1){
      printf("Mass out of bounds ");
      value = M; value_min = IndividualStarPropertiesData.M[0]; value_max = IndividualStarPropertiesData.M[IndividualStarPropertiesData.NumberOfMassBins-1];
    } else if (u < 0 || u > 1){
      printf("Metallicity out of bounds ");
      value = Z; value_min = IndividualStarPropertiesData.Z[0]; value_max = IndividualStarPropertiesData.Z[IndividualStarPropertiesData.NumberOfMetallicityBins-1];
    }
    printf(" with value = %"ESYM" for minimum = %"ESYM" and maximum %"ESYM"\n", value, value_min, value_max);
    return FAIL;
  }
  /* Now apply the coefficients and compute the effective temperature */
  L = (1.0 - t)*(1.0 - u) * IndividualStarPropertiesData.L[i  ][j  ] +
      (1.0 - t)*(      u) * IndividualStarPropertiesData.L[i  ][j+1] +
      (      t)*(      u) * IndividualStarPropertiesData.L[i+1][j+1] +
      (      t)*(1.0 - u) * IndividualStarPropertiesData.L[i+1][j  ] ;

  return SUCCESS;
}

void IndividualStarInterpolateProperties(float &Teff, float &R,
                                         const int &i, const int &j,
                                         const float &M, const float &metallicity){
  float Z = metallicity;

  // WARNING: see statement at beginning of file
  if (Z < IndividualStarPropertiesData.Z[0]){
    Z = IndividualStarPropertiesData.Z[0];
  }

  float t,u;

  t = LinearInterpolationCoefficient(i, M, IndividualStarPropertiesData.M);
  u = LinearInterpolationCoefficient(j, Z, IndividualStarPropertiesData.Z);

  IndividualStarEvaluateInterpolation(Teff, IndividualStarPropertiesData.Teff,
                                      i, j, t, u);
  IndividualStarEvaluateInterpolation(R   , IndividualStarPropertiesData.R   ,
                                      i, j, t, u);

  return;
}

int IndividualStarInterpolateProperties(float &Teff, float &R,
                                        const float &M, const float &metallicity){

  /* ==================================================================
   * IndividualStarInterpolateProperties
   * ==================================================================
   * A. Emerick - March 2016
   *
   * Performs billinear interpolation over star mass and metallicity
   * to compute star effective temperature and radius (for surface
   * gravity) using the PARSEC stellar evolution tracks
   *
   * ==================================================================
   */

  int   i, j; // indexes
  float t, u; // interpolation coefficients

  float Z = metallicity;

  // WARNING: see statement at beginning of file
  if (Z < IndividualStarPropertiesData.Z[0]){
    Z = IndividualStarPropertiesData.Z[0];
  }

  if( LinearInterpolationCoefficients(t, u, i, j, M, Z,
                                      IndividualStarPropertiesData.M, IndividualStarPropertiesData.Z,
                                      IndividualStarPropertiesData.NumberOfMassBins, IndividualStarPropertiesData.NumberOfMetallicityBins) == FAIL){
    /* if interpolation fails, fail here */
    float value, value_min, value_max;
    printf("IndividualStarInterpolateProperties5: Failure in interpolation ");

    if( t < 0 || t > 1){
      printf("Mass out of bounds ");
      value = M; value_min = IndividualStarPropertiesData.M[0]; value_max = IndividualStarPropertiesData.M[IndividualStarPropertiesData.NumberOfMassBins-1];
    } else if (u < 0 || u > 1){
      printf("Metallicity out of bounds ");
      value = M; value_min = IndividualStarPropertiesData.Z[0]; value_max = IndividualStarPropertiesData.Z[IndividualStarPropertiesData.NumberOfMetallicityBins-1];
    }
    printf(" with value = %"ESYM" for minimum = %"ESYM" and maximum %"ESYM"\n", value, value_min, value_max);
    return FAIL;
  }


  /* Now apply the coefficients and compute the effective temperature */
  Teff = (1.0 - t)*(1.0 - u) * IndividualStarPropertiesData.Teff[i  ][j  ] +
         (1.0 - t)*(      u) * IndividualStarPropertiesData.Teff[i  ][j+1] +
         (      t)*(      u) * IndividualStarPropertiesData.Teff[i+1][j+1] +
         (      t)*(1.0 - u) * IndividualStarPropertiesData.Teff[i+1][j  ] ;

  R    = (1.0 - t)*(1.0 - u) * IndividualStarPropertiesData.R[i  ][j  ] +
         (1.0 - t)*(      u) * IndividualStarPropertiesData.R[i  ][j+1] +
         (      t)*(      u) * IndividualStarPropertiesData.R[i+1][j+1] +
         (      t)*(1.0 - u) * IndividualStarPropertiesData.R[i+1][j  ] ;



  return SUCCESS;
}


int IndividualStarInterpolateRadData(float &q0, float &q1,
                                     const float &Teff, const float &g, const float &metallicity){
  /* ===================================================================
   * IndividualStarInterpolateRadData
   * ===================================================================
   * A. Emerick - March 2016
   *
   * Performs trillinear interpolation over effective temperature,
   * surface gravity, and metallicity to compute stellar ionizing fluxes
   * q0 and q1. Method follows billinear interpolation outlined in
   * Numerical Recipes, effectively three first order linear interpolations
   *
   *
   * ===================================================================
   */


  int i, j, k; // i -> Teff , j -> g, k -> Z
  float t, u, v;

  // convert metallicity to solar
  float Z; // Z is in units of solar
  Z = (metallicity) / IndividualStarRadData.Zsolar;

  // WARNING: See warning at beginning of file
  if( Z < IndividualStarRadData.Z[0]){
    Z = IndividualStarRadData.Z[0];
  }

  if( LinearInterpolationCoefficients(t, u, v, i, j, k, Teff, g, Z,
                                      IndividualStarRadData.T, IndividualStarRadData.g, IndividualStarRadData.Z,
                                      IndividualStarRadData.NumberOfTemperatureBins, IndividualStarRadData.NumberOfSGBins, IndividualStarRadData.NumberOfMetallicityBins) == FAIL){
    /* if interpolation fails, fail here */
    float value, value_min, value_max;

    if ( t < 0  || t > 1){
        return FAIL; // Temperature failure is O.K. --- just means do black body
    }

    printf("IndividualStarInterpolateRadData: Failure in interpolation ");

    if( t < 0 || t > 1) {
      printf("Temperature out of bounds ");
      value = Teff; value_min = IndividualStarRadData.T[0]; value_max = IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins-1];
    } else if (u < 0 || u > 1) {
      printf("Surface Gravity out of bounds ");
      value = g; value_min = IndividualStarRadData.g[0]; value_max = IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins-1];
    } else if (v < 0 || v > 1) {
      printf("Metallicity out of bounds ");
      value = Z; value_min = IndividualStarRadData.Z[0]; value_max = IndividualStarRadData.Z[IndividualStarRadData.NumberOfMetallicityBins-1];
    }

    printf(" with value = %"ESYM" for minimum = %"ESYM" and maximum %"ESYM"\n", value, value_min, value_max);
    return FAIL;
  }


  if( (IndividualStarRadDataEvaluateInterpolation(q0, IndividualStarRadData.q0,
                                                  t, u, v, i, j, k) == FAIL) ||
      (IndividualStarRadDataEvaluateInterpolation(q1, IndividualStarRadData.q1,
                                                  t, u, v, i, j, k) == FAIL)){

    printf("IndividualStarRadData: outside sampled grid points, using black body instead\n");
    return FAIL;
  }

  return SUCCESS;
}

int IndividualStarRadDataEvaluateInterpolation(float &y, float **ya[],
                                               const float &t, const float &u, const float &v,
                                               const int &i, const int &j, const int &k){

  /* Since the tabulated data is not fully sampled in surface gravity at every
     temperature, need to make sure that interpolation is not occuring outside
     of the tabulated range. Throw a very loud warning if this is happening    */
  if( ya[i  ][j  ][k  ] == 1.0 ||
      ya[i  ][j+1][k  ] == 1.0 ||
      ya[i+1][j+1][k  ] == 1.0 ||
      ya[i+1][j  ][k  ] == 1.0 ||
      ya[i  ][j  ][k+1] == 1.0 ||
      ya[i  ][j+1][k+1] == 1.0 ||
      ya[i+1][j+1][k+1] == 1.0 ||
      ya[i+1][j  ][k+1] == 1.0   ){
    return FAIL;
  }
  y = (1.0 - t)*(1.0 - u)*(1.0 - v) * ya[i  ][j  ][k  ] +
      (1.0 - t)*(      u)*(1.0 - v) * ya[i  ][j+1][k  ] +
      (      t)*(      u)*(1.0 - v) * ya[i+1][j+1][k  ] +
      (      t)*(1.0 - u)*(1.0 - v) * ya[i+1][j  ][k  ] +
      (1.0 - t)*(1.0 - u)*(      v) * ya[i  ][j  ][k+1] +
      (1.0 - t)*(      u)*(      v) * ya[i  ][j+1][k+1] +
      (      t)*(      u)*(      v) * ya[i+1][j+1][k+1] +
      (      t)*(1.0 - u)*(      v) * ya[i+1][j  ][k+1];


  return SUCCESS;
}

int ComputeBlackBodyFlux(float &flux, const float &Teff, const float &e_min, const float &e_max){
  /* ==============================================================================
   * ComputeBlackBodyFlux
   * ==============================================================================
   * Integral over black body curve to get flux between specified photon energy
   * range.
   *
   * Computation is done by integrating black body curve using series approx
   * ==============================================================================*/

  float x1, x2, A;

  x1 = (e_min) / (kboltz * Teff);
  x2 = (e_max) / (kboltz * Teff);

  BlackBodyFlux(flux, x1, x2);

  A = 2.0 * kboltz*kboltz*kboltz*kboltz * Teff*Teff*Teff*Teff /
           (h_planck*h_planck*h_planck*clight*clight);

  flux *= A;

  return SUCCESS;
}


int ComputeAverageEnergy(float &energy, const float &e_i, const float &Teff){
  /* ==========================================================================
   * ComputeAverageEnergy
   * ==========================================================================
   * Compute the average energy of ionizing photons for a given black body of
   * temperature Teff and for a given ionizing energy of e_i. Energy returned
   * by parameter.
   *
   * Computation is done by integrating black body curve using series approx.
   * ==========================================================================
   */


  // compute the number of photons
  float xmax;

  // convert wavelength to unitless energy
  xmax = e_i /(kboltz*(Teff));

//  printf("ISP: Energy %"ESYM" xmax %"ESYM" Teff %"ESYM"\n",*e_i, xmax, *Teff);
  if(AverageEnergyBlackBody(energy, xmax)==FAIL){
    printf("Warning: Non-convergence in black body integral (summation) for IndividualStar Spectrum\n");
  }

//  printf("ISP: Avg energy unitless %"ESYM"\n", *energy);
  // converto from unitless energy to cgs
  energy = (energy) * (kboltz) * (Teff);

//  printf("ISP: Energy cgs %"ESYM"\n", *energy);

  return SUCCESS;
}

int BlackBodyFlux(float &F, const float &x1, const float &x2){
 /* -----------------------------------------------------------
  * BlackBodyFlux
  * -----------------------------------------------------------
  * Computes the in band flux from the black body curve given
  * unitless wavelength x1 and x2 to integrate over.
  * require that x_2 > x_1
  *
  * Integral is taken as two one sided integrals
  * ------------------------------------------------------------
  */

  float f_1, f_2;
  int error1, error2;

  error1 = BlackBodyFlux(f_1, x1);
  error2 = BlackBodyFlux(f_2, x2);

  F = f_1 - f_2;

 return error1 * error2;
}

int BlackBodyFlux(float &F, const float &x){
 /* -----------------------------------------------------------
  * BlackBodyFlux
  * -----------------------------------------------------------
  * Compute the one sided integral over black body curve to get
  * flux contribution from light with wave number x to infinity
  *
  * ------------------------------------------------------------
  */

  const int max_iterations = 513;
  const int min_iterations = 4;
  const float tolerance = 1.0E-10;

  float difference = huge_number, sum, old_sum;
  sum = old_sum = 0.0;
  int i = 1;

  while((difference > tolerance && i < max_iterations) || (i < min_iterations)){
    old_sum = sum;
    sum += (x*x*x / ((float) i) + 3.0*x*x/((float) i*i) +
           6.0 *x / ((float) i*i*i) + 6.0/((float) i*i*i*i)) * exp(-i*x);

    difference = sum - old_sum;
    i++;
  }

  if ( i >= max_iterations){
    return FAIL;
  }

  F = sum;

  return SUCCESS;
}


int PhotonRadianceBlackBody(float &q, const float &x){
  /* ==========================================================================
   * PhotonRadianceBlackBody
   * ==========================================================================
   * Compute the photon radiance for a black body using unitless wavelength "x",
   * where x*lambda = hc/kT
   *
   * Integral is taken in its series form, which is iterated over until
   * convergence (should only take a few iterations)
   * ==========================================================================
   */

  const int max_iterations = 513;
  const int min_iterations = 4;
  const float tolerance = 1.0E-10;

  float difference, old_sum, sum;
  int i;

  difference = huge_number;
  sum = old_sum = 0.0;
  i = 1;
  // do until tolerance reached, force minimum number of iterations
  while((difference > tolerance && i < max_iterations) || i <min_iterations){
    old_sum = sum;
    sum += ( x*x/((float) i) + 2.0 *x / ((float) i*i) + 2.0/((float) i*i*i) )*exp(-i*x);
    difference = sum - old_sum;
    i++;
  }

  if( i >= max_iterations){
    return FAIL;
  }

  q = sum;

  return SUCCESS;
}

int AverageEnergyBlackBody(float &energy, const float & x){
 //=========================================================================
 // computes the one-sided integral over the black body spectrum
 // to compute the average photon energy for wavelengths lambda and greater,
 // taken as the ratio of the total energy density and num. density of photons
 // using the series approximation of the black body integrals
 // returned unitless energy (E/kT)
 //============================================================================

  const int max_iterations = 513;
  const int min_iterations = 4;
  const float tolerance = 1.0E-10;

  float difference, old_sum, sum, u_dens_summation, n_dens_summation;
  int i;

  difference = huge_number;
  sum = old_sum = 0.0;
  i = 1;
  // compute series approximation of the energy density integral
  // force a minimum number of iterations
  while((difference > tolerance && i < max_iterations) || i < min_iterations){
    old_sum = sum;
    sum += ( x*x*x/((float) i) + 3.0*x*x/((float) i*i)
                     + 6.0*x/((float) i*i*i) + 6.0 / ((float) i*i*i*i)) * exp(-i*x);
    difference = sum - old_sum;
    i++;
  }

  if( i >= max_iterations){
    return FAIL;
  }

  // save value
  u_dens_summation = sum;

  // now compute the series approximation of photon number density integral
  old_sum    = 0.0;
  difference = huge_number;
  sum        = 0.0;
  i = 1;
  // force minimum number of iterations
  while((difference > tolerance && i < max_iterations) || i < min_iterations){
    old_sum = sum;
    sum += (x*x/((float) i) + 2.0*x/((float) i*i) + 2.0/((float) i*i*i)) * exp(-i*x);

    difference = old_sum - sum;
    i++;
  }

  // save value
  n_dens_summation = sum;

  // assign value for unitless energy
  energy = u_dens_summation / n_dens_summation;

//  printf("ISP: nsum %"ESYM" dsum %"ESYM" energy %"ESYM"\n",n_dens_summation, u_dens_summation, *energy);

  if( i >= max_iterations){
    return FAIL;
  }

  return SUCCESS;
}

float LinearInterpolationCoefficient(const int &i,
                                     const float &x1, const float *x1a){
 /* -----------------------------------------------------------
  * LinearInterpolationCoefficient
  * -----------------------------------------------------------
  * If table position is known, compute fractional position
  * between bins.
  * -----------------------------------------------------------
  */

  return (x1 - x1a[i]) / (x1a[i+1] - x1a[i]);
}

int LinearInterpolationCoefficients(float &t, int &i,
                                    const float &x1, const float *x1a, const int &x1a_size){
  /* ------------------------------------------------------------
   * LinearInterpolationCoefficients
   * ------------------------------------------------------------
   * Performs binary search of x1 in assumed sorted array x1a
   * to compute linear interpolation coefficients for that value.
   * ------------------------------------------------------------
   */

  int width, bin_number;

  /* make sure value is within the bounds - FAIL if not */
  if( ( x1 < x1a[0] ) || ( x1 > x1a[x1a_size -1] )){
      t = -1;
      return FAIL;
  }

  i = search_lower_bound((float*)x1a, x1, 0, x1a_size, x1a_size);

  if (i < 0 || i > x1a_size){ ENZO_FAIL("FAILURE IN INTERPOLATION COEFFICIENTS"); }

  t = LinearInterpolationCoefficient(i, x1, x1a);

  return SUCCESS;
}


int LinearInterpolationCoefficients(float &t, float &u, int &i, int &j,
                                    const float &x1, const float &x2,
                                    const float *x1a, const float *x2a,
                                    const int &x1a_size, const int &x2a_size){
  /* ------------------------------------------------------------
   * LinearInterpolationCoefficients
   * ------------------------------------------------------------
   * Overloaded function to do billinear interpolation coefficients
   * ------------------------------------------------------------
   */
  int error1, error2;

  error1 = LinearInterpolationCoefficients(t, i, x1, x1a, x1a_size);
  error2 = LinearInterpolationCoefficients(u, j, x2, x2a, x2a_size);

  return error1 * error2; // FAIL is 0, SUCCESS 1
}

int LinearInterpolationCoefficients(float &t, float &u, float &v, int &i, int &j, int &k,
                                    const float &x1, const float &x2, const float &x3,
                                    const float *x1a, const float *x2a, const float *x3a,
                                    const int &x1a_size, const int &x2a_size, const int &x3a_size){
  /* ------------------------------------------------------------
   * LinearInterpolationCoefficients
   * ------------------------------------------------------------
   * Overloaded function to do trillinear interpolation coefficients
   * ------------------------------------------------------------
   */

  int error1, error2, error3;

  error1 = LinearInterpolationCoefficients(t, i, x1, x1a, x1a_size);
  error2 = LinearInterpolationCoefficients(u, j, x2, x2a, x2a_size);
  error3 = LinearInterpolationCoefficients(v, k, x3, x3a, x3a_size);

  return error1 * error2 * error3; // FAIL is 0, SUCCESS 1
}



float GaussianRandomVariable(void){
/*-----------------------------------------------------------------------------
 GaussianRandomVariable

 Returns a random variable selected over a Gaussian distribution with mean
 of zero and varince of unity using the Box-Muller transform
-----------------------------------------------------------------------------*/
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
