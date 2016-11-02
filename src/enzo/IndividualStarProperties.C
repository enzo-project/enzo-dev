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



#include "StarParticleData.h"
#include "IndividualStarProperties.h"

// some global constants
const double SOLAR_LIFETIME    = 10.0E9 * 3.1536E7 ; // s
const double SOLAR_LUMINOSITY  = 3.9E33            ; // cgs
const double SOLAR_TEFF        = 5777.0            ; // K - Cox 2000
//const double SOLAR_MASS        = 1.989E33           ; // g
const double SOLAR_RADIUS      = 69.63E9           ; // cm
const double SOLAR_METALLICITY = 0.02              ;




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


unsigned_long_int mt_random(void);


float SNIaProbability(const float &current_time, const float &formation_time,
                      const float &lifetime, const float &TimeUnits){


  float  dPdt = IndividualStarSNIaFraction;
  const float hubble_time = 4.408110831E17;


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

int SetWDLifetime(float &WD_lifetime,
                    const float &current_time, const float &formation_time,
                    const float &lifetime, const float &TimeUnits){

  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  float x  = (float) (random_int % max_random) / ((float) (max_random));

  const int size = 1000;
  const float hubble_time = 4.408110831E17;

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

    int width = size / 2;
    int bin_number = size / 2;

    while ( width > 1){
      width /= 2;

      if (x > tabulated_probability[bin_number]){
        bin_number += width;
      } else if (x < tabulated_probability[bin_number]){
        bin_number -= width;
      } else{
        break;
      }
    }

    WD_lifetime = time[bin_number];

    return 1;
  }

  return FAIL;
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
  const double G = 6.6743E-8;

  return G*(mp * SOLAR_MASS) / ( R*R );
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
    if (t < 0){
      return FAIL;
    }

    printf("IndividualStarInterpolateLWFlux: Failure in interpolation");

    if( t < 0) {
      printf("Temperature out of bounds ");
      value = Teff; value_min = IndividualStarRadData.T[0]; value_max = IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins-1];
    } else if (u < 0) {
      printf("Surface Gravity out of bounds ");
      value = g; value_min = IndividualStarRadData.g[0]; value_max = IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins-1];
    } else if (v < 0) {
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

    if ( t < 0){
      return FAIL;
    } // no print statements if temperature is out of bounds -
      // that is O.K. and covered with black body integral

    printf("IndividualStarInterpolateFUVFlux: Failure in interpolation "); 


    if( t < 0) {
      printf("Temperature out of bounds ");
      value = Teff; value_min = IndividualStarRadData.T[0]; value_max = IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins-1];
    } else if (u < 0) {
      printf("Surface Gravity out of bounds ");
      value = g; value_min = IndividualStarRadData.g[0]; value_max = IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins-1];
    } else if (v < 0) {
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

int IndividualStarComputeLWLuminosity(float &L_Lw, Star *cstar){

  float Teff, R, g, LW_flux;
  const double pi = 3.141592653689393;

  int * se_table_pos = cstar->ReturnSETablePosition();
  int * rad_table_pos = cstar->ReturnRadTablePosition();

  IndividualStarInterpolateProperties(Teff, R,
                                      se_table_pos[0], se_table_pos[1],
                                      cstar->ReturnBirthMass(), cstar->ReturnMetallicity());

  g = IndividualStarSurfaceGravity(cstar->ReturnBirthMass(), R);

  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateLWFlux(LW_flux,
                                       rad_table_pos[0], rad_table_pos[1], rad_table_pos[2],
                                       Teff, g, cstar->ReturnMetallicity()) == SUCCESS){

      L_Lw = 4.0 * pi * R * R * LW_flux;

      return SUCCESS; // rates computed from table
    }
  }

  /* if we are here it is because we need to do the black body integration */
  const double lw_emin = 11.2 * 1.6021772E-12 ; // eV -> cgs
  const double lw_emax = 13.6 * 1.6021772E-12 ; // eV -> cgs

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

int IndividualStarComputeFUVLuminosity(float &L_fuv, Star *cstar){

  float Teff, R, g, Fuv;
  const double pi = 3.141592653689393;

  int * se_table_pos = cstar->ReturnSETablePosition();
  int * rad_table_pos = cstar->ReturnRadTablePosition();

  IndividualStarInterpolateProperties(Teff, R,
                                      se_table_pos[0], se_table_pos[1],
                                      cstar->ReturnBirthMass(), cstar->ReturnMetallicity());

  g = IndividualStarSurfaceGravity(cstar->ReturnBirthMass(), R);

  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateFUVFlux(Fuv,
                                        rad_table_pos[0], rad_table_pos[1], rad_table_pos[2],
                                        Teff, g, cstar->ReturnMetallicity()) == SUCCESS){

      L_fuv = 4.0 * pi * R * R * Fuv;

      return SUCCESS; // rates computed from table
    }
  }

  /* if we are here it is because we need to do the black body integration */
  const double fuv_emin =  6.0 * 1.6021772E-12 ; // eV -> cgs
  const double fuv_emax = 13.6 * 1.6021772E-12 ; // eV -> cgs

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
  const double pi = 3.141592653589393;

  IndividualStarInterpolateProperties(Teff, R, mp, metallicity);
  g = IndividualStarSurfaceGravity(mp, R);

  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateFUVFlux(Fuv, Teff, g, metallicity) == SUCCESS){

      L_fuv = 4.0 * pi * R * R * Fuv;

      return SUCCESS; // rates computed from table
    }
  }


  /* if we are here it is because we need to do the black body integration */
  const double fuv_emin =  6.0 * 1.6021772E-12 ; // eV -> cgs
  const double fuv_emax = 13.6 * 1.6021772E-12 ; // eV -> cgs

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

  const double E_ion_HI   = 13.6; // eV
  const double E_ion_He   = 24.587; // eV
  const double k_boltz = 1.380658E-16;
  const double c       = 2.99792458E10;
  const double h       = 6.6260755E-27;
  const double eV_erg  = 6.24150934326E11; // eV in one erg


  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateRadData(q0, q1,
                                        i, j, k,
                                        Teff, g, metallicity) == SUCCESS){
    }

  }

  // if we are here, it means star was outside tabulated data. Use a black
  // body to compute radiation:
  float x;

  x = (E_ion_HI / eV_erg) / (k_boltz * (Teff));
  // compute the black body radiance in unitless numbers
  if( (PhotonRadianceBlackBody(q0, x)) == FAIL){


    ENZO_FAIL("IndividualStarComputeIonizingRates: Summation of black body integral failed to converge\n");
  }

  x = (E_ion_He / eV_erg) / (k_boltz * (Teff));
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


  float A = 2.0 * k_boltz * k_boltz * k_boltz * (Teff*Teff*Teff) /
                               (h*h*h*c*c);

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
  const double E_ion_HI   = 13.6; // eV
  const double E_ion_He   = 24.587; // eV
  const double k_boltz = 1.380658E-16;
  const double c       = 2.99792458E10;
  const double h       = 6.6260755E-27;
  const double eV_erg  = 6.24150934326E11; // eV in one erg


  if (IndividualStarBlackBodyOnly == FALSE){
    if(IndividualStarInterpolateRadData(q0, q1, Teff, g, metallicity) == SUCCESS){
      return SUCCESS; // rates computed from table
    }
  }

  // if we are here, it means star was outside tabulated data. Use a black
  // body to compute radiation:
  float x;

  x = (E_ion_HI / eV_erg) / (k_boltz * (Teff));
  // compute the black body radiance in unitless numbers
  if( (PhotonRadianceBlackBody(q0, x)) == FAIL){
    ENZO_FAIL("IndividualStarComputeIonizingRates: Summation of black body integral failed to converge\n");
  }

  x = (E_ion_He / eV_erg) / (k_boltz * (Teff));
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


  float A = 2.0 * k_boltz * k_boltz * k_boltz * (Teff*Teff*Teff) /
                               (h*h*h*c*c);

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

  // WARNING: see statement at top of file
  if (Z < IndividualStarPropertiesData.Z[0]){
    Z = IndividualStarPropertiesData.Z[0];
  }

  if( LinearInterpolationCoefficients(t, u, i, j, M, Z,
                                      IndividualStarPropertiesData.M, IndividualStarPropertiesData.Z,
                                      IndividualStarPropertiesData.NumberOfMassBins, IndividualStarPropertiesData.NumberOfMetallicityBins) == FAIL){
    /* if interpolation fails, fail here */
    float value, value_min, value_max;
    printf("IndividualStarInterpolateProperties: Failure in interpolation ");

    if( t < 0){
      printf("Mass out of bounds ");
      value = M; value_min = IndividualStarPropertiesData.M[0]; value_max = IndividualStarPropertiesData.M[IndividualStarPropertiesData.NumberOfMassBins-1];
    } else if (u < 0){
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

int IndividualStarEvaluateInterpolation(float &y, float *ya[],
                                        const int &i, const int &j,
                                        const float &t, const float &u){

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
  float Z; // Z is in units of solar
  Z = metallicity;

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
    ENZO_FAIL("IndividualStarInterpolateProperties: Failure in lifetime interpolation, mode must be set to either 1 or 2");
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
    if ( t < 0 ){
      i = -9; j = -9; k = -9;
      return SUCCESS;
    }

    printf("IndividualStarInterpolateRadData: Failure in interpolation ");

    if( t < 0) { 
      printf("Temperature out of bounds ");
      value = Teff; value_min = IndividualStarRadData.T[0]; value_max = IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins-1];
    } else if (u < 0) {
      printf("Surface Gravity out of bounds ");
      value = g; value_min = IndividualStarRadData.g[0]; value_max = IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins-1];
    } else if (v < 0) {
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
    printf("IndividualStarInterpolateProperties: Failure in interpolation ");

    if( t < 0){
      printf("Mass out of bounds ");
      value = M;
      value_min = IndividualStarPropertiesData.M[0];
      value_max = IndividualStarPropertiesData.M[IndividualStarPropertiesData.NumberOfMassBins-1];
    } else if (u < 0){
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
    printf("IndividualStarInterpolateProperties: Failure in interpolation ");

    if( t < 0){ 
      printf("Mass out of bounds ");
      value = M; value_min = IndividualStarPropertiesData.M[0]; value_max = IndividualStarPropertiesData.M[IndividualStarPropertiesData.NumberOfMassBins-1];
    } else if (u < 0){
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
    printf("IndividualStarInterpolateProperties: Failure in interpolation ");

    if( t < 0){ 
      printf("Mass out of bounds ");
      value = M; value_min = IndividualStarPropertiesData.M[0]; value_max = IndividualStarPropertiesData.M[IndividualStarPropertiesData.NumberOfMassBins-1];
    } else if (u < 0){
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
    
    if ( t < 0 ){
        return FAIL; // Temperature failure is O.K. --- just means do black body
    }

    printf("IndividualStarInterpolateRadData: Failure in interpolation ");

    if( t < 0) { 
      printf("Temperature out of bounds ");
      value = Teff; value_min = IndividualStarRadData.T[0]; value_max = IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins-1];
    } else if (u < 0) {
      printf("Surface Gravity out of bounds ");
      value = g; value_min = IndividualStarRadData.g[0]; value_max = IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins-1];
    } else if (v < 0) {
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

  const double c       = 2.99792458E10;
  const double k        = 1.380658E-16;
  const double h       = 6.6260755E-27;

  float x1, x2, A;

  x1 = (e_min) / (k * Teff);
  x2 = (e_max) / (k * Teff);

  BlackBodyFlux(flux, x1, x2);

  A = 2.0 * k*k*k*k * Teff*Teff*Teff*Teff / (h*h*h*c*c);

  flux *= A;

  return SUCCESS;
}


int ComputeAverageEnergy(float *energy, float *e_i, float *Teff){
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

  const double sigma   = 5.6074E-5;
  const double pi      = 3.1415621;
  const double c       = 2.99792458E10;
  const double k_boltz = 1.380658E-16;
  const double h       = 6.6260755E-27;

  // compute the number of photons
  float xmax;

  // convert wavelength to unitless energy
  xmax = *e_i /(k_boltz*(*Teff));

//  printf("ISP: Energy %"ESYM" xmax %"ESYM" Teff %"ESYM"\n",*e_i, xmax, *Teff);
  if(AverageEnergyBlackBody(energy, xmax)==FAIL){
    printf("Warning: Non-convergence in black body integral (summation) for IndividualStar Spectrum\n");
  }

//  printf("ISP: Avg energy unitless %"ESYM"\n", *energy);
  // converto from unitless energy to cgs
  *energy = (*energy) * (k_boltz) * (*Teff);

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

  float difference = 1.0E10, sum, old_sum;
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

  difference = 1.0E10;
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

int AverageEnergyBlackBody(float *energy, float x){
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

  difference = 1.0E10;
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
  difference = 1.0E10;
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
  *energy = u_dens_summation / n_dens_summation;

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

  /* now perform a binary search over the array */
  width = x1a_size / 2;
  i     = x1a_size / 2;

  while (width > 1){
    width /= 2;
    if ( x1 > x1a[i] )
      i += width;
    else if (x1 < x1a[i] )
      i -= width;
    else
      break;
  }

  // ensure we are finding the nearest, value that is less than x1
  if ( x1 < x1a[i] ) i--;
  if ( x1 < x1a[i] ) i--;

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





