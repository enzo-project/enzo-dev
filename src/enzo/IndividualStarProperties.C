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
--------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"

#include "IndividualStarProperties.h"

// some global constants
const double SOLAR_LIFETIME    = 10.0E9 * 3.1536E7 ; // s
const double SOLAR_LUMINOSITY  = 3.9E33            ; // cgs
const double SOLAR_TEFF        = 5777.0            ; // K - Cox 2000
const double SOLAR_MASS        = 1.99E33           ; // g
const double SOLAR_RADIUS      = 69.63E9           ; // cm
const double SOLAR_METALLICITY = 0.02              ;



/* internal function prototypes */
int IndividualStarRadDataEvaluateInterpolation(float &y, float **ya[],
                                               const float &t, const float &u, const float &v,
                                               const int &i, const int &j, const int &k);

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


float IndividualStarLifetime(const float &mp, const float &metallicity){
  /* ========================================================
   * IndividualStarLifetime
   * --------------------------------------------------------
   * A. Emerick - March 2016
   *
   * Compute star lifetime by interpolating tables to obtain L
   * and scale based upon solar lifetime
   *
   * mp is assumed to be in solar masses
   * =======================================================*/

  float L, lifetime;

  if(IndividualStarInterpolateLuminosity(L, mp, metallicity) == FAIL){
    ENZO_FAIL("IndividualStarLifetime: Failed to interpolate luminosity \n");
  }

  // L from above is in solar units already
  lifetime = SOLAR_LIFETIME * (mp) / (L);

  return lifetime;
}

float IndividualStarLuminosity(const float &mp, const float &lifetime){
  /* -----------------------------------------------------
   * IndividualStarLuminosity
   * -----------------------------------------------------
   * A. Emerick - March 2016
   *
   * Since lifetime is computed using the luminosity, return
   * luminosity given mass and lifetime to avoid having to store
   * the value for luminosity in the simulation
   *
   * M in solar, lifetime in cgs
   * ------------------------------------------------------ */

  return SOLAR_LUMINOSITY * (mp)/ ( lifetime / SOLAR_LIFETIME);
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

int IndividualStarInterpolateFUVFlux(float & Fuv, const float &Teff, const float &g, const float &metallicity){



  return SUCCESS;
}

int IndividualStarComputeFUVLuminosity(float &L_fuv, const float &mp, const float &metallicity){

  float Teff, R, g, Fuv;
  const double pi = 3.1415621;

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


  L_fuv = 4.0 * pi * R * R * Fuv;

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

  if( LinearInterpolationCoefficients(t, u, v, i, j, k, Teff, g, Z, 
                                      IndividualStarRadData.T, IndividualStarRadData.g, IndividualStarRadData.Z,
                                      IndividualStarRadData.NumberOfTemperatureBins, IndividualStarRadData.NumberOfSGBins, IndividualStarRadData.NumberOfMetallicityBins) == FAIL){
    /* if interpolation fails, fail here */
    float value, value_min, value_max;
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

  printf("ISP: Energy %"ESYM" xmax %"ESYM" Teff %"ESYM"\n",*e_i, xmax, *Teff);
  if(AverageEnergyBlackBody(energy, xmax)==FAIL){
    printf("Warning: Non-convergence in black body integral (summation) for IndividualStar Spectrum\n");
  }

  printf("ISP: Avg energy unitless %"ESYM"\n", *energy);
  // converto from unitless energy to cgs
  *energy = (*energy) * (k_boltz) * (*Teff);

  printf("ISP: Energy cgs %"ESYM"\n", *energy);

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

  printf("ISP: nsum %"ESYM" dsum %"ESYM" energy %"ESYM"\n",n_dens_summation, u_dens_summation, *energy);

  if( i >= max_iterations){
    return FAIL;
  }

  return SUCCESS;
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

  t = (x1 - x1a[i]) / (x1a[i + 1] - x1a[i]);

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


