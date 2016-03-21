/*--------------------------------------------------------------------------------
/
/ INDIVIDUAL STAR PROPERTIES HELPER ROUTINES
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

float IndividualStarLifetime(float *mp){
  /* ---------------------------------------------------------
   * Compute lifetime given only the mass. This is only used
   * used when radiation for individual stars is not being used
   * as luminosity calculation is not deterministic.
   * value returned in CGS.
   * ---------------------------------------------------------*/
  float luminosity, lifetime;

  luminosity  = IndividualStarLuminosity( mp);
  luminosity /= SOLAR_LUMINOSITY; // convert from cgs to solar

  lifetime   = SOLAR_LIFETIME * (*mp)/(luminosity); // in cgs

  return lifetime;
}

float IndividualStarLuminosity(float *mp){
  /* ------------------------------------------------------------
   * Use the broken power law scaling fits to the M-L
   * relationship in Ekel et. al. 2015 to compute the star's
   * luminosity given the mass. Assume L is Gaussianly distributed
   * about mean relation at fixed mass, using spread taken
   * from Ekel et. al. 2015:
   *
   * log(L [solar]) = alpha * log(M [solar]) + A
   *
   * Luminosity returned in CGS
   * ------------------------------------------------------------- */

  float alpha;     // slope
  float A;         // normalization
  float sigma;     // standard deviation
  float luminosity;

  float rnum;
  int rsign;
  const int max_random = (1<<16);

  // first range is technically 0.38 < M <= 1.05
  if( *mp <= 1.05 ){
    alpha = 4.841;
    A     = -0.026;
    sigma = 0.121;
  } else if (1.05 < *mp && *mp <= 2.40) {
    alpha = 4.328;
    A     = -0.002;
    sigma = 0.108;
  } else if (2.40 < *mp && *mp <= 7.00) {
    alpha = 3.962;
    A     = 0.120;
    sigma = 0.165;
  } else{ // last range is fit over 7 < M <= 32 but extend it
    alpha = 2.726;
    A     = 1.237;
    sigma = 0.158;
  }

  // calculate mean luminosity at given mass
  luminosity = alpha * log10(*mp) + A; // log luminosity

  // now assume gaussian distribution about the mean
  rnum  = (float) (random() % max_random) / (float) (max_random);
  rsign = rnum>0.5 ? 1:-1;
  luminosity = luminosity + rsign * GaussianRandomVariable() * sigma;
  luminosity = POW(10.0, luminosity) * SOLAR_LUMINOSITY;

  return luminosity;
}

float IndividualStarLuminosity(float *mp, float *lifetime){
  /* -----------------------------------------------------
   * Since lifetime is computed using the luminosity, return
   * luminosity given mass and lifetime to avoid having to store
   * the value for luminosity in the simulation
   *
   * M in solar, lifetime in cgs
   * ------------------------------------------------------ */

  return SOLAR_LUMINOSITY * (*mp)/ ( *lifetime / SOLAR_LIFETIME);
}

float IndividualStarRadius(float *mp){
  // computes the stellar radius using mass alone
  // The assumption here is that stellar radius has a one-to-one
  // relationship with mass and does not depend on metallicity
  // The radii are computed assuming a polytropic equation of state
  // scaled to the solar mass / radius. Stars are differentiated
  // based on type of nuclear burning dominance and the associated
  // polytropic value for their E.O.S.... more massive stars are 
  // closer to linear relationship

  float R;

  // assume a polytropic equation of state with p-p dominated burning
  // - in this case n = 4 and R scales as M^( (n-1) / (n+3))
  // this may underestimate radius for most massive stars
  if (*mp <= 0.5){
    R = POW(*mp, 0.9); // fully convective - nearly one-to-one scaling
  } else if (*mp < 2.0){
    R = POW(*mp, 3.0/7.0); // pp chain dominated - n = 4
  } else{
    R = POW(*mp, 15.0/19.0); // CNO dominated - n = 16
  }

  return R * SOLAR_RADIUS;
}

float IndividualStarTeff(float *mp, float *lifetime){
 // given mass and tau, calculate L and R
 // use stefan boltzman to get teff

  float L, R;
  L = IndividualStarLuminosity( mp , lifetime) / SOLAR_LUMINOSITY;
  R = IndividualStarRadius( mp ) / SOLAR_RADIUS;

  return SOLAR_TEFF * POW(L / (R*R), 0.25);
}

float IndividualStarTeff(float *mp, float *L, float *R){
  /* Effective temperature given L and R */

  return SOLAR_TEFF * POW(
                              (*L/SOLAR_LUMINOSITY)
                              / ( (*R)*(*R) /(SOLAR_RADIUS*SOLAR_RADIUS)), 0.25);
}

float IndividualStarSurfaceGravity(float *mp){
  /* given mass, compute radius and get surface gravity (cgs) */

  float R;
  const double G = 6.6743E-8;

  R = IndividualStarRadius( mp );

  return G*(*mp * SOLAR_MASS)/(R*R);
}

float IndividualStarSurfaceGravity(float *mp, float *R){
  /* Compute surface gravity given radius */
  const double G = 6.6743E-8;

  return G*(*mp)/((*R)*(*R));
}

int IndividualStarComputeIonizingRates(float *q0, float *q1, float *Teff,
                                       float *g, float *metallicity){
 /*============================================================
  * IndividualStarComputeIonizingRates
  * ===========================================================
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

  x = (E_ion_HI / eV_erg) / (k_boltz * (*Teff));
  // compute the black body radiance in unitless numbers
  if( (PhotonRadianceBlackBody(q0, x)) == FAIL){
    ENZO_FAIL("IndividualStarComputeIonizingRates: Summation of black body integral failed to converge\n");
  }

  x = (E_ion_He / eV_erg) / (k_boltz * (*Teff));
  if ( PhotonRadianceBlackBody(q1, x) == FAIL ){
    ENZO_FAIL("IndividualStarComputeIonizingRates: Summation of black body integral failed to converge\n");
  }

  float A = 2.0 * k_boltz * k_boltz * k_boltz * (*Teff)*(*Teff)*(*Teff) /
                               (h*h*h*c*c);

   // convert to units of # / s / m-2 / sr-1
  *q0 = A * (*q0);
  *q1 = A * (*q1);

  return SUCCESS;
}

int IndividualStarInterpolateRadData(float *q0, float *q1,
                                     float *Teff, float *g, float *metallicity){
  /* ===================================================================
   * IndividualStarInterpolateRadData
   * ===================================================================
   * Performs trillinear interpolation over effective temperature,
   * surface gravity, and metallicity to compute stellar ionizing fluxes
   * q0 and q1. Method follows billinear interpolation outlined in
   * Numerical Recipes, effectively three first order linear interpolations
   *
   *
   * ===================================================================
   */

  // do a bisect search in each dimension to find index
  int i, j, k; // i -> Teff , j -> g, k -> Z

  int width, bin_number;

  // convert metallicity to solar
  float Z; // Z is in units of solar
  Z = (*metallicity) / SOLAR_METALLICITY;

  // check star values against minimum and maximum tabulated bins
  // FOr now just throw an error. But in the future, if Teff and g are below
  // the gridded values, then just use a black body spectrum as Josh Wall suggested
  // If T or g are ABOVE the gridded table... well... then I have some thinking to do

  if ( *Teff < IndividualStarRadData.T[0] ||
       *Teff > IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins -1]){
    printf("WARNING in IndividualStarInterpolateRadData: Temperature out of bounds - Using black body");
    printf("Teff = %"ESYM" for minimum value %"ESYM" and maximum %"ESYM"\n", *Teff, IndividualStarRadData.T[0], IndividualStarRadData.T[IndividualStarRadData.NumberOfTemperatureBins-1]);
    return -1;
  }
  if ( *g < IndividualStarRadData.g[0] ||
       *g > IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins -1]){
    printf("WARNING in IndividualStarInterpolateRadData: Surface gravity out of bounds - Using black body");
    printf("logg = %"ESYM" for logged minimum value %"ESYM" and maximum %"ESYM"\n", log10(*g), log10(IndividualStarRadData.g[0]), log10(IndividualStarRadData.g[IndividualStarRadData.NumberOfSGBins-1]));
    return -1;
  }
  if ( Z < IndividualStarRadData.Z[0] ||
       Z > IndividualStarRadData.Z[IndividualStarRadData.NumberOfMetallicityBins - 1]){
    printf("WARNING in IndividualStarInterpolateRadData: Metallicity out of bounds - Using black body");
    printf("Z (solar) = %"ESYM" for minimum value %"ESYM" and maximum %"ESYM"\n", Z, IndividualStarRadData.Z[0], IndividualStarRadData.Z[IndividualStarRadData.NumberOfMetallicityBins-1]);
    return -1;
  }

  // binary search over temperature (i)
  width = IndividualStarRadData.NumberOfTemperatureBins / 2;
  i     = IndividualStarRadData.NumberOfTemperatureBins / 2;

  while (width > 1){
    width /=2;
    if ( *Teff > IndividualStarRadData.T[i])
      i += width;
    else if ( *Teff < IndividualStarRadData.T[i])
      i -= width;
    else
      break;
  } // temperature binary search

  // search finds nearest bin - interpolation requires floored nearest bin
  if ( *Teff < IndividualStarRadData.T[i] ) i--;

  // binary search over surface gravity (j)
  width = IndividualStarRadData.NumberOfSGBins / 2;
  j     = IndividualStarRadData.NumberOfSGBins / 2;

  while (width > 1){
    width /= 2;
    if ( *g > IndividualStarRadData.g[j])
      j += width;
    else if ( *g < IndividualStarRadData.g[j])
      j -= width;
    else
      break;
  } // surface gravity binary search

  // search finds closest bin - interpolation requires floored nearest bin
  if ( *g < IndividualStarRadData.g[j]) j--;

  // binary search over metallicity
  width = IndividualStarRadData.NumberOfMetallicityBins / 2;
  k     = IndividualStarRadData.NumberOfMetallicityBins / 2;

  while (width > 1){
    width /=2;
    if ( Z > IndividualStarRadData.Z[k] )
      k += width;
    else if ( Z < IndividualStarRadData.Z[k] )
      k -= width;
    else
      break;
  } // metallicity binary search

  // search finds closest bin - interpolation requires floored nearest bin
  if (Z < IndividualStarRadData.Z[k]) k--;

  /* Now that we've located the point in 3D, compute coefficients */
  float t, u, v;

  t = (*Teff - IndividualStarRadData.T[i]) /
                  (IndividualStarRadData.T[i+1] - IndividualStarRadData.T[i]);
  u = (*g - IndividualStarRadData.g[j]) /
                  (IndividualStarRadData.g[j+1] - IndividualStarRadData.g[j]);
  v = ( Z - IndividualStarRadData.Z[k]) /
                  (IndividualStarRadData.Z[k+1] - IndividualStarRadData.Z[k]);


  /* Since the tabulated data is not fully sampled in surface gravity at every
     temperature, need to make sure that interpolation is not occuring outside
     of the tabulated range. Throw a very loud warning if this is happening    */
  if( IndividualStarRadData.q0[i  ][j  ][k  ] == 1.0 ||
      IndividualStarRadData.q0[i  ][j+1][k  ] == 1.0 ||
      IndividualStarRadData.q0[i+1][j+1][k  ] == 1.0 ||
      IndividualStarRadData.q0[i+1][j  ][k  ] == 1.0 ||
      IndividualStarRadData.q0[i  ][j  ][k+1] == 1.0 ||
      IndividualStarRadData.q0[i  ][j+1][k+1] == 1.0 ||
      IndividualStarRadData.q0[i+1][j+1][k+1] == 1.0 ||
      IndividualStarRadData.q0[i+1][j  ][k+1] == 1.0   ){

    printf("WARNING in IndividaulStarProperties: Values for interpolation outside range:"
           " Teff = %"ESYM" logg = %"ESYM" Z = %"ESYM"\n", *Teff, log10(*g), Z);
    printf("     Nearest neighbor Teff = %"ESYM" logg = %"ESYM" Z = %"ESYM"\n", IndividualStarRadData.T[i], log10(IndividualStarRadData.g[j]), IndividualStarRadData.Z[k]);
    printf(" Using black body instead\n");
    return -1;
  }

  /* Now apply the coefficients and compute the ionizing rate */
  *q0 = (1.0 - t)*(1.0 - u)*(1.0 - v) * IndividualStarRadData.q0[i  ][j  ][k  ] +
        (1.0 - t)*(      u)*(1.0 - v) * IndividualStarRadData.q0[i  ][j+1][k  ] +
        (      t)*(      u)*(1.0 - v) * IndividualStarRadData.q0[i+1][j+1][k  ] +
        (      t)*(1.0 - u)*(1.0 - v) * IndividualStarRadData.q0[i+1][j  ][k  ] +
        (1.0 - t)*(1.0 - u)*(      v) * IndividualStarRadData.q0[i  ][j  ][k+1] +
        (1.0 - t)*(      u)*(      v) * IndividualStarRadData.q0[i  ][j+1][k+1] +
        (      t)*(      u)*(      v) * IndividualStarRadData.q0[i+1][j+1][k+1] +
        (      t)*(1.0 - u)*(      v) * IndividualStarRadData.q0[i+1][j  ][k+1] ;

  /* q1 ionizing flux */
  *q1 = (1.0 - t)*(1.0 - u)*(1.0 - v) * IndividualStarRadData.q1[i  ][j  ][k  ] +
        (1.0 - t)*(      u)*(1.0 - v) * IndividualStarRadData.q1[i  ][j+1][k  ] +
        (      t)*(      u)*(1.0 - v) * IndividualStarRadData.q1[i+1][j+1][k  ] +
        (      t)*(1.0 - u)*(1.0 - v) * IndividualStarRadData.q1[i+1][j  ][k  ] +
        (1.0 - t)*(1.0 - u)*(      v) * IndividualStarRadData.q1[i  ][j  ][k+1] +
        (1.0 - t)*(      u)*(      v) * IndividualStarRadData.q1[i  ][j+1][k+1] +
        (      t)*(      u)*(      v) * IndividualStarRadData.q1[i+1][j+1][k+1] +
        (      t)*(1.0 - u)*(      v) * IndividualStarRadData.q1[i+1][j  ][k+1];


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

int PhotonRadianceBlackBody(float *q, float x){
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

  *q = sum;

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

/*-----------------------------------------------------------------------------
 GaussianRandomVariable

 Returns a random variable selected over a Gaussian distribution with mean
 of zero and varince of unity using the Box-Muller transform
-----------------------------------------------------------------------------*/
float GaussianRandomVariable(void){

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


