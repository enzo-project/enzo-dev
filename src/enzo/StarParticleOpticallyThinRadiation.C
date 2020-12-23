/*--------------------------------------------------------------------------
 * StarParticleOpticallyThinRadiation
 * ------------------------------------------------------------------------
 * Author : A. Emerick
 * Date   : May 2016
 *
 * For each star particle, add contributions of tracked optically thin radiation
 * components. Currently, this is either FUV radiation, leading to photoelectric
 * heating, or LW radiation, leading to H2 dissociation. This computes flux from
 * each star in each cell and the assocaited heating / dissociation rates.
 * To avoid full MxN comparison in cells/stars, uses constant flux calculated at
 * grid center when star's flux varies by <= 10% over grid size.
 *
 * Currently only for individual star particles
 *---------------------------------------------------------------------------------
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "CosmologyParameters.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int FindField(int f, int farray[], int n);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int search_lower_bound(float *arr, float value, int low, int high, int total);


/* internal prototypes */
float ComputeHeatingRateFromDustModel(const float &n_H, const float &n_e,
                                      // const float &T,
                                      const float &Z, const float &G, const float &dx);


float NormalizedDustToGasRatio(const float &Z){
  // Broken power law fit to gas to dust ratio
  // from Remy-Ruyer et. al. 2014, using log metallicity
  // in solar as "x" instead of 12 + log(O/H) as is done
  // in that work. The below gives the right results one
  // would expect assuming a one-to-one with 12+log(O/H)
  // and metallicity

  const float x_t     = -0.73; // 7.96 - 8.69 - threshold from paper
  const float Z_solar = 0.014; // as used in Remy-Ruyer
  const float D_mw    = 6.616595E-3; // MW / solar dust to gas ratio, 1 / 10**(2.21)

  float g_to_d, d_to_g, x;

  x = log10( Z / Z_solar);

  // compute log(Gas / Dust)
  if (x <= x_t){ // low metallicity regime
    g_to_d = 0.68 - 3.08 * x;
  } else{
    g_to_d = 2.21 - 1.00 * x;
  }

  d_to_g = 1.0 / POW(10.0, g_to_d); // convert to actual dust to gas

  // return in MW units

  d_to_g /= D_mw;

  return d_to_g;
}


float ComputeHeatingRateFromDustModel(const float &n_H, const float &n_e,
                                      //const float &T,
                                      const float &Z,
                                      const float &G, const float &dx){
  /* ---------------------------------------------------------------------------------
   * ComputeHeatingRateFromDustModel
   * ---------------------------------------------------------------------------------
   * A. Emerick - May 2016
   *
   * Photoelectric heating rate of dust grains from local FUV flux as computed from
   * models given in Wolfire et. al. 2013 (adopted from Bakes & Tielens 1994). This
   * heating rate is a function of hydrogen and electron number density, temperature,
   * and the FUV flux.
   *
   * If PhotoelectricHeatingDustModel = 1, heating will be scaled with local
   * dust model that is linear with metallicity at high Z, and falls sharply
   * with metallicity at low Z. Local attenuation is computed if dx > 0 (dx <= 0
   * will durn off attenuation approximation.)
   * ----------------------------------------------------------------------------------
   */

  const double G_norm = 1.59E-3; // Habing field normalization - MW flux in cgs
                                 // This is 5.29E-14 erg cm^(-3) times the speed of light
                                 // to convert from energy density to flux density
  const double Z_o    =    0.02; // solar metallicity as defined in Forbes et. al. 2016
  float G_o   = G / G_norm;      // local FUV flux normalized to Habing field
  float epsilon, flux;

//  G_o = fmax(G_o, 0.00324); Background flux threshold now handled in zero background rates routine

  if ( PhotoelectricHeatingDustModelEfficiency > 0.0){
    epsilon = PhotoelectricHeatingDustModelEfficiency;
  } else {
    // power law fit to Wolfire 2003 figure 10b - (hard coded for now)
    // where we are having efficiency scale with local gas density
    epsilon = 0.01488637246 * POW(n_H, 0.235269059) / 1.7; // factor of 1.7 to account for draine G = 1.7 vs. G = 1 Habing
  }
//
//     Taking out full model for now. Would need to be able to compute n_e accurately
//     at all n and T, which we cannot do with just primordial chemistry. Turn this
//     back on if / when a model that includes contributions of C, dust, and PAH
//     ionizations to n_e is used.
//
//  } else {
//    epsilon = 4.9E-2                        / (1.0 + 4.0E-3 * POW(( G_o * POW(T,0.5) / n_e ),0.73)) +
//              3.7E-2 * POW(1.0E-4 * T, 0.7) / (1.0 + 2.0E-4 *    (( G_o * POW(T,0.5) / n_e )     ));
// }

  if (PhotoelectricHeatingDustModel == 0){
    // simplest approximation possible, scaling linerly with metallicity
    // and no self-shielding approximation.
    // Linear scaling valid for Z > 0.1 Z_sun, and ignoring shielding
    // is O.K. for low density / low dust regimes (aka. be careful)

      flux = 1.3E-24 * n_H * epsilon * G_o * (Z / Z_o);
  } else if (PhotoelectricHeatingDustModel == 1){
    // do a slightly better job in high density / high metallicity
    // regimes and do a very approximate, local self-shielding, while
    // using a broken power law dust to gas ratio to better model
    // the different scaling with dust and metallicity below 0.1 Z_sun

    // following Hu et. al. 2016, and Bergin et. al. 2004
    float D = NormalizedDustToGasRatio(Z);

    float attenuation = exp( -1.33E-21 * D * max(dx,0.0) * n_H); // dx and n_H in cgs

    flux = 1.3E-24 * n_H * epsilon * G_o * D * attenuation;

  } else{
    ENZO_FAIL("PhotoelectricHeatingDustModel must be either 0 or 1\n");
  }

  return flux;
}


void grid::ComputeBackgroundFUV(float &G_background){

  const int num_z_bins = 25;
  const float G_background_HM2012[num_z_bins] = { 3.15769E-6, 5.72059E-6, 1.02168E-5, 1.76486E-5,
                                      2.93557E-5, 4.71664E-5, 7.27080E-5, 1.07841E-4,
                                      1.52265E-4, 2.05013E-4, 2.63426E-4, 3.06212E-4,
                                      3.24639E-4, 3.19161E-4, 2.97805E-4, 2.69064E-4,
                                      2.37392E-4, 2.06384E-4, 1.79090E-4, 1.54363E-4,
                                      1.31442E-4, 1.09434E-4, 8.71460E-5, 6.04930E-5,
                                      2.42234E-5};

  const float z_background_HM2012[num_z_bins] = { 0.00000E0 , 1.22020E-1, 2.58930E-1, 4.42540E-1,
                                      5.84890E-1, 7.78280E-1, 9.95260E-1, 1.23870E0 ,
                                      1.51190E0 , 1.81840E0 , 2.16230E0 , 2.54810E0 ,
                                      2.98110E0 , 3.46680E0 , 4.01190E0 , 4.62340E0 ,
                                      5.30960E0 , 6.07950E0 , 6.94330E0 , 7.91250E0 ,
                                      9.00000E0 , 1.02200E1 , 1.15890E1 , 1.31250E1 ,
                                      1.48490E1};

  FLOAT CurrentRedshift = RadiationFieldRedshift;

  if (ComovingCoordinates){
    FLOAT a, dadt, aUnits;

    /* Compute the current redshift (for information only). */
    CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
    aUnits = 1.0 / (1.0 + InitialRedshift);
    CurrentRedshift = (1.0)/(a*aUnits) - 1.0;
  }

  /* Interpolate with current redshift */
  if (!(UseFUVBackground)){
    G_background = 0.0;
  } else if (CurrentRedshift > z_background_HM2012[num_z_bins-1]){
    G_background = 0.0; // turn it off at very high z... no dust anyway
  } else if ( CurrentRedshift <= z_background_HM2012[0]) {
    G_background = G_background_HM2012[0]; // leave at z = 0 value
  } else{
      /* linear interpolate */

    int index = search_lower_bound((float*)z_background_HM2012,
                                   CurrentRedshift, 0, num_z_bins, num_z_bins);

    float t = (CurrentRedshift - z_background_HM2012[index]) /
              (z_background_HM2012[index+1] - z_background_HM2012[index]);

    G_background = (1.0 - t) * G_background_HM2012[index] + t *G_background_HM2012[index+1];

  }

  return;
}


void grid::ZeroPhotoelectricHeatingField(void){
 /* ----------------------------------------------
  * ZeroPhotoelectricHeatingField
  * ----------------------------------------------
  * A. Emerick - May 2016
  *
  * Set photoelectric heating field to zero for all cells.
  * Field is recomputed every timestep.
  * -----------------------------------------------*/

//  const float G_background = 0.00324 * 1.59E-3; // HM2012 FUV background at z = 0

  /* Background FUV Flux density from HM2012 UVB integration */
  const int num_z_bins = 25;
  const float G_background_HM2012[num_z_bins] = { 3.15769E-6, 5.72059E-6, 1.02168E-5, 1.76486E-5,
                                      2.93557E-5, 4.71664E-5, 7.27080E-5, 1.07841E-4,
                                      1.52265E-4, 2.05013E-4, 2.63426E-4, 3.06212E-4,
                                      3.24639E-4, 3.19161E-4, 2.97805E-4, 2.69064E-4,
                                      2.37392E-4, 2.06384E-4, 1.79090E-4, 1.54363E-4,
                                      1.31442E-4, 1.09434E-4, 8.71460E-5, 6.04930E-5,
                                      2.42234E-5};

  const float z_background_HM2012[num_z_bins] = { 0.00000E0 , 1.22020E-1, 2.58930E-1, 4.42540E-1,
                                      5.84890E-1, 7.78280E-1, 9.95260E-1, 1.23870E0 ,
                                      1.51190E0 , 1.81840E0 , 2.16230E0 , 2.54810E0 ,
                                      2.98110E0 , 3.46680E0 , 4.01190E0 , 4.62340E0 ,
                                      5.30960E0 , 6.07950E0 , 6.94330E0 , 7.91250E0 ,
                                      9.00000E0 , 1.02200E1 , 1.15890E1 , 1.31250E1 ,
                                      1.48490E1};


  float G_background = 0.0;

  int PeNum;
  PeNum = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);
  const int FUVRateNum = FindField(FUVRate, this->FieldType, this->NumberOfBaryonFields);

  if ((PeNum < 0) && (FUVRateNum < 0)) return; // not enabled - no need to zero

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  //if(UseUVBackgroundFUVRate){
  if (TRUE){

    // ---- get units ------
    float TemperatureUnits, DensityUnits, LengthUnits,
          VelocityUnits, TimeUnits, EnergyUnits, MassUnits;

    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
             &TimeUnits, &VelocityUnits, this->Time);
    EnergyUnits = DensityUnits * VelocityUnits * VelocityUnits;

    // ---- get field numbers -----
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
        DINum, DIINum, HDINum, ElectronNum, DensNum, MetalNum;
    if (MultiSpecies){
      this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                                  HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

    } else{
      ENZO_FAIL("StarParticleOpticallyThinRadiation: MultiSpeices is required for photoelectrci heating");
    }

    FLOAT CurrentRedshift = RadiationFieldRedshift;

    if (ComovingCoordinates){
      FLOAT a, dadt, aUnits;

      /* Compute the current redshift (for information only). */
      CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt);
      aUnits = 1.0 / (1.0 + InitialRedshift);
      CurrentRedshift = (1.0)/(a*aUnits) - 1.0;
    }

    /* Interpolate with current redshift */
    if (!(UseFUVBackground)){
      G_background = 0.0;
    } else if (CurrentRedshift > z_background_HM2012[num_z_bins-1]){
      G_background = 0.0; // turn it off at very high z... no dust anyway
    } else if ( CurrentRedshift <= z_background_HM2012[0]) {
      G_background = G_background_HM2012[0]; // leave at z = 0 value
    } else{
      /* linear interpolate */

      int index = search_lower_bound((float*)z_background_HM2012,
                                     CurrentRedshift, 0, num_z_bins, num_z_bins);

      float t = (CurrentRedshift - z_background_HM2012[index]) /
                (z_background_HM2012[index+1] - z_background_HM2012[index]);

      G_background = (1.0 - t) * G_background_HM2012[index] + t *G_background_HM2012[index+1];

    }

    // find fields for density, metal density, electron density, and the heating rate
    DensNum     = FindField(Density, this->FieldType, this->NumberOfBaryonFields);
    MetalNum    = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);
    ElectronNum = FindField(ElectronDensity, this->FieldType, this->NumberOfBaryonFields);


    /* get temperature field */
    float *temperature;
    temperature = new float[size];
    if(  this->ComputeTemperatureField(temperature) == FAIL ){
      ENZO_FAIL("Error in compute temperature called from PhotoelectricHeatingFromStar");
    }



    float n_H, n_e, Z;


    // Note inconsistency here (same as defined elsewhere, but different form. 
    //      need to fix this....
    const float FluxConv = EnergyUnits / TimeUnits * LengthUnits;
    const float FluxConv_inv = 1.0 / FluxConv;


    if (G_background > 0.0){
      for (int k = 0; k < GridDimension[2]; k++){
        for(int j = 0; j < GridDimension[1]; j++){
          int index = (k*GridDimension[1] + j)*GridDimension[0];
          for (int i = 0; i < GridDimension[0]; i++, index++){



              n_H = (this->BaryonField[HINum][index] + this->BaryonField[HIINum][index]);

              if ( MultiSpecies > 1){ /* include H2 */
                n_H += this->BaryonField[HMNum][index] +
               0.5 * (this->BaryonField[H2INum][index] + this->BaryonField[H2IINum][index]);
              }
              n_H *= DensityUnits / mh;

              n_e  = this->BaryonField[ElectronNum][index] * DensityUnits / me;
              Z    = this->BaryonField[MetalNum][index] / this->BaryonField[DensNum][index]; // metal dens / dens

              if (FUVRateNum > 0) BaryonField[FUVRateNum][index] = G_background * FluxConv_inv;

              // assign heating rate from model
              if (temperature[index] > IndividualStarFUVTemperatureCutoff){
                BaryonField[PeNum][index]  = 0.0;

              } else {
                BaryonField[PeNum][index]  = ComputeHeatingRateFromDustModel(n_H, n_e,
                                                                  // 100.0, // temperature doesn't matter
                                                                       Z, G_background,
                                                                (this->CellWidth[0][0])*LengthUnits);
                BaryonField[PeNum][index] /= (EnergyUnits / TimeUnits);
              }
          } // end k
        }
      } // end i
    } else { /* background is zero */
      for (int k = 0; k < GridDimension[2]; k++){
        for(int j = 0; j < GridDimension[1]; j++){
          int index = (k*GridDimension[1] + j)*GridDimension[0];
          for (int i = 0; i < GridDimension[0]; i++, index++){
            BaryonField[PeNum][index] = 0.0;
            if (FUVRateNum > 0) BaryonField[FUVRateNum][index] = 0.0;
          }
        }
      }
    }

    delete [] temperature;

  } else{
      for (int k = 0; k < GridDimension[2]; k++){
        for(int j = 0; j < GridDimension[1]; j++){
          int index = (k*GridDimension[1] + j)*GridDimension[0];
          for (int i = 0; i < GridDimension[0]; i++, index++){
            BaryonField[PeNum][index] = 0.0;
            if (FUVRateNum > 0) BaryonField[FUVRateNum][index] = 0.0;
          }
        }
      }

  }

  return;
}
