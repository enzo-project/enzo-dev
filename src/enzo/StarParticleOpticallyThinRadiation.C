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

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int FindField(int f, int farray[], int n);

/* internal prototypes */
float ComputeHeatingRateFromDustModel(const float &n_H, const float &n_e, 
                                      // const float &T,
                                      const float &Z, const float &G, const float &dx);


int StarParticleOpticallyThinRadiation(TopGridData *MetaData,
                                       LevelHierarchyEntry *LevelArray[], int level,
                                       Star* &AllStars){
/* ---------------------------------------------------------------------
 * StarParticleOpticallyThinRadiation
 * --------------------------------------------------------------------
 * Author : A. Emerick
 * Date    : May 2016
 *           Oct 2016
 * Function prepares for computation of local FUV and local LW flux
 * from each stellar component.
 * --------------------------------------------------------------------*/

  if(!STARMAKE_METHOD(INDIVIDUAL_STAR) ||
    !(IndividualStarFUVHeating || IndividualStarLWRadiation)){
    return SUCCESS; // do nothing
  }

  if (IndividualStarOTRadiationMethod != 0) return SUCCESS; // using RT framework

  LevelHierarchyEntry *Temp;

  // loop through all grids and set FUV to zero
  // this is done to ensure no heating if there are no more stars left or none formed
  for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){
    for ( Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel){
      if (! Temp->GridData->isLocal() ){ // only do this for local grids
        continue;
      }
      Temp->GridData->ZeroPhotoelectricHeatingField();
      Temp->GridData->ZeroOTLWRadiationField();
    }
  }


  /* quit if there are no stars */
  if(AllStars == NULL){
    return SUCCESS;
  }

  Star *cstar;

  float *L_fuv, *L_lw, *xs, *ys, *zs, *ts;

  int count = 0, number_of_ot_stars = 0;
  /* figure out how many stars there are that contribute to optically thin radiation */
  for (cstar = AllStars; cstar; cstar = cstar->NextStar){
    if (cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR &&
        cstar->ReturnBirthMass() >= IndividualStarOTRadiationMass){
      count++;
    }
  }
  number_of_ot_stars = count;

  if (number_of_ot_stars == 0){ // return if there are none
    return SUCCESS;
  }

  TIMER_START("StarParticlePhotoelectricHeating");
  /* now, do the FUV heating for all stars */
  L_fuv = new float[number_of_ot_stars];
  L_lw  = new float[number_of_ot_stars];
  xs = new float[number_of_ot_stars];
  ys = new float[number_of_ot_stars];
  zs = new float[number_of_ot_stars];
  ts = new float[number_of_ot_stars];

  count = 0;
  float *star_pos;
  /* loop over all stars and compute their luminosity and position - store to be passed to grid */
  for (cstar = AllStars; cstar; cstar = cstar->NextStar){
    if (cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR &&
        cstar->ReturnBirthMass() >= IndividualStarOTRadiationMass){

      float fuv_luminosity = 0.0, lw_luminosity = 0.0;

      if(IndividualStarFUVHeating)
          cstar->ComputeFUVLuminosity(fuv_luminosity);

      if(IndividualStarLWRadiation)
          cstar->ComputeLWLuminosity(lw_luminosity);

      star_pos = cstar->ReturnPosition();

      L_fuv[count] = fuv_luminosity;
      L_lw[count]  = lw_luminosity;
      xs[count]    = star_pos[0];
      ys[count]    = star_pos[1];
      zs[count]    = star_pos[2];
      ts[count]    = cstar->ReturnBirthTime();

      count++;
    } // if star massive enough
  } // end loop over all stars

  /* loop over each grid passing the star info - heating is assigned on each grid */
  for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){
    for (Temp = LevelArray[l]; Temp; Temp = Temp ->NextGridThisLevel){

      if (!Temp->GridData->isLocal()){ continue ;} // skip if grid is on another proc

      Temp->GridData->AddOpticallyThinRadiationFromStar(L_fuv, L_lw,
                                                        xs, ys, zs, ts, number_of_ot_stars);

    }
  } // end level loop

  delete[] L_fuv;
  delete[] L_lw;
  delete[] xs;
  delete[] ys;
  delete[] zs;
  delete[] ts;

  TIMER_STOP("StarParticlePhotoelectricHeating");
  return SUCCESS;
}

void grid::AddOpticallyThinRadiationFromStar(const float *L_fuv, const float *L_lw,
                                             const float *xs, const float *ys, const float *zs,
                                             const float *ts, const int &max_number_of_ot_stars){
  /* ---------------------------------------------------------------------------------------------------------
   * AddPhotoelectricHeatingFromStar
   * ---------------------------------------------------------------------------------------------------------
   * Author : A. Emerick
   * Date    : May 2016
   *
   * Given a list of star luminosities and positions, computes the local FUV flux in each cell and uses a dust
   * model to compute the local FUV heating rate in each cell. Stored as a baryon field.
   * ---------------------------------------------------------------------------------------------------------
   */

  const double m_e = 9.109E-28; // in g
  const double m_h = 1.673E-24; // in g
  const double pi  = 3.14159265358979;
  const double c_light = 2.99792458E10; // cgs
  const double H2ISigma = 3.71e-18;
  const double LW_energy = 12.8 / 6.241509E11; // LW band energy in erg
  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, TimeUnits, EnergyUnits, MassUnits;

  float dx = this->CellWidth[0][0];

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time);
  MassUnits   = DensityUnits * POW(LengthUnits,3);
  EnergyUnits = DensityUnits * VelocityUnits * VelocityUnits;

  /* get temperature field */
  int size = 1;
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= this->GridDimension[dim];
  }

//  Current model doesn't need temperature field - 
//  float *temperature;
//  temperature = new float[size];
//  if(  this->ComputeTemperatureField(temperature) == FAIL ){
//    ENZO_FAIL("Error in compute temperature called from PhotoelectricHeatingFromStar");
//  }


  /* get multispecies fields */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, ElectronNum, PeNum, DensNum, MetalNum, OTLWkdissH2INum;
  if (MultiSpecies){
    this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                                HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  } else{
    ENZO_FAIL("StarParticleOpticallyThinRadiation: MultiSpeices is required for photoelectrci heating");
  }

  // find fields for density, metal density, electron density, and the heating rate
  DensNum     = FindField(Density, this->FieldType, this->NumberOfBaryonFields);
  MetalNum    = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);
  ElectronNum = FindField(ElectronDensity, this->FieldType, this->NumberOfBaryonFields);
  PeNum       = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);
  OTLWkdissH2INum = FindField(OTLWkdissH2I, this->FieldType, this->NumberOfBaryonFields);

  FLOAT * xstar, *ystar, *zstar;
  float * L_fuv_star, *L_lw_star, *ts_star;
  float fuv_background_flux = 0.0, lw_background_flux = 0.0;
  int number_of_ot_stars;

  xstar = new float[max_number_of_ot_stars]; ystar = new float[max_number_of_ot_stars];
  zstar = new float[max_number_of_ot_stars];
  L_fuv_star = new float [max_number_of_ot_stars];
  L_lw_star = new float[max_number_of_ot_stars];
  ts_star = new float[max_number_of_ot_stars];

  if (IndividualStarApproximateOTRadiation){
     FLOAT xcenter, ycenter, zcenter;
     xcenter = (this->CellLeftEdge[0][0] + this->CellLeftEdge[0][this->GridDimension[0]-1])*0.5 + this->CellWidth[0][0]; // *0.5;
     ycenter = (this->CellLeftEdge[1][0] + this->CellLeftEdge[1][this->GridDimension[1]-1])*0.5 + this->CellWidth[0][0];
     zcenter = (this->CellLeftEdge[2][0] + this->CellLeftEdge[2][this->GridDimension[2]-1])*0.5 + this->CellWidth[0][0];

    /* find out which stars we can approximate, avoiding the NxM comparison */
    int count = 0;
    FLOAT rsqr, rcenter, width; // distance to cell center
    for (int sp = 0; sp < max_number_of_ot_stars; sp++){
      rsqr = (xs[sp] - xcenter)*(xs[sp]-xcenter) +
             (ys[sp] - ycenter)*(ys[sp]-ycenter) + 
             (zs[sp] - zcenter)*(zs[sp]-zcenter);

      // if ratio between distance to cell center and max grid width > 10, we get 20% error
      width = (this->GridDimension[0]*this->CellWidth[0][0])*(this->GridDimension[0]*this->CellWidth[0][0])+
              (this->GridDimension[1]*this->CellWidth[0][0])*(this->GridDimension[1]*this->CellWidth[0][0])+
              (this->GridDimension[2]*this->CellWidth[0][0])*(this->GridDimension[2]*this->CellWidth[0][0]);
      
      float flux_ratio = width / rsqr;

/*
    for (int sp = 0; sp < max_number_of_ot_stars; sp++){
      FLOAT min_rsqr = huge_number, max_rsqr = tiny_number, avg_rsqr = 0.0;
      // check rsqr values at corner points and find max variation over grid 
      for(int k = 0; k < this->GridDimension[2]; k += this->GridDimension[2]-1){
        FLOAT zcell = this->CellLeftEdge[2][k] + 0.5 * this->CellWidth[2][k];

        for(int j =0; j < this->GridDimension[1]; j += this->GridDimension[1]-1){
        FLOAT ycell = this->CellLeftEdge[1][j] + 0.5 * this->CellWidth[1][j];

          for(int i =0; i < this->GridDimension[0]; i += this->GridDimension[0]-1){
            FLOAT xcell = this->CellLeftEdge[0][i] + 0.5 * this->CellWidth[0][i];

            FLOAT rsqr = (xs[sp] - xcell)*(xs[sp]-xcell) +
                         (ys[sp] - ycell)*(ys[sp]-ycell) +
                         (zs[sp] - zcell)*(zs[sp]-zcell);

            if(rsqr < min_rsqr)
              min_rsqr = rsqr;
            if(rsqr > max_rsqr)
              max_rsqr = rsqr;
            avg_rsqr += rsqr;
          }
        }
      } // end z loop
*/    

//      float inv_avg_rsqr = 1.0 / (avg_rsqr / (8.0)); // there are 8 corner points

      // compute flux difference in min and max cell
//      float flux_ratio = fmax( fabs( (1.0/min_rsqr - inv_avg_rsqr)/inv_avg_rsqr), 
//                               fabs( (1.0/max_rsqr- inv_avg_rsqr)/inv_avg_rsqr));

      // make shorter threshold name
      if (flux_ratio <= IndividualStarApproximateOTThreshold){
        /* we can approximate this star - add to background radiation level */

//        FLOAT rsqr = (xs[sp] - xcenter)*(xs[sp] - xcenter) +
//                    (ys[sp] - ycenter)*(ys[sp] - ycenter) +
//                     (zs[sp] - zcenter)*(zs[sp] - zcenter);

        fuv_background_flux += L_fuv[sp] / (4.0 * pi * rsqr * LengthUnits * LengthUnits);
        lw_background_flux  += L_lw[sp]  / (4.0 * pi * rsqr * LengthUnits * LengthUnits);

      } else{
        /* we cannot approximate this star - add to list of stars */
        xstar[count] = xs[sp];
        ystar[count] = ys[sp];
        zstar[count] = zs[sp];
        L_fuv_star[count] = L_fuv[sp];
        L_lw_star[count]  = L_lw[sp];
        ts_star[count]      = ts[sp];
        count++;
      }
    } // end loop over stars

    number_of_ot_stars = count - 1;
  } else{
    /* just use full star arrays - this copy over is not ideal */
    for(int i = 0; i < max_number_of_ot_stars; i++){
      xstar[i] = xs[i];
      ystar[i] = ys[i];
      zstar[i] = zs[i];
      L_fuv_star[i] = L_fuv[i];
      L_lw_star[i]  = L_lw[i];
      ts_star[i] = ts[i];
    }

    number_of_ot_stars = max_number_of_ot_stars;
  }


  /* loop over every cell, sum flux contribution from each star in each cell */
  for(int k = 0; k < this->GridDimension[2]; k++){

    FLOAT zcell = this->CellLeftEdge[2][k] + 0.5*this->CellWidth[2][k];

    for(int j = 0; j < this->GridDimension[1]; j++){
      FLOAT ycell = this->CellLeftEdge[1][j] + 0.5*this->CellWidth[1][j];

      for(int i = 0; i < this->GridDimension[0]; i++){

        int index = i + (j + k * this->GridDimension[1])* this->GridDimension[0];

        FLOAT xcell = this->CellLeftEdge[0][i] + 0.5*this->CellWidth[0][i];

        /* if the cell is below the temperature threshold for dust to exist, apply heating */
        //if( temperature[index] < IndividualStarFUVTemperatureCutoff){
          float local_fuv_flux = fuv_background_flux;
          float local_lw_flux  = lw_background_flux;
          FLOAT rsqr;

          /* find total local flux due to all stars */
          for (int sp = 0; sp < number_of_ot_stars; sp++){

            rsqr = (xcell - xstar[sp])*(xcell - xstar[sp]) +
                   (ycell - ystar[sp])*(ycell - ystar[sp]) +
                   (zcell - zstar[sp])*(zcell - zstar[sp]);
            rsqr = fmax(rsqr, 0.0625*dx*dx); // minimum separation of 1/4 cell to avoid divide by zero issues

            float speed = (sqrt(rsqr) * LengthUnits) / ((this->Time - ts_star[sp]) * TimeUnits);

            if ( speed <= c_light ){
                    local_fuv_flux += L_fuv_star[sp] / (4.0 * pi * rsqr * LengthUnits * LengthUnits);
                    local_lw_flux  += L_lw_star[sp]  / (4.0 * pi * rsqr * LengthUnits * LengthUnits);
            }

          }

          if(IndividualStarFUVHeating){
            float n_H, n_e, Z;

            n_H = (this->BaryonField[HINum][index] + this->BaryonField[HIINum][index]);

            if ( MultiSpecies > 1){ /* include H2 */
              n_H += this->BaryonField[HMNum][index] +
                        0.5 * (this->BaryonField[H2INum][index] + this->BaryonField[H2IINum][index]);
            }


            n_H *= DensityUnits / m_h;

            n_e  = this->BaryonField[ElectronNum][index] * DensityUnits / m_e;

            Z    = this->BaryonField[MetalNum][index] / this->BaryonField[DensNum][index]; // metal dens / dens

            // assign heating rate from model - adds to existing background (if present)
            BaryonField[PeNum][index]  += ComputeHeatingRateFromDustModel(n_H, n_e, 
                                                                        // temperature[index],
                                                                         Z, local_fuv_flux,
                                                                         dx*LengthUnits) / (EnergyUnits/TimeUnits);
//            BaryonField[PeNum][index] /= (EnergyUnits / TimeUnits);
          } // end PE heating

          if(IndividualStarLWRadiation){
            BaryonField[OTLWkdissH2INum][index] = H2ISigma * local_lw_flux * TimeUnits / ( LW_energy) ;
          }

        //} else {
        //  BaryonField[PeNum][index] = 0.0; // no heating above temperature threshold
        //}

      } // end loop over i
    } // j
  } // k


  // delete[] temperature;

//  if (IndividualStarApproximateOTRadiation){
    delete [] xstar;
    delete [] ystar;
    delete [] zstar;
    delete [] L_fuv_star;
    delete [] L_lw_star;
    delete [] ts_star;
//  }

//  delete[] approximate_radiation;
}

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
   * and the FUV flux. This works well for MW like galaxies but its abilitiy to handle
   * low metetallicity galaxies is questionable.
   *
   * CURRENTLY DOES NOT DO METAL SCALIGN - THIS NEEDS TO BE PUT IN
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
    epsilon = 0.01488637246 * POW(n_H, 0.235269059);
  }
//
// AJE - 6/17/16
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

    float attenuation = exp( -1.33E-21 * D * dx * n_H); // dx and n_H in cgs

    flux = 1.3E-24 * n_H * epsilon * G_o * D * attenuation;

  } else{
    ENZO_FAIL("PhotoelectricHeatingDustModel must be either 0 or 1\n");
  }

  return flux;
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

  const float G_background = 0.00324 * 1.59E-3; // HM2012 FUV background at z = 0

  int PeNum;
  PeNum = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);

  if (PeNum < 0) return; // not enabled - no need to zero

  int size = 1;
  for (int dim = 0; dim < this->GridRank; dim++){
    size *= this->GridDimension[dim];
  }


  //if(UseUVBackgroundFUVRate){
  if (TRUE){
    const double m_e = 9.109E-28; // in g
    const double m_h = 1.673E-24; // in g

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

    // find fields for density, metal density, electron density, and the heating rate
    DensNum     = FindField(Density, this->FieldType, this->NumberOfBaryonFields);
    MetalNum    = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);
    ElectronNum = FindField(ElectronDensity, this->FieldType, this->NumberOfBaryonFields);

    float n_H, n_e, Z;

    for( int i = 0; i < size; i++){ // apply background heating rate
      n_H = (this->BaryonField[HINum][i] + this->BaryonField[HIINum][i]);

      if ( MultiSpecies > 1){ /* include H2 */
         n_H += this->BaryonField[HMNum][i] +
           0.5 * (this->BaryonField[H2INum][i] + this->BaryonField[H2IINum][i]);
      }


      n_H *= DensityUnits / m_h;

      n_e  = this->BaryonField[ElectronNum][i] * DensityUnits / m_e;

      Z    = this->BaryonField[MetalNum][i] / this->BaryonField[DensNum][i]; // metal dens / dens

      // assign heating rate from model
      BaryonField[PeNum][i]  = ComputeHeatingRateFromDustModel(n_H, n_e, 
                                                              // 100.0, // temperature doesn't matter
                                                                   Z, G_background,
                                                            (this->CellWidth[0][0])*LengthUnits);
      BaryonField[PeNum][i] /= (EnergyUnits / TimeUnits);

    } // end loop

  } else{

    for(int i = 0; i < size; i++){
      this->BaryonField[PeNum][i] = 0.0;
    }
  }

  return;
}



void grid::ZeroOTLWRadiationField(void){
 /* ----------------------------------------------
  * ZeroOTLWRadiationField
  * ----------------------------------------------
  * A. Emerick - OCt 2016
  *
  * Set H2I dissociation field from optically thin
  * LW to zero - recompute field each timestep
  * -----------------------------------------------*/

  int OTLWkdissH2INum;
  OTLWkdissH2INum = FindField(OTLWkdissH2I, this->FieldType, this->NumberOfBaryonFields);

  if (OTLWkdissH2INum < 0) return; // not enabled - no need to zero

  int size = 1;
  for (int dim = 0; dim < this->GridRank; dim++){
    size *= this->GridDimension[dim];
  }

  for(int i = 0; i < size; i++){
    this->BaryonField[OTLWkdissH2INum][i] = 0.0;
  }

}

