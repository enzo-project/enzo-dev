/*--------------------------------------------------------------------------
 * StarParticlePhotoelectricHeating
 * ------------------------------------------------------------------------
 * Author : A. Emerick
 * Date   : May 2016
 *
 * For each star particle with an available photoelectric heating flux model,
 * compute the FUV flux from each star in each grid cell and compute the resulting
 * local photoelectric heating.
 *
 * Currently only works for individual star particles
 *---------------------------------------------------------------------------------
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int FindField(int f, int farray[], int n);

int IndividualStarComputeFUVLuminosity(float &L_fuv, const float &mp, const float &metallicity);


/* internal prototypes */
float ComputeHeatingRateFromDustModel(const float &n_H, const float &n_e, const float &Z,
                                      const float &T, const float &G);


int StarParticlePhotoelectricHeating(TopGridData *MetaData,
                                     LevelHierarchyEntry *LevelArray[], int level,
                                     Star* &AllStars){
/* ---------------------------------------------------------------------
 * StarParticlePhotoelectricHeating
 * --------------------------------------------------------------------
 * Author : A. Emerick
 * Date    : May 2016
 * Function prepares for the calculation of the local FUV heating field.
 * First checks if there are FUV stars around. Then loops over all
 * and computes their FUV properties. Then hands this off to each grid
 * to compute the local FUV heating rate in each grid.
 * --------------------------------------------------------------------*/

  if(!STARMAKE_METHOD(INDIVIDUAL_STAR) || !IndividualStarFUVHeating){
    return SUCCESS; // do nothing
  }

  LevelHierarchyEntry *Temp;

  // loop through all grids and set FUV to zero
  // this is done to ensure no heating if there are no more stars left or none formed
  for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){
    for ( Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel){
      if (! Temp->GridData->isLocal() ){ // only do this for local grids
        continue;
      }

      Temp->GridData->ZeroPhotoelectricHeatingField();
    }
  }


  /* quit if there are no stars */
  if(AllStars == NULL){
    return SUCCESS;
  }

  Star *cstar;

  float *Ls, *xs, *ys, *zs;

  int count = 0, number_of_fuv_stars = 0;
  /* figure out how many stars there are that contribute to FUV heating */
  for (cstar = AllStars; cstar; cstar = cstar->NextStar){
    if (cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR && cstar->ReturnBirthMass() > IndividualStarFUVMinimumMass){
      count++;
    }
  }
  number_of_fuv_stars = count;

  if (number_of_fuv_stars == 0){ // return if there are none
    return SUCCESS;
  }

  /* now, do the FUV heating for all stars */
  Ls = new float[number_of_fuv_stars];
  xs = new float[number_of_fuv_stars];
  ys = new float[number_of_fuv_stars];
  zs = new float[number_of_fuv_stars];

  count = 0;
  float *star_pos;
  /* loop over all stars and compute their FUV luminosity and position - store to be passed to grid */
  for (cstar = AllStars; cstar; cstar = cstar->NextStar){
    if (cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR && cstar->ReturnBirthMass() > IndividualStarFUVMinimumMass){

      float L;
      IndividualStarComputeFUVLuminosity(L, cstar->ReturnBirthMass(), cstar->ReturnMetallicity());

      star_pos = cstar->ReturnPosition();

      Ls[count] = L;
      xs[count] = star_pos[0];
      ys[count] = star_pos[1];
      zs[count] = star_pos[2];

      count++;
    } // if star massive enough
  } // end loop over all stars

  /* loop over each grid passing the star info - heating is assigned on each grid */
  for (int l = level; l < MAX_DEPTH_OF_HIERARCHY; l ++){
    for (Temp = LevelArray[l]; Temp; Temp = Temp ->NextGridThisLevel){

      if (!Temp->GridData->isLocal()){ continue ;} // skip if grid is on another proc

      Temp->GridData->AddPhotoelectricHeatingFromStar(Ls, xs, ys, zs, number_of_fuv_stars);

    }
  } // end level loop

  delete[] Ls;
  delete[] xs;
  delete[] ys;
  delete[] zs;

  return SUCCESS;
}

void grid::AddPhotoelectricHeatingFromStar(const float *Ls, const float *xs, const float *ys, const float *zs,
                                           const int &number_of_fuv_stars){
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
  const double pi  = 3.1415621;
  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, TimeUnits, EnergyUnits, MassUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time);
  MassUnits   = DensityUnits * POW(LengthUnits,3);
  EnergyUnits = DensityUnits * VelocityUnits * VelocityUnits;

  /* get temperature field */
  float *temperature;
  int size = 1;
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= this->GridDimension[dim];
  }

  temperature = new float[size];

  if(  this->ComputeTemperatureField(temperature) == FAIL ){
    ENZO_FAIL("Error in compute temperature called from PhotoelectricHeatingFromStar");
  }

  /* get multispecies fields */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, ElectronNum, PeNum, DensNum, MetalNum;
  if (MultiSpecies){
    this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                                HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  } else{
    ENZO_FAIL("StarParticlePhotoelectricHeating: MultiSpeices is required for photoelectrci heating");
  }

  // find fields for density, metal density, electron density, and the heating rate
  DensNum     = FindField(Density, this->FieldType, this->NumberOfBaryonFields);
  MetalNum    = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);
  ElectronNum = FindField(ElectronDensity, this->FieldType, this->NumberOfBaryonFields);
  PeNum       = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);


  /* loop over every cell, sum flux contribution from each star in each cell */
  for(int k = 0; k < this->GridDimension[2]; k++){

    FLOAT zcell = this->CellLeftEdge[2][k] + 0.5*this->CellWidth[2][k];

    for(int j = 0; j < this->GridDimension[1]; j++){
      FLOAT ycell = this->CellLeftEdge[1][j] + 0.5*this->CellWidth[1][j];

      for(int i = 0; i < this->GridDimension[0]; i++){

        int index = i + (j + k * this->GridDimension[1])* this->GridDimension[0];

        FLOAT xcell = this->CellLeftEdge[0][i] + 0.5*this->CellWidth[0][i];

        /* if the cell is below the temperature threshold for dust to exist, apply heating */
        if( temperature[index] < IndividualStarFUVTemperatureCutoff){
          float local_flux = 0.0;
          FLOAT rsqr;

          /* find total local flux due to all stars */
          for (int sp = 0; sp < number_of_fuv_stars; sp++){

            rsqr = (xcell - xs[sp])*(xcell - xs[sp]) +
                   (ycell - ys[sp])*(ycell - ys[sp]) +
                   (zcell - zs[sp])*(zcell - zs[sp]);

            local_flux += Ls[sp] / (4.0 * pi * rsqr * LengthUnits * LengthUnits);

          }


          float n_H, n_e, Z;

          n_H = (this->BaryonField[HINum][index] + this->BaryonField[HIINum][index]);

          if ( MultiSpecies > 1){ /* include H2 */
            n_H += this->BaryonField[HMNum][index] +
                      0.5 * (this->BaryonField[H2INum][index] + this->BaryonField[H2IINum][index]);
          }


          n_H *= DensityUnits / m_h;

          n_e  = this->BaryonField[ElectronNum][index] * DensityUnits / m_e;

          Z    = this->BaryonField[MetalNum][index] / this->BaryonField[DensNum][index]; // metal dens / dens


          // assign heating rate from model
          BaryonField[PeNum][index] = ComputeHeatingRateFromDustModel(n_H, n_e, Z, temperature[index], local_flux);

        } else {
          BaryonField[PeNum][index] = 0.0; // no heating above temperature threshold
        }

      } // end loop over i
    } // j
  } // k


  delete[] temperature;

}

float ComputeHeatingRateFromDustModel(const float &n_H, const float &n_e,
                                      const float &T,   const float &Z,
                                      const float &G){
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

  const double Z_o    =    0.02; // solar metallicity as defined in Forbes et. al. 2016

  float G_o   = G / G_norm;      // local FUV flux normalized to Habing field
  float epsilon, flux;


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

  flux = 1.3E-24 * n_H * epsilon * G_o * (Z / Z_o);

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

  int PeNum;
  PeNum = FindField(PeHeatingRate, this->FieldType, this->NumberOfBaryonFields);

  int size = 1;
  for (int dim = 0; dim < this->GridRank; dim++){
    size *= this->GridDimension[dim];
  }

  for(int i = 0; i < size; i++){
    this->BaryonField[PeNum][i] = 0.0;
  }

}



