#include <stdlib.h>
#include <stdio.h>
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

#include "StellarYieldsRoutines.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);


int IndividualStarGetSETablePosition (int &i, int &j, const float &M, const float &metallicity);
int IndividualStarGetRadTablePosition(int &i, int &j, int &k,
                                      const float &Teff, const float &g, const float &metallicity);
int StellarYieldsGetYieldTablePosition(int &i, int &j, const float &M, const float &metallicity);


float IndividualStarSurfaceGravity(const float &mp, const float &R);

void IndividualStarInterpolateProperties(float &Teff, float &R,
                                         const int &i, const int &j,
                                         const float &M, const float &metallicity);

// need to put functions here:
// IndividualStar etc.
float IndividualStarSurfaceGravity(const float &mp, const float &R);
void IndividualStarComputeIonizingRates(float &q0, float &q1,
                                        const int &i, const int &j, const int &k,
                                        const float &Teff, const float &g,
                                        const float &metallicity);
int  IndividualStarComputeLWLuminosity(float &L_LW, const int &i,
                                       const int &j, const int &k,
                                       const float &Teff, const float &R,
                                       const float &g, const float &Z);
int  IndividualStarComputeFUVLuminosity(float &L_FUV, const int &i,
                                       const int &j, const int &k,
                                       const float &Teff, const float &R,
                                       const float &g, const float &Z);
int  IndividualStarInterpolateLifetime(float & tau, const int &i,
                                       const int &j, const float &M,
                                       const float &metallicity, const int &mode);
float StellarYieldsInterpolateYield(int yield_type, const int &i,
                                    const int &j ,const float &M,
                                    const float &metallicity,
                                    int atomic_number);


/* Here I have overloaded all of the interpolation functions for the individual
   star routines in order to be a little more self-consistent. Each routine
   now 100% handles the interpolation without having to pass parameters,
   while checking to make sure the table positions are already computed for
   each star. Doing it this way ensures that the table positions will always
   be set when needed, that they are only interpolated when actually needed,
   and that this is only ever done ONCE for each star in each timestep
*/

void Star::InterpolateProperties(void){

  // set this if not set already
  this->AssertInterpolationPositions(1);

  IndividualStarInterpolateProperties(this->Teff, this->Radius,
                                      this->se_table_position[0],
                                      this->se_table_position[1],
                                      this->BirthMass, this->Metallicity);

  this->SurfaceGravity = IndividualStarSurfaceGravity(this->BirthMass,
                                                      this->Radius);
  return;
}

void Star::AssertStellarProperties(void){

  if ( (this->Teff <= 0.0) ||
       (this->Radius <= 0.0) ||
       (this->SurfaceGravity <= 0.0)){

         this->InterpolateProperties();

       }
  return;
}

void Star::ComputeIonizingRates(float &Q0, float &Q1){

  this->AssertStellarProperties();       // need Teff, g
  this->AssertInterpolationPositions(2); // need

  IndividualStarComputeIonizingRates(Q0, Q1,
                                     this->rad_table_position[0],
                                     this->rad_table_position[1],
                                     this->rad_table_position[2],
                                     this->Teff, this->SurfaceGravity,
                                     this->Metallicity);
  return;
}

void Star::ComputeLWLuminosity(float &L_LW){

  this->AssertStellarProperties();
  this->AssertInterpolationPositions(2);

  IndividualStarComputeLWLuminosity(L_LW,
      this->rad_table_position[0], this->rad_table_position[1], this->rad_table_position[2],
      this->Teff, this->Radius, this->SurfaceGravity, this->Metallicity);
  return;
}

void Star::ComputeFUVLuminosity(float &L_FUV){

  this->AssertStellarProperties();
  this->AssertInterpolationPositions(2);

  IndividualStarComputeFUVLuminosity(L_FUV,
      this->rad_table_position[0], this->rad_table_position[1], this->rad_table_position[2],
      this->Teff, this->Radius, this->SurfaceGravity, this->Metallicity);

  return;
}

int Star::InterpolateLifetime(float &tau, const int & mode ){

  this->AssertInterpolationPositions(1);

  return IndividualStarInterpolateLifetime(tau, this->se_table_position[0],
                                    this->se_table_position[1],
                                    this->BirthMass, this->Metallicity,
                                    mode);
}

double Star::ReturnSurfaceGravity(void){
  this->AssertStellarProperties();
  return this->SurfaceGravity;
}

double Star::ReturnEffectiveTemperature(void){
  this->AssertStellarProperties();
  return this->Teff;
}

double Star::ReturnRadius(void){
  this->AssertStellarProperties();
  return this->Radius;
}

float Star::InterpolateYield(int yield_type, int atomic_number){

  this->AssertInterpolationPositions(3);

  if (  ( this->type == PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII) ||
        ((IndividualStarPopIIIFormation) && (this->Metallicity < PopIIIMetalCriticalFraction))){

    return StellarYieldsInterpolatePopIIIYield(this->yield_table_position[0],
                                               this->BirthMass, atomic_number);
  } else{

    return StellarYieldsInterpolateYield(yield_type,
                                         this->yield_table_position[0],
                                         this->yield_table_position[1],
                                         this->BirthMass, this->Metallicity,
                                         atomic_number);
  }
}

///////////////////////////////////////////////////////////////////////////////

void Star::AssignInterpolationTablePositions(void){
  /* Sets table positions for each star */

  if (this->TablePositionsAssigned()){ return;} // don't reassign!!

  if( ABS(this->type) == PARTICLE_TYPE_INDIVIDUAL_STAR ){

    if ( !( this->FeedbackFlag == NO_FEEDBACK ) ||
          ( this->Mass <= IndividualStarSNIIMassCutoff)) {
        this->AssignSETablePosition();
    }

    /* if we have any radiation on, assign radiation properties */
    if( (RadiativeTransfer && (this->BirthMass >= IndividualStarRadiationMinimumMass)) ||
        ((IndividualStarFUVHeating || IndividualStarLWRadiation) && (this->BirthMass >= IndividualStarOTRadiationMass))){

      this->AssignRadTablePosition();
    }

    // set yield table position if following yields and if feedback is ON
    if (IndividualStarFollowStellarYields && !(this->FeedbackFlag == NO_FEEDBACK) ){
      this->AssignYieldTablePosition();
    }
  }

  return;
}

int Star::TablePositionsAssigned(void){
  int result = 0;

  if ( (this->se_table_position[0] >= 0) &&
       (this->se_table_position[1] >= 0) &&
       (this->rad_table_position[0] >= 0) &&
       (this->rad_table_position[1] >= 0) &&
       (this->rad_table_position[2] >= 0) &&
       (this->yield_table_position[0] >= 0) &&
       (this->yield_table_position[1] >= 0)) result = 1;

  return result;
}

void Star::AssignSETablePosition(void){

  //  if (this->se_table_position[0] >= 0 && this->se_table_position[0] >= 0) return;

  /* find stellar properties - used for winds, SN, and radiation */
  if ( (ABS(this->type) == PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII)  ||
       (ABS(this->type) == PARTICLE_TYPE_INDIVIDUAL_STAR_UNRESOLVED) ||
       (IndividualStarPopIIIFormation && (this->Metallicity < PopIIIMetalCriticalFraction))){
    // We do not need these for PopIII stars. Set separately
    // using hard-coded routines in pop3_properties.F and
    // ComputePhotonRates. Set to INT_UNDEFINED so it is > 0
    // and routine doesn't keep trying to set it to a value.

    this->se_table_position[0] = -1;
    this->se_table_position[1] = -1;
  } else {
    IndividualStarGetSETablePosition(this->se_table_position[0],
                                     this->se_table_position[1],
                                     this->BirthMass, this->Metallicity);
  }

  return;
}

void Star::AssignRadTablePosition(void){

  //if (this->rad_table_position[0] >= 0 &&
  //    this->rad_table_position[1] >= 0 &&
  //   this->rad_table_position[2] >= 0){ return;}

  this->AssertInterpolationPositions(1);

  float Teff, R;
  IndividualStarInterpolateProperties(Teff,R,
                                      this->se_table_position[0],
                                      this->se_table_position[1],
                                      this->BirthMass, this->Metallicity);

  float g = IndividualStarSurfaceGravity(this->BirthMass, R);

  IndividualStarGetRadTablePosition(this->rad_table_position[0],
                                    this->rad_table_position[1],
                                    this->rad_table_position[2],
                                    Teff, g, this->Metallicity);


  return;
}

void Star::AssignYieldTablePosition(void){

  //if (this->yield_table_position[0] >= 0 &&
  //    this->yield_table_position[1] >= 0){ return;}

  if ((ABS(this->type) == IndividualStarPopIII)   ||
      ((IndividualStarPopIIIFormation) && (this->Metallicity <= PopIIIMetalCriticalFraction))) {
    StellarYieldsGetPopIIIYieldTablePosition(this->yield_table_position[0], this->BirthMass);
    this->yield_table_position[1] = 0;
  } else {
    StellarYieldsGetYieldTablePosition(this->yield_table_position[0],
                                       this->yield_table_position[1],
                                       this->BirthMass, this->Metallicity);
  }

  return;
}

void Star::AssertInterpolationPositions(void){
  this->AssertInterpolationPositions(1);
  this->AssertInterpolationPositions(2);
  this->AssertInterpolationPositions(3);

  return;
}

void Star::AssertInterpolationPositions(int mode){
 // mode = 1: se    table
 //        2: rad   table
 //        3: yield table

//  if (this->type < 0) return; // do not do anything

  switch(mode){
    case 1:
      if (this->se_table_position[0] < 0 || this->se_table_position[1] < 0){

        if(IndividualStarSaveTablePositions){
          printf("STAR: Table positions not being set correctly. Setting SE table.\n");
        }
        this->AssignSETablePosition();
      }
      break;

    case 2:
      // -9 is used as a flag to use black body calculation instead of table
      if( (this->rad_table_position[0] < 0 && this->rad_table_position[0] > -2) &&
           this->BirthMass > IndividualStarRadiationMinimumMass &&
           this->type == PARTICLE_TYPE_INDIVIDUAL_STAR){  // only interpolate if actually needed
        if (IndividualStarSaveTablePositions){
          printf("STAR: Table positions not being set correctly. Setting Rad table.\n");
        }

        this->AssignRadTablePosition();
      }
      break;

    case 3:
      if (this->yield_table_position[0] < 0 || this->yield_table_position[1] < 0){
        if (IndividualStarSaveTablePositions){
          printf("STAR: Table positions not being set correctly. Setting Yield table.\n");
        }
        this->AssignYieldTablePosition();
      }
      break;
    default:
      ENZO_FAIL("Star::AssertInterpolationPositions: Need to select a mode = 1, 2, 3\n");
  }

  return;
}
