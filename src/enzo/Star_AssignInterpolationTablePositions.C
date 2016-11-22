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



void Star::AssignInterpolationTablePositions(void){
  /* Sets table positions for each star */

  if( abs(this->type) == PARTICLE_TYPE_INDIVIDUAL_STAR ){

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

void Star::AssignSETablePosition(void){

  /* find stellar properties - used for winds, SN, and radiation */
  IndividualStarGetSETablePosition(this->se_table_position[0],
                                   this->se_table_position[1],
                                   this->BirthMass, this->Metallicity);
  return;
}

void Star::AssignRadTablePosition(void){

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

  StellarYieldsGetYieldTablePosition(this->yield_table_position[0],
                                     this->yield_table_position[1],
                                     this->BirthMass, this->Metallicity);

  return;
}

void Star::AssertInterpolationPositions(int mode){
 // mode = 1: se    table
 //        2: rad   table
 //        3: yield table

  switch(mode){
    case 1:
      if (this->se_table_position[0] < 0 || this->se_table_position[1] < 0){
        this->AssignSETablePosition();
      }
      break;

    case 2:
      for(int i = 0; i < 3; i++){
        if(this->rad_table_position[i] < 0 && this->rad_table_position[i] > -2){
          this->AssignRadTablePosition();
          continue;
        }
      }
      break;

    case 3:
      if (this->yield_table_position[0] < 0 || this->yield_table_position[1] < 0){
        this->AssignYieldTablePosition();
      }
      break;
    default:
      ENZO_FAIL("Star::AssertInterpolationPositions: Need to select a mode = 1, 2, 3\n");
  }

  return;
}
