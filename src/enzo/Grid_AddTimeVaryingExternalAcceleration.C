#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"
 
/* function prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
              float *TemperatureUnits, float *TimeUnits,
              float *VelocityUnits, double *MassUnits, FLOAT Time);

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int search_lower_bound(float *arr, float value, int low, int high, int total);

// In IndividualStarProperties.C
float LinearInterpolationCoefficient(const int &i, const float &x1, const float *x1a);

int InitializeTimeVaryingExternalAcceleration(float Time);


int grid::AddTimeVaryingExternalAcceleration(void){

  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
    VelocityUnits = 1, AccelUnits = 1, MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
  AccelUnits = LengthUnits/TimeUnits/TimeUnits;

  const float myr = 3.15576E13;

  float time_myr = this->Time * TimeUnits / myr;

  if( this->Time < ExternalGravityTimeOn ||
      this->Time > ExternalGravityTimeOff ){
    return SUCCESS;
  }

  // make sure everything is initialized 
  InitializeTimeVaryingExternalAcceleration(Time);

  FLOAT a = 1.0, dadt = 1.0;
  if (ComovingCoordinates){
    if(CosmologyComputeExpansionFactor(this->Time + 0.5*this->dtFixed, &a, &dadt) == FAIL){
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");
    }
  }

  /* find position of external potential source */
  int index = 0;
  float coeff;

  index = search_lower_bound(ExternalGravityTime, time_myr, 0, ExternalGravityNumberofTimePoints - 1,
                                                               EXTERNAL_GRAVITY_ENTRIES);
  coeff = LinearInterpolationCoefficient(index, time_myr, ExternalGravityTime);

  // find coordinates and normalize to box frame coordinates, which is what grid cell and particle positions will be in
  // ExternalGravityPosition is initally in coordinates relative to the center of the secondary galaxy
  ExternalGravityPosition[0] = ((1.0 - coeff)*ExternalGravityTimePositions[0][index] +
                               (      coeff)*ExternalGravityTimePositions[0][index+1])*Mpc/LengthUnits + DiskGravityPosition[0];
  ExternalGravityPosition[1] = ((1.0 - coeff)*ExternalGravityTimePositions[1][index] +
                               (      coeff)*ExternalGravityTimePositions[1][index+1])*Mpc/LengthUnits + DiskGravityPosition[1];
  ExternalGravityPosition[2] = ((1.0 - coeff)*ExternalGravityTimePositions[2][index] +
                               (      coeff)*ExternalGravityTimePositions[2][index+1])*Mpc/LengthUnits + DiskGravityPosition[2];


  FLOAT x, y, z, xpos, ypos, zpos, r, rsquared;
  double accel;

  for (int dim = 0; dim < GridRank; dim++) {
    int n = 0;

    for (int k = 0; k < GridDimension[2]; k++){
      if (GridRank > 2)
        zpos = this->CellLeftEdge[2][k]+0.5*this->CellWidth[2][k] - ExternalGravityPosition[2];
      if (dim == 2 && HydroMethod == Zeus_Hydro)
        zpos -= 0.5*this->CellWidth[2][k];

      for( int j = 0; j < GridDimension[1]; j++){
        if (GridRank > 1)
          ypos = this->CellLeftEdge[1][j] + 0.5*this->CellWidth[1][j] - ExternalGravityPosition[1];
        if (dim == 1 && HydroMethod == Zeus_Hydro)
          ypos -= 0.5*this->CellWidth[1][j];

        for (int i = 0; i < GridDimension[0]; i++, n++){
          xpos = this->CellLeftEdge[0][i] + 0.5*this->CellWidth[0][i] - ExternalGravityPosition[0];
          if (dim == 0 && HydroMethod == Zeus_Hydro)
            xpos -= 0.5*this->CellWidth[0][i];

          /* compute position */
          rsquared = (xpos*xpos + ypos*ypos + zpos*zpos) * LengthUnits * LengthUnits;
          r        = sqrt(rsquared);

          /* Worry only about a sperical potential for now - can add disk and things
             later following DiskGravity methods - */

          // ExternalGravityDensity -> dm central density in cgs
          // ExternalGravityRadius  -> dm halo scale radius in Mpc

          if (ExternalGravity == 2){ // NFW
            ENZO_FAIL("Time varying NFW external gravity potential not yet implemented\n");
          } else if (ExternalGravity == 3){ // spherical burkert profile

            accel = pi*(GravConst)*ExternalGravityDensity*POW(ExternalGravityRadius*Mpc,3)/(rsquared)
                      *(- 2.0 * atan(r/ExternalGravityRadius/Mpc)
                        + 2.0 *  log(1.0 + r/ExternalGravityRadius/Mpc)
                        +        log(1.0 + rsquared / POW(ExternalGravityRadius/Mpc,2))
                       );

            accel = ( r ==0.0?0.0:fabs(accel) / (r/LengthUnits) / AccelUnits);

          } // end potential type check

          if (dim == 0)
            AccelerationField[0][n] -= (accel * xpos);
          if (dim == 1)
            AccelerationField[1][n] -= (accel * ypos);
          if (dim == 2)
            AccelerationField[2][n] -= (accel * zpos);

        } // end i
      } // end j
    } // end k



  } // loop over dimension


  /* add in particle accelerations if any exist */

  if( NumberOfParticles > 0 && ParticleAcceleration[0] != NULL){
    for (int i = 0; i < NumberOfParticles; i++){
      x = ParticlePosition[0][i] + 0.5*this->dtFixed*ParticleVelocity[0][i]/a;
      y = ParticlePosition[1][i] + 0.5*this->dtFixed*ParticleVelocity[1][i]/a;
      z = ParticlePosition[2][i] + 0.5*this->dtFixed*ParticleVelocity[2][i]/a;

      // recenter particle coordinates and get r
      xpos = x - ExternalGravityPosition[0];
      ypos = y - ExternalGravityPosition[1];
      zpos = z - ExternalGravityPosition[2];

      rsquared = (xpos*xpos + ypos*ypos + zpos*zpos) * LengthUnits * LengthUnits;
      r        = sqrt(rsquared);


      if (ExternalGravity == 2){
        ENZO_FAIL("NFW time varying external gravity not yet working for particles\n");

      } else if (ExternalGravity == 3){
        accel = pi*(GravConst)*ExternalGravityDensity*POW(ExternalGravityRadius*Mpc,3)/(rsquared)
                  *(- 2.0 * atan(r/ExternalGravityRadius/Mpc)
                    + 2.0 *  log(1.0 + r/ExternalGravityRadius/Mpc)
                    +        log(1.0 + rsquared / POW(ExternalGravityRadius/Mpc,2))
                   );


        accel = ( r ==0.0?0.0:fabs(accel) / (r/LengthUnits) / AccelUnits);

      } // end potential type

      ParticleAcceleration[0][i] -= ( accel*xpos);
      ParticleAcceleration[1][i] -= ( accel*ypos);
      ParticleAcceleration[2][i] -= ( accel*zpos);

    } // end particle loop

  } // end particle check

  return SUCCESS;
}
