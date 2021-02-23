/**********************************************************************
/
/
/  written by: Andrew Emerick
/  date:       May, 2018
/  modified1:
/
/  PURPOSE:
/          Enforce a maximum temperature
/
/  RETURNS:
/     SUCCESS or FAIL
***********************************************************************/

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

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);


int grid::ApplyTemperatureLimit(void){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (TemperatureLimit < 0)
    return SUCCESS;

  // else, go through and instill density / temperature thresh
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass unit

  /* obtain baryon field indexes */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }



  int size = 1;
  for (int i = 0; i < GridRank; i++)
    size *= GridDimension[i];

  float * temperature;
  temperature = new float[size];

  if(this->ComputeTemperatureField(temperature) == FAIL){
    ENZO_FAIL("Error in compute temeprature called from apply tempeature limit");
  }

  // don't mess with halo gas on the low end
  float low_density_threshold  = 1.0E-30 / DensityUnits;

  // don't mess with actual SN injection sites
  float high_density_threshold = 1.0E-26 / DensityUnits;

  for(int i = 0; i < size; i ++){
    if( (temperature[i] > TemperatureLimit) &&
        ((BaryonField[DensNum][i] > low_density_threshold) &&
         (BaryonField[DensNum][i] < high_density_threshold))  ){
      float factor = min(temperature[i] / TemperatureLimit,
                         TemperatureLimitFactor);
      float old_dens = BaryonField[DensNum][i];

      BaryonField[DensNum][i] *= factor;

      if (DualEnergyFormalism){
         float old_kinetic = (BaryonField[TENum][i] - BaryonField[GENum][i])*old_dens;
         float thermal = old_dens*BaryonField[GENum][i];

         BaryonField[GENum][i] = thermal / BaryonField[DensNum][i];


         float kinetic = 0.0;
         for (int dim = 0; dim < GridRank; dim++){
           BaryonField[Vel1Num+dim][i] *= (old_dens / BaryonField[DensNum][i]);
           kinetic += 0.5 * BaryonField[Vel1Num+dim][i] * BaryonField[Vel1Num+dim][i];
         }
         BaryonField[TENum][i] = (thermal+kinetic)/BaryonField[DensNum][i];
      } else{
         BaryonField[TENum][i] /= factor;
      }
    }
  }

  delete [] temperature;

  return SUCCESS;
}
