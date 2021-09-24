/***********************************************************************
/
/ Calculates the accretion rate of cold gas within R_Cooling.
/ Accretion rate is M_cold/time_delay
/ Returns the accretion rate in code units
/
/  written by: Greg Meece
/  date:       April 2015
/
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "ActiveParticle_AGNParticle.h"

#define NO_DEBUG_AP

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

float grid::AGNParticleGetColdMassRate(ActiveParticleType* ThisParticle) {

   /* Return if this doesn't involve us */
   if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

   // Cast this particle to an AGN particle so that we can access AGN specific
   // properties.
   ActiveParticleType_AGNParticle* tp = static_cast <ActiveParticleType_AGNParticle*>(ThisParticle);

   float max_radius = max(tp -> CoolingRadius, tp -> FeedbackRadius);

   float xsink = tp -> pos[0];
   float ysink = tp -> pos[1];
   float zsink = tp -> pos[2];

   if ((GridLeftEdge[0]    > xsink + max_radius) ||
         (GridLeftEdge[1]  > ysink + max_radius) ||
         (GridLeftEdge[2]  > zsink + max_radius) ||
         (GridRightEdge[0] < xsink - max_radius) ||
         (GridRightEdge[1] < ysink - max_radius) ||
         (GridRightEdge[2] < zsink - max_radius))
      return SUCCESS;

   // Figure out how to access grid quantities
   int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;

   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                        Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
      }

   // Get the units
   float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
         VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
   float MassUnits;

   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, Time);
   MassUnits = DensityUnits * pow(LengthUnits, 3.0);

   float cell_volume = pow(CellWidth[0][0], 3.0); //Code
   float radius; //Code

   int index;

   float cm_per_kpc = 3.08567758e21;
   float seconds_per_year = 365.25 * 24.0 * 3600.0;
   float cell_mass;
   float mdot;

   // Calculate the delay time in code units
   float cold_delay;
   cold_delay = AGNParticleColdMassDelay * seconds_per_year / TimeUnits;

   float* temperature = new float[GridDimension[0] * GridDimension[1] * GridDimension[2]];
   ComputeTemperatureField(temperature);

   // Calculate and sum Mdot over all cells
   mdot = 0.0;
   float acc_frac = 1.0; //acc_frac accounts for gas fraction within the accretion zone 
                          // to accrete on SMBH. Added by Deovrat Prasad.

   for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
         for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {
            index = k * GridDimension[1] * GridDimension[0]
                  + j * GridDimension[0] + i;

            radius = pow(CellLeftEdge[0][i] + CellWidth[0][i] * 0.5 - xsink, 2.0)
                   + pow(CellLeftEdge[1][j] + CellWidth[1][j] * 0.5 - ysink, 2.0)
                   + pow(CellLeftEdge[2][k] + CellWidth[2][k] * 0.5 - zsink, 2.0);
            radius = sqrt(radius);
            

            if (radius < tp -> CoolingRadius) {
               cell_mass = BaryonField[DensNum][index] * cell_volume; //code

               if (temperature[index] < AGNParticleColdMassTemp) {
                  mdot += acc_frac*BaryonField[DensNum][index] * cell_volume
                          * (dtFixed / cold_delay);

                  // Condensed mass drops out
                  BaryonField[DensNum][index] -= acc_frac*BaryonField[DensNum][index]
                                                 * (dtFixed / cold_delay);
               }
            } // End mdot calculation
          } // End loop over i
       } // End loop over j
    } // End loop over k

   // Make mdot be condensation rate per time unit
   mdot /= dtFixed;

   delete[] temperature;
   return mdot;
}
