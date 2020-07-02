/***********************************************************************
/
/ Calculates and applies AGN feedback.
/
/  written by: Greg Meece
/  date:       January 2015
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

int FindField(int field, int farray[], int numfields);

/* Do AGN Feedback
 * This function manages the process of AGN feedback.
 * The only parameter it takes is the AGN particle itself.
 * The return value is success or failure.
 */
int grid::DoAGNFeedback( ActiveParticleType* ThisParticle) {

    if (MyProcessorNumber != ProcessorNumber)
       return SUCCESS;
     
    if (debug)
       printf("Entering DoAGNFeedback. \n");
    // Cast this particle to an AGN particle so that we can access AGN specific
    // properties.
    ActiveParticleType_AGNParticle* tp = static_cast <ActiveParticleType_AGNParticle*>(ThisParticle);

    // Get the units
    float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
         VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
    float MassUnits;

    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, Time);
    MassUnits = DensityUnits * pow(LengthUnits, 3.0);

    // Make sure the particle is on this grid.
    float cr, fr;
    cr = tp -> CoolingRadius; //* kpc / LengthUnits;
    fr = tp -> FeedbackRadius; // * kpc / LengthUnits;

   float max_radius = max(cr, fr);

   FLOAT xsink = tp -> ReturnPosition()[0];
   FLOAT ysink = tp -> ReturnPosition()[1];
   FLOAT zsink = tp -> ReturnPosition()[2];

   if ((GridLeftEdge[0]    > xsink + max_radius) ||
         (GridLeftEdge[1]  > ysink + max_radius) ||
         (GridLeftEdge[2]  > zsink + max_radius) ||
         (GridRightEdge[0] < xsink - max_radius) ||
         (GridRightEdge[1] < ysink - max_radius) ||
         (GridRightEdge[2] < zsink - max_radius))
      return SUCCESS;

   // If we have made it here, the particle is on this grid
   //
   //    // Figure out how to access grid quantities
   int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;

   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                        Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
      }

   int ColorNum = FindField(SNColour, FieldType, NumberOfBaryonFields);

   // Some physical constants

   float c_cgs = 2.99792458e+10;
   float g_cgs = 6.67259e-8;
   float mbh = 1.0e+8 * 1.99e+33;
   float mp = 1.67e-24;
   float sigma_t = 6.6525e-25;

   // By how much can mass, energy, and momentum conservation be violated?
   //    // This is the fractional change
   float change_tolerance = 1.0e-4;

   // Do Jet precession- this could be its own function
   // Theta = 0 -> Jet pointing in positive y direction
   // Theta = pi/2 -> Jet pointing in positive x direction
   //float t_theta = 0.017;
   float t_theta = 0.01; // Li+Bryan 2014 value

   tp -> JetTheta = 2.0 * M_PI * fmod(Time, t_theta) / t_theta;
   //tp -> JetPhi = (20.0/360.0) * 2.0 * M_PI;

   tp -> JetPhi = AGNPrecessionAngleRad;

   if (debug)
      printf("Jet theta, phi: %"FSYM", %"FSYM" (radians)\n", tp -> JetTheta, tp -> JetPhi);

   // Step 1: Figure out Mdot and remove mass from the grid
   // mdot is in units of code mass/time
   // Each of these functions also removes mdot*dt from the grid

   float mdot;
   /*if (AGNParticleFeedbackType == 0){
      mdot = AGNParticleGetPrecipitationMass(ThisParticle);
   }else */
   if (AGNParticleFeedbackType == 1){
      //mdot = AGNParticleGetColdMassRate(ThisParticle);
      mdot = 1.0*2.0e33*TimeUnits/(3.15e7*MassUnits);
/*
   else if (AGNParticleFeedbackType == 2)
      mdot = AGNParticleGetFixedPowerMass(ThisParticle);

   else if (AGNParticleFeedbackType == 3)
      mdot = AGNParticleGetBondiMassRate(ThisParticle);

   else if (AGNParticleFeedbackType == 4)
      mdot = AGNParticleGetAccretionDelayMass(ThisParticle);

   else if (AGNParticleFeedbackType == 5)
      mdot = AGNParticleGetStaggeredPowerMass(ThisParticle);

   else if (AGNParticleFeedbackType == 6)
      mdot = AGNParticleGetBondiMassRate(ThisParticle);

   // Exit for unknown heating methods
*/
   } else {
      mdot = 0.0;
      ENZO_FAIL("Error! Unimplemented Accretion rate function.");
      } 

   float heating_rate;
   float heating_rate_cgs;
   float ledd;

   float cell_volume = pow(CellWidth[0][0], 3.0); //Code
   float speed_of_light_code = c_cgs / (LengthUnits / TimeUnits);

   // Total heating rate in code units (energy/time)
  
   heating_rate = tp -> FeedbackEfficiency * mdot * pow(speed_of_light_code, 2.0);

   // Print out the heating rate in code, cgs, and as a fraction of the
   //    // Eddington luminosity
   heating_rate_cgs = tp -> FeedbackEfficiency * (mdot * (MassUnits / TimeUnits)) * pow(c_cgs, 2.0);
   
   ledd = (4.0 * M_PI * g_cgs * mbh * mp * c_cgs) / sigma_t;
   
   printf("Heating rate: %"GSYM" (%"GSYM" erg/s) (%"GSYM" L_edd) (dt=%"GSYM")\n",
   heating_rate, heating_rate_cgs, heating_rate_cgs / ledd, dtFixed);
   
   // Write edot to a file
   fprintf(AGNEdotFile, "%4.4"ISYM" %8.10"GSYM" %8.10"GSYM"\n", tp -> ReturnID(), Time, heating_rate);

   // Store the heating rate (in code). This is not used anywhere in the code-
   //    // it just allows it to be written out with each datadump.
   tp -> Edot = heating_rate;

   // Step 3: Apply heating
/*   if (AGNParticleFeedbackWeight == 0)
      AGNParticleVolumeWeightedJet(ThisParticle, mdot, heating_rate);

   else if (AGNParticleFeedbackWeight == 1)
      AGNParticleMassWeightedJet(ThisParticle, mdot, heating_rate);

   else if (AGNParticleFeedbackWeight == 2)
      AGNParticleDiskJet(ThisParticle, mdot, heating_rate);

   else if (AGNParticleFeedbackWeight == 3)
      AGNParticleSphereFeedback(ThisParticle, mdot, heating_rate);
   
   else */ 
   if (AGNParticleFeedbackWeight == 4)
      AGNParticleCylinderFeedback(ThisParticle, mdot, heating_rate);
   
   else {
      printf("Error! Unimplemented AGN feedback option\n");
      exit(1);
      } 

   // To avoid the timestep going to 0, check if any cells have a temperature above 10^9 K.
   // If they do, average their mass, momentum, and energy with the coldest adjacent cell.
   // Do the same if the density is below 10^-30 g/cc.
   bool use_temp_averaging = true;

   if (use_temp_averaging) {
      float* temperature = new float[GridDimension[0] * GridDimension[1] * GridDimension[2]];
      ComputeTemperatureField(temperature);

      int index; // index of this cell
      int cni; // index of coldest neighbor
      int ti; // test index

      float min_dt = 1.0;
      float cgs_dens;

      for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
         for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
            for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {

               index = k * GridDimension[1] * GridDimension[0] + j * GridDimension[0] + i;
               cgs_dens = BaryonField[DensNum][index] * DensityUnits;

               if (temperature[index] > 1.0e+9 || cgs_dens < 7.0e-30) {
                  //printf("Applying temperature averaging...\n");
                  float cell_volume;
                  cell_volume = CellWidth[0][0] * CellWidth[1][0] * CellWidth[2][0];
                  /*// left
                  cni = k * GridDimension[1] * GridDimension[0] + j * GridDimension[0] + i - 1;
                  //right
                  ti = k * GridDimension[1] * GridDimension[0] + j * GridDimension[0] + i + 1;
                  if (temperature[i] < temperature[cni])
                     cni = ti;

                  // front
                  ti = k * GridDimension[1] * GridDimension[0] + (j+1) * GridDimension[0] + i;
                  if (temperature[i] < temperature[cni])
                     cni = ti;

                  // back
                  ti = k * GridDimension[1] * GridDimension[0] + (j-1) * GridDimension[0] + i;
                  if (temperature[i] < temperature[cni])
                     cni = ti;

                  // top
                  ti = (k+1) * GridDimension[1] * GridDimension[0] + j * GridDimension[0] + i;
                  if (temperature[i] < temperature[cni])
                     cni = ti;

                  // bottom
                  ti = (k-1) * GridDimension[1] * GridDimension[0] + j * GridDimension[0] + i;
                  if (temperature[i] < temperature[cni])
                     cni = ti;
                  */
                  // added by Deovrat Prasad for picking out the coldest adjacent cell 
                  cni = index;
                  //left
                  ti = k * GridDimension[1] * GridDimension[0] + j * GridDimension[0] + i - 1;
                  if (temperature[ti] < temperature[cni])
                     cni = ti;
                  //right
                  ti = k * GridDimension[1] * GridDimension[0] + j * GridDimension[0] + i + 1;
                  if (temperature[ti] < temperature[cni])
                     cni = ti;
		  //front
		  ti = k * GridDimension[1] * GridDimension[0] + (j+1) * GridDimension[0] + i;
                  if (temperature[ti] < temperature[cni])
                     cni = ti;
 		  //back
                  ti = k * GridDimension[1] * GridDimension[0] + (j-1) * GridDimension[0] + i;
                  if (temperature[ti] < temperature[cni])
                     cni = ti; 		                    
                  //top
                  ti = (k+1) * GridDimension[1] * GridDimension[0] + j * GridDimension[0] + i;
                  if (temperature[ti] < temperature[cni])
                     cni = ti;
		  //bottom
                  ti = (k-1) * GridDimension[1] * GridDimension[0] + j * GridDimension[0] + i;
                  if (temperature[ti] < temperature[cni])
                     cni = ti;		  

                  // Average mass, momentum, and energy

                  float total_mass, total_px, total_py, total_pz, total_ge, total_ke;

                  total_mass = BaryonField[DensNum][index] + BaryonField[DensNum][cni];
                  total_mass *= cell_volume;

                  total_px = BaryonField[DensNum][index] * BaryonField[Vel1Num][index]
                           + BaryonField[DensNum][cni] * BaryonField[Vel1Num][cni];
                  total_px *= cell_volume;

                  total_py = BaryonField[DensNum][index] * BaryonField[Vel2Num][index]
                           + BaryonField[DensNum][cni] * BaryonField[Vel2Num][cni];
                  total_py *= cell_volume;

                  total_pz = BaryonField[DensNum][index] * BaryonField[Vel3Num][index]
                           + BaryonField[DensNum][cni] * BaryonField[Vel3Num][cni];
                  total_pz *= cell_volume;

                  total_ke = 0.5 * BaryonField[DensNum][index] * (pow(BaryonField[Vel1Num][index], 2.0)
                           + pow(BaryonField[Vel1Num][index], 2.0) + pow(BaryonField[Vel1Num][index], 2.0));
                  total_ke += 0.5 * BaryonField[DensNum][cni] * (pow(BaryonField[Vel1Num][cni], 2.0)
                           + pow(BaryonField[Vel1Num][cni], 2.0)+ pow(BaryonField[Vel1Num][cni], 2.0));
                  total_ke *= cell_volume;

                  if (HydroMethod==Zeus_Hydro) {
                     total_ge = BaryonField[TENum][index] * BaryonField[DensNum][index] + BaryonField[TENum][cni] * BaryonField[DensNum][cni];
                     total_ge *= cell_volume;
                     }

                  else if (HydroMethod == PPM_DirectEuler &&DualEnergyFormalism) {
                     total_ge = BaryonField[GENum][index] * BaryonField[DensNum][index] + BaryonField[GENum][cni] * BaryonField[DensNum][cni];
                     total_ge *= cell_volume;
                     }

                  else if (HydroMethod == PPM_DirectEuler && !DualEnergyFormalism) {
                     total_ge = BaryonField[TENum][index] * BaryonField[DensNum][index] + BaryonField[TENum][cni] * BaryonField[DensNum][cni];
                     total_ge *= cell_volume;
                     total_ge -= total_ke;
                     }

                  // Set cell quantities
                  BaryonField[DensNum][index] = total_mass * (0.5 / cell_volume);
                  BaryonField[DensNum][cni] = total_mass * (0.5 / cell_volume);

                  BaryonField[Vel1Num][index] = total_px * (0.5 / total_mass);
                  BaryonField[Vel2Num][index] = total_py * (0.5 / total_mass);
                  BaryonField[Vel3Num][index] = total_pz * (0.5 / total_mass);

                  BaryonField[Vel1Num][cni] = total_px * (0.5 / total_mass);
                  BaryonField[Vel2Num][cni] = total_py * (0.5 / total_mass);
                  BaryonField[Vel3Num][cni] = total_pz * (0.5 / total_mass);

                  if (HydroMethod==Zeus_Hydro) {
                     BaryonField[TENum][index] = total_ge * (0.5 / total_mass);
                     BaryonField[TENum][cni] = total_ge * (0.5 / total_mass);
                     }

                  else if (HydroMethod == PPM_DirectEuler &&DualEnergyFormalism) {
                     BaryonField[GENum][index] = total_ge * (0.5 / total_mass);
                     BaryonField[GENum][cni] = total_ge * (0.5 / total_mass);

                     BaryonField[TENum][index] = (total_ge + total_ke) * (0.5 / total_mass);
                     BaryonField[TENum][cni] = (total_ge + total_ke) * (0.5 / total_mass);
                     }

                  else if (HydroMethod == PPM_DirectEuler && !DualEnergyFormalism) {
                     BaryonField[TENum][index] = (total_ge + total_ke) * (0.5 / total_mass);
                     BaryonField[TENum][cni] = (total_ge + total_ke) * (0.5 / total_mass);
                     }
                  } // End if T > 10^9 K */
               } // End loop over x
            } // End loop over y
         } // End loop over z
                 
      delete [] temperature;
      } 
      if (debug)
         printf ("Leaving DoAGNFeedback [%"ISYM"]\n", MyProcessorNumber );

      return SUCCESS;
}

#undef DEBUG_AP
