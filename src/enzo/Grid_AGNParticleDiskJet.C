/***********************************************************************
/
/ Calculates and applies AGN feedback for a jet.
/ Feedback emenates from a disk of cells
/
/ Parameters:
/    ThisParticle: The AGN particle
/    mdot: The mass accretion rate in mass/time, code units
/    heating_rate: Total jet power in energy/time, code units
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

// Inputs:
//    mdot is the mass accretion rate in units of code mass/code time
//    heating_rate is the total jet heating rate in units of code
//                 energy/code_time
int grid::AGNParticleDiskJet(ActiveParticleType* ThisParticle, float mdot, float heating_rate) {

   /* Return if this doesn't involve us */
   if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;
   if (MyProcessorNumber == ProcessorNumber)
      printf("proceeding with the disk jet.\n");
   // Cast this particle to an AGN particle so that we can access AGN specific
   // properties.
   ActiveParticleType_AGNParticle* tp = static_cast <ActiveParticleType_AGNParticle*>(ThisParticle);

   float xsink = tp -> ReturnPosition()[0];
   float ysink = tp -> ReturnPosition()[1];
   float zsink = tp -> ReturnPosition()[2];

   // Figure out how to access grid quantities
   int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;

   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                        Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
      }

   // Get the units
   float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
         VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
   float MassUnits, EnergyUnits;

   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, Time);
   MassUnits = DensityUnits * pow(LengthUnits, 3.0);
   EnergyUnits = MassUnits * pow(LengthUnits / TimeUnits, 2.0);

   float code_disk_radius, code_disk_distance;
   code_disk_radius = AGNParticleDiskRadius;
   code_disk_distance = AGNParticleDiskDistance;

   // Declare a bunch of variables
   float heating_this_cell, a, x, y, z, phi, radius_in_plane,
         edot_this_cell, ke_this_cell, ge_this_cell, velocity_this_cell, radius;
   float mdot_this_cell, cell_mass;

   float cell_volume = pow(CellWidth[0][0], 3.0); //Code

   // Total number of cells (including ghosts) in this grid.
   int size = 1;

   for (int i = 0; i < GridRank; i++)
      size *= GridDimension[i];

   // Sum up the amount of energy being injected- should equal heating_rate!
   float injected_heating = 0.0;
   int index;

   // Calculate the jet velocity
   float jet_velocity;
   jet_velocity = sqrt(2.0 * tp -> FeedbackEfficiency * tp -> KineticFraction) * 2.99792458e+10; //cm/s
   jet_velocity /= VelocityUnits; // code

   // Loop over both disks and flag cells in the injection zone
   //  0 = Not in injection zone
   //  1 = In injection zone, facing up
   // -1 = In injection zone, facing down
   int* injection_cells = new int[size];

   for (int i = 0; i < size; i++)
      injection_cells[i] = 0;

   for (int facing = 0; facing <= 1; facing++) {
      // Figure out the position of the center of the disk
      float dpos[3];

      if (facing == 0) {
         dpos[0] = tp -> ReturnPosition()[0] + code_disk_distance * cos(tp -> JetTheta) * sin(tp -> JetPhi);
         dpos[1] = tp -> ReturnPosition()[1] + code_disk_distance * sin(tp -> JetTheta) * sin(tp -> JetPhi);
         dpos[2] = tp -> ReturnPosition()[2] + code_disk_distance * cos(tp -> JetPhi);
         }
      else {
         dpos[0] = tp -> ReturnPosition()[0] - code_disk_distance * cos(tp -> JetTheta) * sin(tp -> JetPhi);
         dpos[1] = tp -> ReturnPosition()[1] - code_disk_distance * sin(tp -> JetTheta) * sin(tp -> JetPhi);
         dpos[2] = tp -> ReturnPosition()[2] - code_disk_distance * cos(tp -> JetPhi);
         }
   
      // Find properties of the plane in which the disk lies
      // da*x + db*y + dz * z - d0 = 0
      float da, db, dc, d0;
   
      da = cos(tp -> JetTheta) * sin(tp -> JetPhi);
      db = sin(tp -> JetTheta) * sin(tp -> JetPhi);
      dc = cos(tp -> JetPhi);

      d0 = da * dpos[0] + db * dpos[1] + dc * dpos[2];
      //if (debug) 
      //   printf("d0=%"GSYM"\n", d0);
      // Sum up the volume taken up by cells intersecting the disk
      float disk_volume = 0.0;
   
      for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
         for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
            for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {
               index = k * GridDimension[1] * GridDimension[0]
                     + j * GridDimension[0] + i;
   
               // Check whether the disk plane intersects this cell
               int above = 0;
               int below = 0;
               float l;
               float min_dist = huge_number;
   
               for (int c0 = 0; c0 <= 1; c0++) {
                  for (int c1 = 0; c1 <= 1; c1++) {
                     for (int c2 = 0; c2 <= 1; c2++) {
                        // Position of this corner
                        x = CellLeftEdge[0][i] + (float)c0 * CellWidth[0][i];
                        y = CellLeftEdge[1][j] + (float)c1 * CellWidth[1][j];
                        z = CellLeftEdge[2][k] + (float)c2 * CellWidth[2][k];
   
                        // Distance of corner from center of disk
                        //radius = sqrt( pow(x - dpos[0], 2.0)
                        //             + pow(y - dpos[1], 2.0)
                        //             + pow(z - dpos[2], 2.0));
                        radius = sqrt(pow (x - tp->ReturnPosition()[0], 2.0)
                                     +pow (y - tp->ReturnPosition()[1], 2.0)
                                     +pow (z - tp->ReturnPosition()[2], 2.0));  
   
                        min_dist = min(radius, min_dist);
   
                        l = da * x + db * y + dc * z - d0;
                         
                        if (l > 0.0)
                           above++;
                        else
                           below++;
                        } //c2
                     } //c1
                  } //c0
   
               // If all points are on the same side
               if (above != 0 && below != 0) {
                  //printf("both above and below non-zero.\n");
                  // If at least one corner is within the radius of the disk
                  // Note that this will miss cases where the disk just grazes one
                  // of the sides of the cell. These cases should be few, however,
                  // and the behavior of the jet should not be affected.
                  if (min_dist < code_disk_radius) {
                     if (facing == 0){
                        injection_cells[index] = 1;
                     }else{
                        injection_cells[index] = -1;
                     }
		   }
                 }
               } // End loop over i
            } // End loop over j
         } // End loop over k
      } // End loop over both disks

   // Calculate total mass, energy, and volume of the injection zone
   int total_cells;
   float total_mass, total_energy, total_volume;

   total_cells = 0;
   total_mass = 0.0;
   total_energy = 0.0;
   total_volume = 0.0;

   for (int i = 0; i < size; i++) {
      if (injection_cells[i] != 0) {
         total_cells++;

         cell_mass = BaryonField[DensNum][i] * cell_volume;
         total_mass += cell_mass;

         total_volume += cell_volume;

         // ZEUS: TENum is internal energy
         // Total energy calculated by adding kinetic
         if (HydroMethod==Zeus_Hydro) {
            total_energy += BaryonField[TENum][i] * cell_mass;
            total_energy += 0.5 * cell_mass * (
                            pow(BaryonField[Vel1Num][i], 2.0)
                            + pow(BaryonField[Vel2Num][i], 2.0)
                            + pow(BaryonField[Vel3Num][i], 2.0));
            }

         // PPM and DEF: TENum is total energy, GENum is internal energy
         else if (HydroMethod == PPM_DirectEuler && DualEnergyFormalism) {
            total_energy += BaryonField[TENum][i]*cell_mass;
            }

         // PPM without DEF: TENum is total energy,
         //  internal energy calcualted by subtracting kinetic energy.
         else if (HydroMethod == PPM_DirectEuler && !DualEnergyFormalism) {
            total_energy += BaryonField[TENum][i]*cell_mass;
            }
         }
         //if (debug)
         //   printf("total volume = %"GSYM"\n", total_volume);
      } // End summing mass, volume, and energy

   // Check if the energy to be injected is high enough to boost the gas above
   // the threshold temperature.
   float to_inject;
   float threshold_energy;

   to_inject = heating_rate * dtFixed + tp -> StoredEnergy; //code energy

   threshold_energy = (3.0/2.0) * kboltz * AGNParticleInjectionTemperature
                      * (total_mass * MassUnits) / (Mu * mh); //cgs energy

   threshold_energy /= EnergyUnits; // code energy

   // If stored energy plus new energy is less than the threshold, store the
   // energy and mass, don't apply heating yet
   if (to_inject < threshold_energy) {
      tp -> StoredEnergy += heating_rate * dtFixed;
      tp -> StoredMass += mdot * dtFixed;
      }

   else {
      // Add in the energy from the jet and the stored energy
      total_mass += mdot * dtFixed + tp -> StoredMass;
      total_energy += heating_rate * dtFixed + tp -> StoredEnergy;

      tp -> StoredEnergy = 0.0;
      tp -> StoredMass = 0.0;

      // Calculate properties of gas in the injection zone
      float cell_dens;
      float cell_ge;
      float cell_ke;
      float cell_vel;

      cell_dens = total_mass/total_volume;
      cell_ge = (1.0 - tp -> KineticFraction) * (total_energy/total_mass);
      cell_ke = tp -> KineticFraction * (total_energy / total_mass);
      cell_vel = sqrt(2.0 * cell_ke);

      // Make sure that the gas is at least 10^4 K
      float min_temp = 1.0e4;
      float min_ge = min_temp;
      min_ge /= (TemperatureUnits * (Gamma - 1.0) * Mu);

      if (cell_ge < min_ge)
         cell_ge = min_ge;

      printf("total_mass, cell_den, cell_ge: %"GSYM" %"GSYM" %"GSYM"\n", total_mass, cell_dens, cell_ge);

      // Apply feedback to all cells in the injection zone
      for (int i = 0; i < size; i++) {
         if (injection_cells[i] != 0) {
            // Set the density
            BaryonField[DensNum][i] = cell_dens;

            // Velocity if facing the jet
            if (injection_cells[i] == 1) {
               BaryonField[Vel1Num][i] = cell_vel * sin(tp -> JetTheta) * sin(tp -> JetPhi);
               BaryonField[Vel2Num][i] = cell_vel * cos(tp -> JetTheta) * sin(tp -> JetPhi);
               BaryonField[Vel3Num][i] = cell_vel * cos(tp -> JetPhi);
               }

            // Velocity on other side
            else if (injection_cells[i] == -1) {
               BaryonField[Vel1Num][i] = -cell_vel * sin(tp -> JetTheta) * sin(tp -> JetPhi);
               BaryonField[Vel2Num][i] = -cell_vel * cos(tp -> JetTheta) * sin(tp -> JetPhi);
               BaryonField[Vel3Num][i] = -cell_vel * cos(tp -> JetPhi);
               }

            // ZEUS: TENum is internal energy
            if (HydroMethod == Zeus_Hydro) {
               BaryonField[TENum][i] = cell_ge;
               }

            // PPM and DEF: TENum is total energy, GENum is internal energy
            else if (HydroMethod == PPM_DirectEuler && DualEnergyFormalism) {
               BaryonField[TENum][i] = cell_ge + cell_ke;
               BaryonField[GENum][i] = cell_ge;
               }

            // PPM without DEF: TENum is total energy,
            //  internal energy calcualted by subtracting kinetic energy.
            else if (HydroMethod == PPM_DirectEuler && !DualEnergyFormalism) {
               BaryonField[TENum][i] = cell_ge + cell_ke;
               }
            }
         } // Done injecting energy
      }

   // Clean up memory
   delete[] injection_cells;

   if (debug)
      printf("Leaving AGNParticleDiskJet[%"ISYM"]\n", MyProcessorNumber);

   return SUCCESS;
   }

#undef DEBUG_AP
