/***********************************************************************
/
/ Applies AGN feedback in a cylinder.
/ For use in low resolution/cosmological simulations
/
/ Parameters:
/    ThisParticle: The AGN particle
/    mdot: The mass accretion rate in mass/time, code units
/    heating_rate: Total jet power in energy/time, code units
/
/  written by: Greg Meece
/  date:       February 2016
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
#include "AGN_Volume_Cylinder.h"
#include "AGN_Volume_Grid.h"

#define NO_DEBUG_AP

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

// Inputs:
//    mdot is the mass accretion rate in units of code mass/code time
//    heating_rate is the total jet heating rate in units of code
//                 energy/code_time
int grid::AGNParticleCylinderFeedback(ActiveParticleType* ThisParticle, float mdot, float heating_rate) {

   /* Return if this doesn't involve us */
   if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

   // Cast this particle to an AGN particle so that we can access AGN specific
   // properties.
   ActiveParticleType_AGNParticle* tp = static_cast <ActiveParticleType_AGNParticle*>(ThisParticle);

   float xsink = tp -> pos[0];
   float ysink = tp -> pos[1];
   float zsink = tp -> pos[2];

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
   printf("Time in code units (Get units test): %"GSYM"\n", Time);

   // Get the volume of a cell
   float cell_volume;
   cell_volume = CellWidth[0][0] * CellWidth[1][0] * CellWidth[2][0];

   // Do the feedback here...
   //
   printf("Doing AGN Cylindrical feedback!\n");

   // Check if the jet is on
   if (heating_rate <= tiny_number)
      return SUCCESS;

   float jet_heating_rate, shock_heating_rate;
   jet_heating_rate = heating_rate * (1.0 - AGNParticleShockFraction);
   shock_heating_rate = heating_rate * AGNParticleShockFraction;

   // Need to do feedback for two cylinders- one in each direction.
   for (int facing = 0; facing < 2; facing++) {
      // Create the AGN Cylinder
      float rad, height;
      float base[3]; // Position of center of bottom
      float n[3]; // Unit vector pointing along the cylinder

      rad = AGNCylinderRadiusKpc * kpc_cm / LengthUnits; // code units
      height = AGNCylinderHeightKpc * kpc_cm / LengthUnits; // code units

      /* 
      The unit vector of the jet is currently set with respect to z-axis of 
      the simulation. Also the jet precession is set with respect to z-axis.
      This needs to be modified in case the jet is to be injected along some 
      arbitrary axis.. 
      */
      //n[0] = cos(tp -> JetTheta) * sin(tp -> JetPhi);
      //n[1] = sin(tp -> JetTheta) * sin(tp -> JetPhi);
      //n[2] = cos(tp -> JetPhi);

/*    Calculating the unit vector of the jet alignment  based on the net 
 *    angular momenutum of the cold gas within the cooling radius.
 *    Added by Deovrat Prasad.
 */ 
      float xx,yy,zz, radius;
      float l_tot, Lx, Ly, Lz, c_mass, c_volume;
      int tot_ind;
      for (int k = 0; k < GridDimension[2]; k++) {
         for (int j = 0; j < GridDimension[1]; j++) {
            for (int i = 0; i < GridDimension[0]; i++) {
               tot_ind = k * GridDimension[1] * GridDimension[0] + j *
                   GridDimension[1] + i;

               xx = ( CellLeftEdge[0][i] + CellWidth[0][i] * 0.5 ) - xsink;
               yy = ( CellLeftEdge[1][j] + CellWidth[1][j] * 0.5 ) - ysink;
               zz = ( CellLeftEdge[2][k] + CellWidth[2][k] * 0.5 ) - zsink;

               radius = pow(xx, 2.0) + pow(yy, 2.0) + pow(zz, 2.0);
               radius = sqrt(radius);

               if( radius < tp->CoolingRadius ){
                  c_volume = CellWidth[0][i]*CellWidth[1][j]*CellWidth[2][k];
                  c_mass = BaryonField[DensNum][tot_ind]*c_volume;
                  Lx += c_mass*(BaryonField[Vel3Num][tot_ind]*yy - BaryonField[Vel2Num][tot_ind]*zz);
                  Ly += c_mass*(BaryonField[Vel1Num][tot_ind]*zz - BaryonField[Vel3Num][tot_ind]*xx);
                  Lz += c_mass*(BaryonField[Vel2Num][tot_ind]*xx - BaryonField[Vel1Num][tot_ind]*yy);
               }  //if loop

               //l_tot = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
               //n[0] = Lx/l_tot;
               //n[1] = Ly/l_tot;
               //n[2] = Lz/l_tot;
            } //i loop
         } //j loop
      } //k loop

      l_tot = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);

      n[0] = Lx/l_tot;
      n[1] = Ly/l_tot;
      n[2] = Lz/l_tot;
      if (facing == 0) {
         n[0] *=-1.0;
         n[1] *=-1.0;
         n[2] *=-1.0;
      }
      base[0] = n[0] * AGNCylinderDistanceKpc * kpc_cm / LengthUnits;
      base[1] = n[1] * AGNCylinderDistanceKpc * kpc_cm / LengthUnits;
      base[2] = n[2] * AGNCylinderDistanceKpc * kpc_cm / LengthUnits;

      base[0] += xsink;
      base[1] += ysink;
      base[2] += zsink;

      AGN_Volume_Cylinder* cyl = new AGN_Volume_Cylinder(base, n, rad, height);

      // Create a grid for computing the volume
      AGN_Volume_Grid* vol_grid = new AGN_Volume_Grid(GridDimension, GridLeftEdge, GridRightEdge);

      // Compute the volume of each cell enclosed by the cylinder.
      float*** vol = vol_grid -> get_intersection_volume(cyl, 0, 2, 3);

      // Find the total volume on the grid enclosed by the cylinder.
      // Check that it is reasonably close to the actual volume.
      float a_vol, g_vol;
      a_vol = M_PI * rad * rad * height; // Analytic volume
      g_vol = 0.0;

      for (int k = 0; k < GridDimension[2]; k++)
         for (int j = 0; j < GridDimension[1]; j++)
            for (int i = 0; i < GridDimension[0]; i++)
               g_vol += vol[k][j][i];

      //printf("Volume (analytic, total): %"GSYM" %"GSYM"\n", a_vol, g_vol);
    
      if (fabs(a_vol-g_vol) / a_vol > 0.1)
         printf("Warning! AGN Cylinder volume is off by more than 10%. This is probably a bad thing...\n");

      // Add energy
      float old_cell_mass;
      float old_cell_px, old_cell_py, old_cell_pz;
      float old_cell_ge, old_cell_ke;

      float cell_vel;
      float delta_cell_mass;
      float delta_cell_px, delta_cell_py, delta_cell_pz;
      float delta_cell_ge, delta_cell_ke;

      float new_cell_mass;
      float new_cell_px, new_cell_py, new_cell_pz;
      float new_cell_ge, new_cell_ke;
  
      int index;

      for (int k = 0; k < GridDimension[2]; k++) {
         for (int j = 0; j < GridDimension[1]; j++) {
            for (int i = 0; i < GridDimension[0]; i++) {
               index = k * GridDimension[1] * GridDimension[0] + j *
                  GridDimension[1] + i;

               // Skip ahead if this cell isn't in the cylinder.
               if (vol[k][j][i] == 0.0)
                  continue;

               // Calculate old quantities
               old_cell_mass = BaryonField[DensNum][index] * cell_volume;
               //printf("index, rho, v: %"ISYM", %"GSYM", %"GSYM"\n", index, BaryonField[DensNum][index], cell_volume);

               old_cell_px = old_cell_mass * BaryonField[Vel1Num][index];
               old_cell_py = old_cell_mass * BaryonField[Vel2Num][index];
               old_cell_pz = old_cell_mass * BaryonField[Vel3Num][index];

               // ZEUS: TENum is internal energy
               if (HydroMethod == Zeus_Hydro) {
                  old_cell_ge = BaryonField[TENum][index] * old_cell_mass;
                  old_cell_ke = 0.5 * old_cell_mass * (
                           pow(BaryonField[Vel1Num][index], 2.0) + 
                           pow(BaryonField[Vel2Num][index], 2.0) + 
                           pow(BaryonField[Vel3Num][index], 2.0));
                  }

               // PPM and DEF: TENum is total energy, GENum is internal energy
               else if (HydroMethod == PPM_DirectEuler && DualEnergyFormalism) {
                  old_cell_ke = (BaryonField[TENum][index] - BaryonField[GENum][index]) * old_cell_mass;
                  old_cell_ge = BaryonField[GENum][index] * old_cell_mass;
                  }

               // PPM without DEF: TENum is total energy,
               //  internal energy calcualted by subtracting kinetic energy.
               else if (HydroMethod == PPM_DirectEuler && !DualEnergyFormalism) {
                  old_cell_ke = 0.5 * old_cell_mass * (
                           pow(BaryonField[Vel1Num][index], 2.0) + 
                           pow(BaryonField[Vel2Num][index], 2.0) + 
                           pow(BaryonField[Vel3Num][index], 2.0));
                  old_cell_ge = (BaryonField[TENum][index] * old_cell_mass) - old_cell_ke;
                  }

               // Calculate changes to quantities
               delta_cell_mass = mdot * dtFixed * vol[k][j][i] / g_vol;

               delta_cell_ge = jet_heating_rate * dtFixed;
               delta_cell_ge *= vol[k][j][i] / g_vol;
               delta_cell_ke = delta_cell_ge;

               delta_cell_ge *= (1.0 - tp -> KineticFraction);
               delta_cell_ke *= tp -> KineticFraction;

               cell_vel = sqrt(2.0 * delta_cell_ke / delta_cell_mass);
               delta_cell_px = delta_cell_mass * cell_vel * n[0];
               delta_cell_py = delta_cell_mass * cell_vel * n[1];
               delta_cell_pz = delta_cell_mass * cell_vel * n[2];

               // Calculate new quantities
               new_cell_mass = old_cell_mass + delta_cell_mass;

               new_cell_px = old_cell_px + delta_cell_px;
               new_cell_py = old_cell_py + delta_cell_py;
               new_cell_pz = old_cell_pz + delta_cell_pz;

               new_cell_ge = old_cell_ge + delta_cell_ge;
               new_cell_ke = old_cell_ke + delta_cell_ke;

               //printf("old: m: %"GSYM"\n", old_cell_mass);
               //printf("old: px, py, pz: %"GSYM", %"GSYM", %"GSYM"\n", old_cell_px, old_cell_py, old_cell_pz);
               //printf("old: ge, ke: %"GSYM", %"GSYM"\n\n", old_cell_ge, old_cell_ke);

               // Update grid quantities
               BaryonField[DensNum][index] = new_cell_mass / cell_volume;

               BaryonField[Vel1Num][index] = new_cell_px / new_cell_mass;
               BaryonField[Vel2Num][index] = new_cell_py / new_cell_mass;
               BaryonField[Vel3Num][index] = new_cell_pz / new_cell_mass;

               // ZEUS: TENum is internal energy
               if (HydroMethod == Zeus_Hydro) {
                  BaryonField[TENum][index] = new_cell_ge / new_cell_mass;
                  }

               // PPM and DEF: TENum is total energy, GENum is internal energy
               else if (HydroMethod == PPM_DirectEuler && DualEnergyFormalism) {
                  BaryonField[TENum][index] = (new_cell_ge + new_cell_ke) / new_cell_mass;
                  BaryonField[GENum][index] = new_cell_ge / new_cell_mass;
                  }

               // PPM without DEF: TENum is total energy,
               //  internal energy calcualted by subtracting kinetic energy.
               else if (HydroMethod == PPM_DirectEuler && !DualEnergyFormalism) {
                  BaryonField[TENum][index] = (new_cell_ge + new_cell_ke) / new_cell_mass;
                  }

               }
            }
         }
      
      // Done adding feedback from this cylinder.
      // Clean up memory
      for (int j = 0; j < GridDimension[2]; j++) {
         for (int i = 0; i < GridDimension[1]; i++) {
            delete[] vol[j][i];
            }
         delete[] vol[j];
         }

      delete[] vol;

      } // End of this cylinder.

   // Add in shock feedback
   //

   // If the size of the feedback sphere is smaller than the smallest cell size,
   // inject energy into the central cells.
   //
   // feedback_radius is the radius of the feedback sphere in code units
   float feedback_radius, min_dx;

   min_dx = min(CellWidth[0][0], CellWidth[1][0]);
   min_dx = min(CellWidth[2][0], min_dx);

   feedback_radius = 2.1*min_dx;

   // Figure out the mass within the feedback radius
   float total_mass;
   float x, y, z, radius;
   int index;

   total_mass = 0.0;

   for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
         for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {
            index = k * GridDimension[1] * GridDimension[0]
                  + j * GridDimension[0] + i;

            x = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i];
            y = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j];
            z = CellLeftEdge[2][k] + 0.5 * CellWidth[2][k];

            radius = sqrt(
                          pow(x - xsink, 2.0)
                        + pow(y - ysink, 2.0)
                        + pow(z - zsink, 2.0));

            if (radius < feedback_radius) {
               total_mass += BaryonField[DensNum][index] * cell_volume;
               }
            }
         }
      }

   // Check if the energy to be injected is high enough to boost the gas above
   // the threshold temperature.
   float to_inject;
   float threshold_energy;

   to_inject = shock_heating_rate * dtFixed + tp -> StoredEnergy; //code energy

   threshold_energy = (3.0/2.0) * kboltz * AGNParticleInjectionTemperature
                      * (total_mass * MassUnits) / (Mu * mh); //cgs energy

   threshold_energy /= EnergyUnits; // code energy

   // If stored energy plus new energy is less than the threshold, store the
   // energy and mass, don't apply heating yet
   //if (to_inject < threshold_energy) {
   float myr = 1.0e6 * 3.154e+7;
   myr /= TimeUnits;

   if (Time < tp -> TimeOfLastShock + 10.0 * myr) {
      tp -> StoredEnergy += shock_heating_rate * dtFixed;
      tp -> StoredMass += mdot * dtFixed;
      }

   // If the total energy is above the threshold, inject it along with the
   // stored energy.
   else {
      tp -> TimeOfLastShock += 10.0 * myr;

      // Calculate total energy to inject

      shock_heating_rate = shock_heating_rate * dtFixed + tp -> StoredEnergy; // code energy
      shock_heating_rate /= dtFixed; // code (energy/time)

      mdot = mdot * dtFixed + tp -> StoredMass; //code mass
      mdot /= dtFixed; //code (mass/time)

      tp -> StoredEnergy = 0.0;
      tp -> StoredMass = 0.0;

      // Apply feedback
      float cell_mass;
      float cell_mass_new;
      float cell_TE, cell_GE;
      float edot_this_cell, mdot_this_cell;
   
      for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
         for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
            for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {
               index = k * GridDimension[1] * GridDimension[0]
                     + j * GridDimension[0] + i;
   
               x = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i];
               y = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j];
               z = CellLeftEdge[2][k] + 0.5 * CellWidth[2][k];
   
               radius = sqrt(
                             pow(x - xsink, 2.0)
                           + pow(y - ysink, 2.0)
                           + pow(z - zsink, 2.0));
   
               if (radius < feedback_radius) {
                  cell_mass = BaryonField[DensNum][index] * cell_volume;
                  cell_TE = BaryonField[TENum][index] * cell_mass;
   
                  if (HydroMethod == PPM_DirectEuler)
                     cell_GE = BaryonField[GENum][index];
   
                  mdot_this_cell = mdot;
                  mdot_this_cell *= cell_mass / total_mass;
                  mdot_this_cell *= dtFixed;
   
                  cell_mass_new = cell_mass + mdot_this_cell;
   
                  edot_this_cell = shock_heating_rate;
                  edot_this_cell *= cell_mass / total_mass;
                  edot_this_cell *= dtFixed;
   
                  // Update cell quantities
                  BaryonField[DensNum][index] = cell_mass_new / cell_volume;
                  
                  // ZEUS: TENum is internal energy
                  // Total energy calculated by adding kinetic
                  if (HydroMethod==Zeus_Hydro) {
                     BaryonField[TENum][index] = (cell_TE + edot_this_cell) / cell_mass_new;
                     }
   
                  // PPM and DEF: TENum is total energy, GENum is internal energy
                  else if (HydroMethod == PPM_DirectEuler &&DualEnergyFormalism) {
                     BaryonField[TENum][index] = (cell_TE + edot_this_cell) / cell_mass_new;
                     }
   
                  // PPM without DEF: TENum is total energy,
                  //  internal energy calculated by subtracting kinetic energy.
                  else if (HydroMethod == PPM_DirectEuler && !DualEnergyFormalism) {
                     BaryonField[TENum][index] = (cell_TE + edot_this_cell) / cell_mass_new;
                     BaryonField[GENum][index] = (cell_GE + edot_this_cell) / cell_mass_new;
                     }
   
                  //printf("Dens, TE: %"GSYM", %"GSYM"\n", BaryonField[DensNum][index], BaryonField[TENum][index]);
   
                  } // End if within feedback radius
   
               } // End loop over i
            } // End loop over j
         } // End loop over k

      } // End if energy is above threshold energy
   if (debug)
      printf("Leaving AGNParticleCylinderFeedback[%"ISYM"]\n", MyProcessorNumber);

   return SUCCESS;
   }

#undef DEBUG_AP
