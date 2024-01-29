/***********************************************************************
/
/ Calculates and applies AGN feedback for a jet.
/ Feedback is mass weighted- equal energy added to each particle, so that each
/ cell has the same change in velocity/temperature
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

int grid::AGNParticleMassWeightedJet(ActiveParticleType* ThisParticle, float mdot, float heating_rate) {

   /* Return if this doesn't involve us */
   if (MyProcessorNumber != ProcessorNumber) {
      if (debug)
         printf("Leaving DoAGNFeedback [%"ISYM"]\n", MyProcessorNumber);

      return SUCCESS;
      }

   // Cast this particle to an AGN particle so that we can access AGN specific
   // properties.
   ActiveParticleType_AGNParticle* tp = static_cast <ActiveParticleType_AGNParticle*>(ThisParticle);

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

   // Declare a bunch of variables
   float heating_this_cell, a, x, y, z, phi, radius_in_plane,
         edot_this_cell, kinetic_this_cell, velocity_this_cell, radius;
   float mdot_this_cell, cell_mass;

   float cell_volume = pow(CellWidth[0][0], 3.0); //Code

   // Apply feedback
   // Precompute a normalization factor for the jet
   a = (4.0 * M_PI / 3.0) * (1.0 - cos(tp -> JetAngle)) * pow(tp -> FeedbackRadius,
   3.0);
   a = 1.0 / a;

   // Subsampling uses an nxnxn cube of points.
   // If n=1, use the cell center
   int n_subsamples = 5;

   // Sum up the amount of energy being injected- should equal heating_rate!
   float injected_heating = 0.0;
   int index;

   // Keep track of the total mass
   float total_mass = 0.0;

   // Go through each cell and add up the total mass
   for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
         for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {
            index = k * GridDimension[1] * GridDimension[0]
                  + j * GridDimension[0] + i;

            // Calculate the radius of this cell relative to the particle.
            x = CellLeftEdge[0][i] + CellWidth[0][i] * 0.5;// - tp->ReturnPosition()[0];
            y = CellLeftEdge[1][j] + CellWidth[1][j] * 0.5;// - tp->ReturnPosition()[1];
            z = CellLeftEdge[2][k] + CellWidth[2][k] * 0.5;// - tp->ReturnPosition()[2];

            radius = sqrt(x*x + y*y + z*z);
            //printf("Radius=%"GSYM"\n", radius);

            // Check if the cell is possibly within the feedback radius.
            // Do this before doing subsampling.
            if (radius < tp -> FeedbackRadius + CellWidth[0][i]) {
               for (int xsub = 0; xsub < n_subsamples; xsub++) {
                  for (int ysub = 0; ysub < n_subsamples; ysub++) {
                     for (int zsub = 0; zsub < n_subsamples; zsub++) {
                        // Calculate the radius of this subsample
                        // Calculate offsets as a fraction of cell width
                        float xoff, yoff, zoff;

                        if (n_subsamples > 1) {
                           xoff = float(xsub)/float(n_subsamples - 1);
                           yoff = float(ysub)/float(n_subsamples - 1);
                           zoff = float(zsub)/float(n_subsamples - 1);

                           x = CellLeftEdge[0][i] + CellWidth[0][i] * xoff;// - tp->ReturnPosition()[0];
                           y = CellLeftEdge[1][j] + CellWidth[1][j] * yoff;// - tp->ReturnPosition()[1];
                           z = CellLeftEdge[2][k] + CellWidth[2][k] * zoff;// - tp->ReturnPosition()[2];

                           radius = sqrt(x*x + y*y + z*z);
                           }

                        // if n_subsamples == 1, only do one subsample!
                        if (n_subsamples == 1) {
                           xsub = 1;
                           ysub = 1;
                           zsub = 1;
                           }

                        // Angle is calculated relative to the jet.
                        float xhat, yhat, zhat;
                        float jxhat, jyhat, jzhat;
                        float xdotj;

                        xhat = x / radius;
                        yhat = y / radius;
                        zhat = z / radius;

                        jxhat = sin(tp -> JetTheta) * sin(tp -> JetPhi);
                        jyhat = cos(tp -> JetTheta) * sin(tp -> JetPhi);
                        jzhat = cos(tp -> JetPhi);

                        xdotj = xhat * jxhat + yhat * jyhat + zhat * jzhat;

                        // Phi is the angle beween this cell and the jet
                        phi = acos(xdotj);

                        if (phi > M_PI / 2.0) {
                           phi = M_PI - phi;
                           }

                        float magjhat = sqrt(jxhat * jxhat + jyhat * jyhat + jzhat * jzhat);

                        // Is this cell in the feedback zone?
                        if ((phi < tp -> JetAngle) && (radius < tp -> FeedbackRadius)) {
                           cell_mass = BaryonField[DensNum][index] * cell_volume; //code
                           total_mass += cell_mass / (float)(n_subsamples * n_subsamples * n_subsamples);
                           } //End feedback

                        } // End x subsample loop
                     } // End y subsample loop
                  } // End z subsample loop

               } // End within feedback radius

            } // End loop over i
         } // End loop over j
      } // End loop over k

   // Go through each cell and apply feedback
   for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
         for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {
            index = k * GridDimension[1] * GridDimension[0]
                  + j * GridDimension[0] + i;

            // Calculate the radius of this cell relative to the particle.
            x = CellLeftEdge[0][i] + CellWidth[0][i] * 0.5; //- tp->ReturnPosition()[0];
            y = CellLeftEdge[1][j] + CellWidth[1][j] * 0.5; //- tp->ReturnPosition()[1];
            z = CellLeftEdge[2][k] + CellWidth[2][k] * 0.5; //- tp->ReturnPosition()[2];

            radius = sqrt(x*x + y*y + z*z);

            // Check if the cell is possibly within the feedback radius.
            // Do this before doing subsampling.
            if (radius < tp -> FeedbackRadius + CellWidth[0][i]) {
               for (int xsub = 0; xsub < n_subsamples; xsub++) {
                  for (int ysub = 0; ysub < n_subsamples; ysub++) {
                     for (int zsub = 0; zsub < n_subsamples; zsub++) {
                        // Calculate the radius of this subsample
                        // Calculate offsets as a fraction of cell width
                        float xoff, yoff, zoff;

                        if (n_subsamples > 1) {
                           xoff = float(xsub)/float(n_subsamples - 1);
                           yoff = float(ysub)/float(n_subsamples - 1);
                           zoff = float(zsub)/float(n_subsamples - 1);

                           x = CellLeftEdge[0][i] + CellWidth[0][i] * xoff; //- tp->ReturnPosition()[0];
                           y = CellLeftEdge[1][j] + CellWidth[1][j] * yoff; //- tp->ReturnPosition()[1];
                           z = CellLeftEdge[2][k] + CellWidth[2][k] * zoff; //- tp->ReturnPosition()[2];

                           radius = sqrt(x*x + y*y + z*z);
                           }

                        // if n_subsamples == 1, only do one subsample!
                        if (n_subsamples == 1) {
                           xsub = 1;
                           ysub = 1;
                           zsub = 1;
                           }

                        // Is the radial vector of this cell aligned with the
                        // jet? If so, facing=true. Otherwise, false
	                bool facing = true;

                        // Angle is calculated relative to the jet.
                        float xhat, yhat, zhat;
                        float jxhat, jyhat, jzhat;
                        float xdotj;

                        xhat = x / radius;
                        yhat = y / radius;
                        zhat = z / radius;

                        jxhat = sin(tp -> JetTheta) * sin(tp -> JetPhi);
                        jyhat = cos(tp -> JetTheta) * sin(tp -> JetPhi);
                        jzhat = cos(tp -> JetPhi);

                        xdotj = xhat * jxhat + yhat * jyhat + zhat * jzhat;

                        // Phi is the angle beween this cell and the jet
                        phi = acos(xdotj);

                        if (phi > M_PI / 2.0) {
                           facing = false;
                           phi = M_PI - phi;
                           }

                        float magjhat = sqrt(jxhat * jxhat + jyhat * jyhat + jzhat * jzhat);

                        // Is this cell in the feedback zone?
                        if ((phi < tp -> JetAngle) && (radius < tp -> FeedbackRadius)) {
                           cell_mass = BaryonField[DensNum][index] * cell_volume; //code

                           edot_this_cell = heating_rate;
                           edot_this_cell *= cell_mass / total_mass;
                           edot_this_cell /= float(n_subsamples * n_subsamples * n_subsamples);

                           injected_heating += edot_this_cell;

                           mdot_this_cell = mdot;
                           mdot_this_cell *= cell_mass / total_mass;
                           mdot_this_cell /= float(n_subsamples * n_subsamples * n_subsamples);

                           // Calculate thermal energy to be added to this cell
                           heating_this_cell = edot_this_cell * (1.0 - tp -> KineticFraction);

                           // Add mass to this cell
                           BaryonField[DensNum][index] += mdot_this_cell * dtFixed / cell_volume;

                           // Calculate the new cell mass and add energy

                           // Add thermal energy to this cell
                           BaryonField[TENum][index] += heating_this_cell * dtFixed / cell_mass;

                           if (DualEnergyFormalism)
                              BaryonField[GENum][index] += heating_this_cell * dtFixed / cell_mass;

                           // Calculate the kinetic energy and velocity to add
                           // to this cell.
                           kinetic_this_cell = tp -> KineticFraction * edot_this_cell * dtFixed / cell_mass;
                           velocity_this_cell = sqrt(2.0 * kinetic_this_cell);

	                   if (!facing)
	                      velocity_this_cell *= -1.0;

                           // Add velocity to this cell
                           BaryonField[TENum][index] += kinetic_this_cell;

                           BaryonField[Vel1Num][index] += velocity_this_cell * sin(tp -> JetTheta) * sin(tp -> JetPhi);
                           BaryonField[Vel2Num][index] += velocity_this_cell * cos(tp -> JetTheta) * sin(tp -> JetPhi);
                           BaryonField[Vel3Num][index] += velocity_this_cell * cos(tp -> JetPhi);
                           } //End feedback

                        } // End x subsample loop
                     } // End y subsample loop
                  } // End z subsample loop

               } // End within feedback radius

            } // End loop over i
         } // End loop over j
      } // End loop over k

   // Print the expected energy input and the actual
   if (debug) {
      printf("Injected energy:: expected=%6.8"FSYM", actual=%6.8"FSYM" (%4.6"FSYM"\%)\n",
             heating_rate, injected_heating, injected_heating/heating_rate * 100.0);
      }

   if (debug)
      printf("Leaving AGNMassWeightedJet[%"ISYM"]\n", MyProcessorNumber);

   return SUCCESS;
   }

#undef DEBUG_AP
