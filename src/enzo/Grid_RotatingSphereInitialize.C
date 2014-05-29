/***********************************************************************
/
/   GRID CLASS (INITIALIZE THE GRID FOR SEDOV BLAST WAVE TEST)
/
/   written by: Greg Meece
/   date:          April 2014
/   modified1:   
/
/   PURPOSE: Sets up the grid for the RotatingSphereTestProblem. For details,
/            see Meece 2014 in ApJ. If you are using turbulence, your run directory
/            must include the file 'turbulence.in', which can be generated using the
/            file 'turbulence_generatory.py', which should be in 
/            run/Hydro/Hydro-3D/RotatingSphere
/
/            Note that setting up with turbulence can take a minute or so. If you are
/            trying to initialize with a large number of grids, the initialization might
/            take longer than you would like. There are ways around this (such as making a
/            singleton class to hold the turbulence) which would be faster but a little more
/            complex. Email me (meecegre@msu.edu) if you have questions.
/
/   RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);

// Function prototypes
float get_gas_density(float r, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift);
float get_drhodr(float r, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift);
float get_gas_temperature(float r, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature);
float get_critical_density(float redshift);

float integrate_temperature(float r1, float r2, float T1, float dr, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature);

float rk4(float (*dydt)(float r, float T, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature), float y, float t, float dt, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature);

float dTdr(float r, float T, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature);

float get_grav_accel(float r);

float* integrate_mass_energy(float r_final,
                             float core_radius,
                             float core_density,
                             float core_exponent,
                             float outer_exponent,
                             float redshift);

void mass_energy_rk4(float r,
                     float* (*derivs)(float r,
                                      float* mass_energy,
                                      float core_radius,
                                      float core_density,
                                      float core_exponent,
                                      float outer_exponent,
                                      float redshift),
                     float* mass_energy,
                     float core_radius,
                     float core_density,
                     float core_exponent,
                     float outer_exponent,
                     float redshift,
                     float dr);

float* mass_energy_derivs(float r,
                          float* mass_energy,
                          float core_radius,
                          float core_density,
                          float core_exponent,
                          float outer_exponent,
                          float redshift);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

// Physical Constants
const float VIRIAL_COEFFICIENT = 200.0; // Use r_200
const float SOLAR_MASS_IN_GRAMS = 1.9891e+33;

const float G_CGS = 6.67259e-8;
const float AMU_CGS = 1.6605402e-24;
const float KB_CGS = 1.380658e-16;
const float CM_PER_KM = 1.0e5;
const float CM_PER_MEGAPARSEC = 3.085677581e24;

// Cosmological parameters, used for computing the critical density
// and setting up the NFW halo. The setup should not be very sensitive
// to changes in these parameters.
const float HUBBLE_CONSTANT_NOW = 70.0; // km/s/Mpc
const float OMEGA_MATTER = 0.3;
const float OMEGA_LAMBDA = 0.7;

// Simulation constants
// Shouldn't need to change these unless the setup isn't working right.
const int N_RADIAL_BINS = 2048;
const int N_SAMPLES_PER_BIN = 10;
const float MIN_BIN_RADIUS = 1.0e-1; // Code units
const float MAX_BIN_RADIUS = 2.0e3; // Code units
const float RS_INTEGRATION_INTERVAL = 1.0e-2; // Code units

int grid::RotatingSphereInitializeGrid(float RotatingSphereNFWMass,
                                       float RotatingSphereNFWConcentration,
                                       float RotatingSphereCoreRadius,
                                       float RotatingSphereCentralDensity,
                                       float RotatingSphereCoreDensityExponent,
                                       float RotatingSphereOuterDensityExponent,
                                       float RotatingSphereExternalTemperature,
                                       float RotatingSphereSpinParameter,
                                       float RotatingSphereAngularMomentumExponent,
                                       int RotatingSphereUseTurbulence,
                                       float RotatingSphereTurbulenceRMS,
                                       float RotatingSphereRedshift)
{
   if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

   if(debug){
      printf("Entering RotatingSphereInitializeGrid\n");
      fflush(stdout);
      }
 
   // Figure out grid quantities and how to access fields 
   int size = 1;

   for (int dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

   int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
   int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, MetalNum;

   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                      Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
      }

   int MetallicityField = FALSE;

   if (MultiSpecies) {
      if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
         ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
         }
      }

   if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
         != -1)
      MetallicityField = TRUE;
   else
      MetalNum = 0;

   // Get the units
   float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits; 
   double MassUnits;
   
   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	    &TimeUnits, &VelocityUnits, &MassUnits, Time);

   printf("Getting ready to set up the sphere.\n");

   // Set the gravitational constant that Enzo uses.
   // This is 4 pi G_Cgs in code units.
   GravitationalConstant = 4.0 * M_PI * G_CGS * DensityUnits * pow(TimeUnits, 2.0);

   // Set up the NFW halo.
   float g_code, nfw_mvir_code, critical_density_code, c, nfw_scale_density, nfw_scale_radius;
   float nfw_v_circular, nfw_temp, nfw_cs;

   g_code = G_CGS * DensityUnits * pow(TimeUnits, 2.0);

   nfw_mvir_code = RotatingSphereNFWMass * SOLAR_MASS_IN_GRAMS / MassUnits;
   critical_density_code = get_critical_density(RotatingSphereRedshift);

   c = RotatingSphereNFWConcentration;

   nfw_scale_density = (VIRIAL_COEFFICIENT / 3.0) * critical_density_code * pow(c, 3.0) / (log(1.0+c) - c/(1.0+c));
   nfw_scale_radius = nfw_mvir_code / (4.0 * M_PI * nfw_scale_density * (log(1.0+c) - (c/(1.0+c))));
   nfw_scale_radius = pow(nfw_scale_radius, 1.0/3.0);

   // These are in code units. At the end of initialization,
   // they get changed to cgs units, since that is what is used in 
   // Grid_ComputeAccelerationFieldExternal.
   PointSourceGravity = 2;
   PointSourceGravityConstant = 4.0 * M_PI * nfw_scale_density * pow(nfw_scale_radius, 3.0) * (log(2.0) - 0.5);
   PointSourceGravityCoreRadius = nfw_scale_radius;

   nfw_v_circular = sqrt(g_code * nfw_mvir_code / (nfw_scale_radius * c));
   nfw_cs = nfw_v_circular * sqrt(Gamma/2.0);

   printf("Finished setting up the NFW halo.\n");

   /* Divide the sphere into radial bins.
    * Set the temperature, density, and rotational velocity in each one.
    * All values are in code units.
    */

   float bin_left_edge[N_RADIAL_BINS];
   float bin_right_edge[N_RADIAL_BINS];
   float bin_center[N_RADIAL_BINS];
   float bin_density[N_RADIAL_BINS];
   float bin_temperature[N_RADIAL_BINS];
   float bin_angular_velocity[N_RADIAL_BINS];
   
   // Set up the radial bins.
   float log_r_min = log10(MIN_BIN_RADIUS);
   float log_r_max = log10(MAX_BIN_RADIUS);

   float bin_step = (log_r_max - log_r_min) / (float)(N_RADIAL_BINS);
   float this_bin_log_radius = log_r_min;

   for (int i = 0; i < N_RADIAL_BINS; i++) {
      bin_left_edge[i] = pow(10.0, this_bin_log_radius);
      bin_right_edge[i] = pow(10.0, this_bin_log_radius + bin_step);
      bin_center[i] = 0.5 * (pow(10.0, this_bin_log_radius) + pow(10.0, this_bin_log_radius + bin_step));
      this_bin_log_radius += bin_step;
      }

   printf("Done setting up radial bins.\n");

   // Set the density.
   for (int i = 0; i < N_RADIAL_BINS; i++) {
      float dens_sum = 0.0;
      float rad;

      for (int s = 0; s < N_SAMPLES_PER_BIN; s++) {
         rad = bin_left_edge[i] + (bin_right_edge[i] - bin_left_edge[i]) * (float)s/((float)N_SAMPLES_PER_BIN);
         dens_sum += get_gas_density(rad,
                                     RotatingSphereCoreRadius,
                                     RotatingSphereCentralDensity,
                                     RotatingSphereCoreDensityExponent,
                                     RotatingSphereOuterDensityExponent,
                                     RotatingSphereRedshift);
         }

      bin_density[i]  = dens_sum / (float)N_SAMPLES_PER_BIN;
      }

   printf("Done setting up the density.\n");

   // Set the temperature
   for (int i = 0; i < N_RADIAL_BINS; i++) {
      float temp_sum = 0.0;
      float rad;

      for (int s = 0; s < N_SAMPLES_PER_BIN; s++) {
        rad = bin_left_edge[i] + (bin_right_edge[i] - bin_left_edge[i]) * (float)s/((float)N_SAMPLES_PER_BIN);
        temp_sum += get_gas_temperature(rad,
                                         RotatingSphereCoreRadius,
                                         RotatingSphereCentralDensity,
                                         RotatingSphereCoreDensityExponent,
                                         RotatingSphereOuterDensityExponent,
                                         RotatingSphereRedshift,
                                         RotatingSphereExternalTemperature);
         }
      bin_temperature[i] = temp_sum / (float) N_SAMPLES_PER_BIN;
      }

   printf("Done setting up the temperature.\n");

   // Set up the rotation.
   float* mass_energy;
   float r_sphere = RotatingSphereCoreRadius * pow(RotatingSphereCentralDensity / get_critical_density(RotatingSphereRedshift), 1.0 / RotatingSphereOuterDensityExponent);

   mass_energy = integrate_mass_energy(r_sphere,
                                       RotatingSphereCoreRadius,
                                       RotatingSphereCentralDensity,
                                       RotatingSphereCoreDensityExponent,
                                       RotatingSphereOuterDensityExponent,
                                       RotatingSphereRedshift);

   float total_angular_momentum = RotatingSphereSpinParameter * g_code * pow(mass_energy[0], 5.0/2.0) / pow(fabs(mass_energy[2]), 0.5);
   float total_gas_mass = mass_energy[0];

   for (int i = 0; i < N_RADIAL_BINS; i++) {
      if (bin_center[i] < r_sphere) {
         mass_energy = integrate_mass_energy(bin_center[i],
                                             RotatingSphereCoreRadius,
                                             RotatingSphereCentralDensity,
                                             RotatingSphereCoreDensityExponent,
                                             RotatingSphereOuterDensityExponent,
                                             RotatingSphereRedshift);
         bin_angular_velocity[i] = total_angular_momentum
                                   * ((RotatingSphereAngularMomentumExponent + 1.0) / total_gas_mass)
                                   * pow(mass_energy[0] / total_gas_mass, RotatingSphereAngularMomentumExponent)
                                   * pow(bin_center[i], -2.0) * (3.0/2.0);
         }
      else {
         bin_angular_velocity[i] = 0.0;
         }
      }

   printf("Done setting up rotation.\n");

   // Set up the turblence field.
   int pert_size_x;
   int pert_size_y;
   int pert_size_z;

   float*** turbulence_field_vx;
   float*** turbulence_field_vy;
   float*** turbulence_field_vz;

   if (RotatingSphereUseTurbulence) {
      FILE* inf;
      inf = fopen("turbulence.in", "r");

      if (inf == NULL) {
         printf("Could not open file 'turbulence.in', which is needed for turbulence.\n");
         printf("If it is not here, try running the file 'turbulence_generator.py'.\n");
         exit(1);
         }
      
      // Read in comments and the first two lines (the dimension and the grid size)
      int lines_read = 0;
      size_t line_length = 80;
      
      int pert_dim;

      char* line = NULL;

      char* pert_dim_x_s = new char[line_length];
      char* pert_dim_y_s = new char[line_length];
      char* pert_dim_z_s = new char[line_length];

      while (lines_read < 2) {
         getline(&line, &line_length, inf);
         
         if (line[0] != '#') {
            if (lines_read == 0) {
               sscanf(line, "%i", &pert_dim);
               lines_read ++;
               }

            else {
               sscanf(line, "%s %s %s", pert_dim_x_s, pert_dim_y_s, pert_dim_z_s);
               lines_read ++;
               }
            }
         }

      pert_size_x = atoi(pert_dim_x_s);
      pert_size_y = atoi(pert_dim_y_s);
      pert_size_z = atoi(pert_dim_z_s);

      printf("Reading in a %i dim perturbation grid of size %i x %i x %i\n", pert_dim, pert_size_x, pert_size_y, pert_size_z);

      // Create an array to hold the turbulence
      turbulence_field_vx = new float**[pert_size_z];
      turbulence_field_vy = new float**[pert_size_z];
      turbulence_field_vz = new float**[pert_size_z];

      for (int j = 0; j < pert_size_z; j++) {
         turbulence_field_vx[j] = new float*[pert_size_y];
         turbulence_field_vy[j] = new float*[pert_size_y];
         turbulence_field_vz[j] = new float*[pert_size_y];

         for (int i = 0; i < pert_size_y; i++){
            turbulence_field_vx[j][i] = new float[pert_size_x];
            turbulence_field_vy[j][i] = new float[pert_size_x];
            turbulence_field_vz[j][i] = new float[pert_size_x];
            }
         }

      for (int k = 0; k < pert_size_z; k++)
         for (int j = 0; j < pert_size_y; j++)
            for (int i = 0; i < pert_size_x; i++) {
               turbulence_field_vx[k][j][i] = 0.0;
               turbulence_field_vy[k][j][i] = 0.0;
               turbulence_field_vz[k][j][i] = 0.0;
               }

      // Read in the turbulence.
      char* x_s = new char[line_length];
      char* y_s = new char[line_length];
      char* z_s = new char[line_length];

      char* vx_s = new char[line_length];
      char* vy_s = new char[line_length];
      char* vz_s = new char[line_length];

      while (getline(&line, &line_length, inf) != -1) {
         int x, y, z;
         float vx, vy, vz;

         sscanf(line, "%s %s %s %s %s %s", x_s, y_s, z_s, vx_s, vy_s, vz_s);
         x = atoi(x_s);
         y = atoi(y_s);
         z = atoi(z_s);

         vx = atof(vx_s);
         vy = atof(vy_s);
         vz = atof(vz_s);

         turbulence_field_vx[z][y][x] = vx;
         turbulence_field_vy[z][y][x] = vy;
         turbulence_field_vz[z][y][x] = vz;
         }
      
      printf("Done reading in turbulence.\n");

      // Normalize the turbulent field so that the RMS velocity
      // is some fraction of the halo sound speed.
      float ssum = 0.0;

      for (int k = 0; k < pert_size_z; k++)
         for (int j = 0; j < pert_size_y; j++)
            for (int i = 0; i < pert_size_x; i++) {
               ssum += pow(turbulence_field_vx[k][j][i], 2.0)
                      + pow(turbulence_field_vy[k][j][i], 2.0)
                      + pow(turbulence_field_vz[k][j][i], 2.0);
               }
      ssum /= (float)(pert_size_z * pert_size_y * pert_size_x);
      ssum = sqrt(ssum);

      for (int k = 0; k < pert_size_z; k++)
         for (int j = 0; j < pert_size_y; j++)
            for (int i = 0; i < pert_size_x; i++) {
               turbulence_field_vx[k][j][i] *= (RotatingSphereTurbulenceRMS * nfw_cs / ssum);
               turbulence_field_vy[k][j][i] *= (RotatingSphereTurbulenceRMS * nfw_cs / ssum);
               turbulence_field_vz[k][j][i] *= (RotatingSphereTurbulenceRMS * nfw_cs / ssum);
               }

      // Convert the NFW variables to cgs units. This is confusing,
      // but the acceleration code in Grid_ComputeAccelerationFieldExternal
      // takes these as cgs rather than code.
      PointSourceGravityConstant *= MassUnits;
      PointSourceGravityCoreRadius *= LengthUnits;

      // Free some memory.
      delete [] x_s;
      delete [] y_s;
      delete [] z_s;

      delete [] vx_s;
      delete [] vy_s;
      delete [] vz_s;

      delete [] pert_dim_x_s;
      delete [] pert_dim_y_s;
      delete [] pert_dim_z_s;

      printf("Finished reading in turbulence.\n");
      } // End if RotatingSphereUseTurbulence

   // Set grid values.
   printf("Setting grid values. This might take a minute.\n");

   int cell_index, bin_index;
   float x, y, z;
   float radius, radius_in_plane;

   for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
         for (int i = 0; i < GridDimension[0]; i++) {
            cell_index = k * (GridDimension[1] * GridDimension[0]) + j * GridDimension[0] + i;

            // x, y, z, radius in code units
	    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	    y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	    z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

            radius = sqrt(x*x + y*y + z*z);
            radius_in_plane = sqrt(x*x + y*y);

            float sin_theta, cos_theta, sin_phi;

            sin_theta = x / radius_in_plane;
            cos_theta = y / radius_in_plane;
            sin_phi = radius_in_plane / radius;

            // Figure out which radial bin contains this cell
            bin_index = 0;

            for (int n = 1; n < N_RADIAL_BINS; n++) {
               if (bin_left_edge[n-1] < radius && bin_left_edge[n] > radius)
                  bin_index = n;
               }

            // Set the density
            BaryonField[DensNum][cell_index] = bin_density[bin_index];

            // Set the temperature
            BaryonField[TENum][cell_index] = bin_temperature[bin_index] / (TemperatureUnits * (Gamma - 1.0) * Mu);

            if (DualEnergyFormalism)
               BaryonField[GENum][cell_index] = bin_temperature[bin_index] / (TemperatureUnits * (Gamma - 1.0) * Mu);

            // Set the velocity
            BaryonField[Vel1Num][cell_index] =  (bin_angular_velocity[bin_index] * bin_center[bin_index]) * cos_theta * sin_phi;
            BaryonField[Vel2Num][cell_index] = -(bin_angular_velocity[bin_index] * bin_center[bin_index]) * sin_theta * sin_phi;
            BaryonField[Vel3Num][cell_index] = 0.0;

            // Add turbulence
            if (RotatingSphereUseTurbulence) {
               int pert_index[3];

               pert_index[0] = floor((float)pert_size_x * (x - DomainLeftEdge[0]) / (DomainRightEdge[0] - DomainLeftEdge[0]));
               pert_index[1] = floor((float)pert_size_y * (y - DomainLeftEdge[1]) / (DomainRightEdge[1] - DomainLeftEdge[1]));
               pert_index[2] = floor((float)pert_size_z * (z - DomainLeftEdge[2]) / (DomainRightEdge[2] - DomainLeftEdge[2]));

               //printf("x, y, z: %e %e %e\n", x, y, z);
               //printf("size x, size y, size z: %i %i %i\n", pert_size_x, pert_size_y, pert_size_z);

               if (pert_index[0] < 0)
                  pert_index[0] = 0;
               if (pert_index[1] < 0)
                  pert_index[1] = 0;
               if (pert_index[2] < 0)
                  pert_index[2] = 0;

               if (pert_index[0] > pert_size_x - 1)
                  pert_index[0] = pert_size_x - 1;
               if (pert_index[1] > pert_size_y - 1)
                  pert_index[1] = pert_size_y - 1;
               if (pert_index[2] > pert_size_z - 1)
                  pert_index[2] = pert_size_z - 1;

               //printf("pert x, pert y, pert z: %i %i %i\n", pert_index[0], pert_index[1], pert_index[2]);

               BaryonField[Vel1Num][cell_index] += turbulence_field_vx[pert_index[2]][pert_index[1]][pert_index[0]];
               BaryonField[Vel2Num][cell_index] += turbulence_field_vy[pert_index[2]][pert_index[1]][pert_index[0]];
               BaryonField[Vel3Num][cell_index] += turbulence_field_vz[pert_index[2]][pert_index[1]][pert_index[0]];
               }

               // Set up the chemistry.
               if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE)
                 BaryonField[MetalNum][cell_index] = BaryonField[DensNum][cell_index]*TestProblemData.MetallicityField_Fraction;

               /* set multispecies values --- EVERYWHERE, not just inside the sphere radius! */

               if(TestProblemData.MultiSpecies) {

                 BaryonField[HIINum][cell_index] = TestProblemData.HII_Fraction * 
                   TestProblemData.HydrogenFractionByMass * BaryonField[DensNum][cell_index];
                     
                 BaryonField[HeIINum][cell_index] = TestProblemData.HeII_Fraction *
                   BaryonField[DensNum][cell_index] * (1.0-TestProblemData.HydrogenFractionByMass);
                     
                 BaryonField[HeIIINum][cell_index] = TestProblemData.HeIII_Fraction *
                   BaryonField[DensNum][cell_index] * (1.0-TestProblemData.HydrogenFractionByMass);

                 BaryonField[HeINum][cell_index] = 
                   (1.0 - TestProblemData.HydrogenFractionByMass)*BaryonField[DensNum][cell_index] -
                   BaryonField[HeIINum][cell_index] - BaryonField[HeIIINum][cell_index];
                     
                 if(TestProblemData.MultiSpecies > 1){
                   BaryonField[HMNum][cell_index] = TestProblemData.HM_Fraction *
                     BaryonField[HIINum][cell_index];
               	
                   BaryonField[H2INum][cell_index] = TestProblemData.H2I_Fraction *
                     BaryonField[0][cell_index] * TestProblemData.HydrogenFractionByMass;
               	
                   BaryonField[H2IINum][cell_index] = TestProblemData.H2II_Fraction * 2.0 *
                     BaryonField[HIINum][cell_index];
                 }

                 // HI density is calculated by subtracting off the various ionized fractions
                 // from the total
                 BaryonField[HINum][cell_index] = TestProblemData.HydrogenFractionByMass*BaryonField[0][cell_index]
                   - BaryonField[HIINum][cell_index];
                 if (MultiSpecies > 1)
                   BaryonField[HINum][cell_index] -= (BaryonField[HMNum][cell_index] + BaryonField[H2IINum][cell_index]
               				      + BaryonField[H2INum][cell_index]);

                 // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
                 // density for convenience) is calculated by summing up all of the ionized species.
                 // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
                 // calculating mass density, not number density (because the BaryonField values are 4x as
                 // heavy for helium for a single electron)
                 BaryonField[DeNum][cell_index] = BaryonField[HIINum][cell_index] +
                   0.25*BaryonField[HeIINum][cell_index] + 0.5*BaryonField[HeIIINum][cell_index];
                 if (MultiSpecies > 1)
                   BaryonField[DeNum][cell_index] += 0.5*BaryonField[H2IINum][cell_index] -
                     BaryonField[HMNum][cell_index];
                     
                 // Set deuterium species (assumed to be a negligible fraction of the total, so not
                 // counted in the conservation)
                 if(TestProblemData.MultiSpecies > 2){
                   BaryonField[DINum ][cell_index]  = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HINum][cell_index];
                   BaryonField[DIINum][cell_index] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HIINum][cell_index];
                   BaryonField[HDINum][cell_index] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * BaryonField[H2INum][cell_index];
                 }

               } // if(TestProblemData.MultiSpecies)	

            } //End loop over i
         } //End loop over j
      } //End loop over k

   printf("Finished setting grid values.\n");

   // Free the turbulent velocity array
   if (RotatingSphereUseTurbulence) {
      for (int j = 0; j < pert_size_z; j++) {
         for (int i = 0; i < pert_size_y; i++) {
            delete[] turbulence_field_vx[j][i];
            delete[] turbulence_field_vy[j][i];
            delete[] turbulence_field_vz[j][i];
            }

         delete[] turbulence_field_vx[j];
         delete[] turbulence_field_vy[j];
         delete[] turbulence_field_vz[j];
         }

      delete[] turbulence_field_vx;
      delete[] turbulence_field_vy;
      delete[] turbulence_field_vz;
      }

   printf("All finished initializing this grid.\n");

   // Done with initialization
   if(debug){
      printf("Exiting RotatingSphereInitialize\n");
      fflush(stdout);
      }

   return SUCCESS;
}

/* Get the density at radius r
 *
 * Inputs:
 *   r: radius in code units
 *   core_radius: core_radius in code units
 *   core_density: core_density in code units
 *   core_exponent: The core density exponent
 *   outer_exponent: The power law exponent outside of the core
 *   redshift: The redshift
 *
 * Returns:
 *   The density at r in code units
 */
float get_gas_density(float r, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift) {
   float r_hat, dens, critical_density;

   r_hat = r / core_radius;

   dens = core_density * pow(r_hat, -core_exponent) * pow(1.0 + r_hat, core_exponent - outer_exponent);
   critical_density = get_critical_density(redshift);

   if (dens < critical_density)
      dens = critical_density;

   return dens;
}

/* Get the density derivative at radius r
 *
 * Inputs:
 *   r: radius in code units
 *   core_radius: core_radius in code units
 *   core_density: core_density in code units
 *   core_exponent: The core density exponent
 *   outer_exponent: The power law exponent outside of the core
 *   redshift: The redshift
 *
 * Returns:
 *   The density derivative d_rho/dr at r in code units
 */
float get_drhodr(float r, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift) {
   float r_hat, density, critical_density;
   r_hat = r / core_radius;

   // Need to calculate density to check if we are outside the sphere.
   density = get_gas_density(r, core_radius, core_density, core_exponent, outer_exponent, redshift);
   critical_density = get_critical_density(redshift);

   if (density > critical_density) {
      return -(core_density / core_radius)
            * pow(r_hat, -core_exponent) * pow(1.0 + r_hat, core_exponent - outer_exponent)
            * (core_exponent / r_hat + (outer_exponent - core_exponent) / (1.0 + r_hat));
      }

   else
      return 0;
}

/* Get the temperature at radius r
 *
 * Inputs:
 *   r: radius in code units
 *   core_radius: core_radius in code units
 *   core_density: core_density in code units
 *   core_exponent: The core density exponent
 *   outer_exponent: The power law exponent outside of the core
 *   redshift: The redshift
 *   exterior_temperature: The temperature outside the sphere
 *
 * Returns:
 *   The temperature in K
 */
float get_gas_temperature(float r, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature) {
   float critical_density, core_boundary_temperature, core_outer_density, density, sphere_radius;

   // Figure out the sphere radius
   critical_density = get_critical_density(redshift);
   sphere_radius = core_radius * pow(core_density / critical_density, 1.0/outer_exponent);

   // Inside the core, integrate to find the temperature for HSE.
   core_outer_density = get_gas_density(core_radius, core_radius, core_density, core_exponent, outer_exponent, redshift);
   core_boundary_temperature = exterior_temperature * pow(core_outer_density / critical_density, Gamma - 1.0);

   if (r < core_radius) {
      return integrate_temperature(core_radius, r, core_boundary_temperature, RS_INTEGRATION_INTERVAL, core_radius, core_density, core_exponent, outer_exponent, redshift, exterior_temperature);
      }

   // Outside of the core, temperature increases adiabatically
   else if (r < sphere_radius) {
      density = get_gas_density(r, core_radius, core_density, core_exponent, outer_exponent, redshift);
      return exterior_temperature * pow(density / critical_density, Gamma - 1.0);
      }

   // Outside of the sphere, temperature is constant.
   else {
      return exterior_temperature;
      }
}

/* Uses numerical integration to integrate the temperature derivative from
 * radius r1 to radius r2. Returns the temperature for HSE.
 *
 * Inputs:
 *   r1: The starting radius
 *   r2: The end radius
 *   T1: The temperature at r1
 *   dr: The radial step.
 *
 * Returns:
 *   The temperature at r2 in K
 */
float integrate_temperature(float r1, float r2, float T1, float dr, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature) {
   float T = T1;
   float r = r1;

   // Check that the input is valid
   if (r1 < r2) {
      printf("Error! r2 must be less than r1.");
      return 0;
      }

   // Do the integration
   while (r > r2) {
      T = rk4(&dTdr, T, r, -dr, core_radius, core_density, core_exponent, outer_exponent, redshift, exterior_temperature);
      r -= dr;
      }

   return T;
}

/* RK4 function for numerical integration
 * Arguments:
 *    dydt: A function pointer to the derivative function
 *    y: The dependent variable at the start of the interval y(t)
 *    t: The independent variable at the beginning of the interval
 *    dt: The step size
 *
 * Returns:
 *    y(t+dt)
 */
float rk4(float (*dydt)(float r, float T, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature), float T, float r, float dr, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature) {
   float k1, k2, k3, k4;

   k1 = (*dydt)(r, T, core_radius, core_density, core_exponent, outer_exponent, redshift, exterior_temperature);
   k2 = (*dydt)(r + 0.5*dr, T + 0.5*dr * k1, core_radius, core_density, core_exponent, outer_exponent, redshift, exterior_temperature);
   k3 = (*dydt)(r + 0.5*dr, T + 0.5*dr * k2, core_radius, core_density, core_exponent, outer_exponent, redshift, exterior_temperature);
   k4 = (*dydt)(r + dr, T + dr * k3, core_radius, core_density, core_exponent, outer_exponent, redshift, exterior_temperature);

   return T + dr * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
}

/* Radial derivative of the temperature
 *
 * Inputs:
 *   r: The radius in code units
 *   T: The temperature at that radius in K
 *   core_radius: core_radius in code units
 *   core_density: core_density in code units
 *   core_exponent: The core density exponent
 *   outer_exponent: The power law exponent outside of the core
 *   redshift: The redshift
 *   exterior_temperature: The temperature outside the sphere
 *
 * Returns:
 *   dT/dr at r
 */
float dTdr(float r, float T, float core_radius, float core_density, float core_exponent, float outer_exponent, float redshift, float exterior_temperature) {
   float density, grav_accel, drhodr;
   float amu_code, kb_code;

   // Get the units
   float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits; 
   double MassUnits;
   FLOAT Time;
   
   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	    &TimeUnits, &VelocityUnits, &MassUnits, Time);

   // Convert the atomic mass unit (in g) to code units
   amu_code = AMU_CGS / (DensityUnits * pow(LengthUnits, 3.0));

   // Convert kb (in cm^2 g s^-2 K^-1) to code units
   kb_code = KB_CGS / (pow(LengthUnits, 5.0) * DensityUnits / pow(TimeUnits, 2.0));

   density = get_gas_density(r, core_radius, core_density, core_exponent, outer_exponent, redshift);
   grav_accel = get_grav_accel(r);
   drhodr =  get_drhodr(r, core_radius, core_density, core_exponent, outer_exponent, redshift);

   return (1.0 / density) * (Mu * amu_code * density * grav_accel / kb_code - T * drhodr);
}

/* Get the gravitational acceleration at radius r due to the dark matter.
 * Currently does not include potential due to gas, which is valid if the gas
 * density is much less than the dark matter density.
 *
 * Inputs:
 *   r: The radius in code units
 *
 * Returns:
 *   The acceleration in code units.
 */
float get_grav_accel(float r) {
   float dm_mass;

   // Get the units
   float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits; 
   double MassUnits;
   FLOAT Time;
   
   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	    &TimeUnits, &VelocityUnits, &MassUnits, Time);

   // Calculate G (in cm^3 g^-1 s^-2) code units
   float g_code;
   g_code = G_CGS * DensityUnits * pow(TimeUnits, 2.0);

   dm_mass = PointSourceGravityConstant
             * (log(1.0 + r / PointSourceGravityCoreRadius) - r / (r + PointSourceGravityCoreRadius))
                / (log(2.0) - 0.5);

   return -g_code * dm_mass / pow(r, 2.0);
}

/* Computes the critical density of the universe at a given redshift.
 *
 * Inputs:
 *   redshift: The redshift
 *
 * Returns:
 *   The critical density at r in code units
 */
float get_critical_density(float redshift) {
   // Get the units
   float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits; 
   double MassUnits;
   FLOAT Time;
   
   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	    &TimeUnits, &VelocityUnits, &MassUnits, Time);

   // Calculate G (in cm^3 g^-1 s^-2) code units
   float g_code, h_code;
   g_code = G_CGS / (1.0 / (DensityUnits * pow(TimeUnits, 2.0)));

   h_code = HUBBLE_CONSTANT_NOW * (CM_PER_KM / CM_PER_MEGAPARSEC); // Now in cgs
   h_code *= TimeUnits; // Now in code units

   float E = sqrt(OMEGA_MATTER * pow(1.0 + redshift, 3.0) + OMEGA_LAMBDA);

   return 3.0 * pow(h_code * E, 2.0) / (8.0 * M_PI * g_code);
}

/* Integrates the gas density and dark matter profile to find the gas mass,
 * the total mass, and the binding energy of the gas enclosed within a given
 * radius.
 *
 * Inputs:
 *
 * Returns:
 *   An array containing the gas mass, total mass, and binding energy of the
 *   gas.
 */
float* integrate_mass_energy(float r_final,
                             float core_radius,
                             float core_density,
                             float core_exponent,
                             float outer_exponent,
                             float redshift)
{
   float dr = RS_INTEGRATION_INTERVAL;
   float r = RS_INTEGRATION_INTERVAL;

   float* mass_energy = new float[3];

   mass_energy[0] = (4.0 * M_PI * core_density) * pow(r, 3.0) * pow(r/core_radius, -core_exponent) / (3.0 - core_exponent);
   mass_energy[1] = PointSourceGravityConstant * (log((r + PointSourceGravityCoreRadius) / PointSourceGravityCoreRadius) - r / (r+PointSourceGravityCoreRadius)) / (log(2) - 0.5);
   mass_energy[2] = 0.0; // Not important

   while(r < r_final) {
      if (r + dr > r_final)
         dr = r_final - r;

      mass_energy_rk4(r, &mass_energy_derivs, mass_energy, core_radius, core_density, core_exponent, outer_exponent, redshift, dr);
      r += dr;
      }

   return mass_energy;
}

/* This is the rk4 function for the mass/energy integration.
 * Unlike the other rk4 function, this one solves a system of two
 * equations simultaniously. I'm sure that someone could come up with
 * a more compact and/or more general way to do this, but this method
 * works.
 *
 * Inputs:
 *   r: The outer limit of integration.
 *   derivs: A function which returns an array of derivatives.
 *   mass_energy: The array of gas mass, total mass, and binding energy
 *                enclosed within r.
 *   core_radius: The radius of the core in code units
 *   core_density: The scale density of the core in code units.
 *   core_exponent: The power law index in the core
 *   outer_exponent: The power law index outside of the core
 *   redshift: The redshift of the simulation
 *   dr: The integration interval
 *
 * Returns:
 *   Does not return- modifies the array in place.
 */

void mass_energy_rk4(float r,
                     float* (*derivs)(float r,
                                      float* mass_energy,
                                      float core_radius,
                                      float core_density,
                                      float core_exponent,
                                      float outer_exponent,
                                      float redshift),
                     float* mass_energy,
                     float core_radius,
                     float core_density,
                     float core_exponent,
                     float outer_exponent,
                     float redshift,
                     float dr)
{
   float temp_mass_energy[3];

   float k1[3];
   float k2[3];
   float k3[3];
   float k4[3];

   float* temp_derivs;

   temp_mass_energy[0] = mass_energy[0];
   temp_mass_energy[1] = mass_energy[1];
   temp_mass_energy[2] = mass_energy[2];

   // Calculate k1
   temp_derivs =  (*derivs)(r, temp_mass_energy, core_radius, core_density, core_exponent, outer_exponent, redshift);

   k1[0] = temp_derivs[0];
   k1[1] = temp_derivs[1];
   k1[2] = temp_derivs[2];

   delete [] temp_derivs;

   // Calculate k2
   temp_mass_energy[0] = mass_energy[0] + 0.5 * dr * k1[0];
   temp_mass_energy[1] = mass_energy[1] + 0.5 * dr * k1[1];
   temp_mass_energy[2] = mass_energy[2] + 0.5 * dr * k1[2];

   temp_derivs = (*derivs)(r + 0.5 * dr, temp_mass_energy, core_radius, core_density, core_exponent, outer_exponent, redshift);

   k2[0] = temp_derivs[0];
   k2[1] = temp_derivs[1];
   k2[2] = temp_derivs[2];

   delete [] temp_derivs;

   // Calculate k3
   temp_mass_energy[0] = mass_energy[0] + 0.5 * dr * k2[0];
   temp_mass_energy[1] = mass_energy[1] + 0.5 * dr * k2[1];
   temp_mass_energy[2] = mass_energy[2] + 0.5 * dr * k2[2];

   temp_derivs = (*derivs)(r + 0.5 * dr, temp_mass_energy, core_radius, core_density, core_exponent, outer_exponent, redshift);

   k3[0] = temp_derivs[0];
   k3[1] = temp_derivs[1];
   k3[2] = temp_derivs[2];

   delete [] temp_derivs;

   // Calculate k4
   temp_mass_energy[0] = mass_energy[0] + dr * k3[0];
   temp_mass_energy[1] = mass_energy[1] + dr * k3[1];
   temp_mass_energy[2] = mass_energy[2] + dr * k3[2];

   temp_derivs = (*derivs)(r + dr, temp_mass_energy, core_radius, core_density, core_exponent, outer_exponent, redshift);

   k4[0] = temp_derivs[0];
   k4[1] = temp_derivs[1];
   k4[2] = temp_derivs[2];

   delete [] temp_derivs;

   // Compute the new value of the masses and energies.
   mass_energy[0] += (1.0/6.0) * dr * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]);
   mass_energy[1] += (1.0/6.0) * dr * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]);
   mass_energy[2] += (1.0/6.0) * dr * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]);
}

/* Calculates the gas mass, total mass, and binding energy derivatives at radius r
 *
 * Inputs:
 *   r: The radius in code units
 *   mass_energy: Array containing gas mass, total mass, and binding energy in code units.
 *   core_radius: core_radius in code units
 *   core_density: core_density in code units
 *   core_exponent: The core density exponent
 *   outer_exponent: The power law exponent outside of the core
 *   redshift: The redshift
 *
 * Returns:
 *   An array of the derivatives at r.
 */
float* mass_energy_derivs(float r,
                          float* mass_energy,
                          float core_radius,
                          float core_density,
                          float core_exponent,
                          float outer_exponent,
                          float redshift)
{
   // Get the units
   float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits; 
   double MassUnits;
   FLOAT Time;
   
   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	    &TimeUnits, &VelocityUnits, &MassUnits, Time);

   // Calculate G (in cm^3 g^-1 s^-2) code units
   float g_code;
   g_code = G_CGS * DensityUnits * pow(TimeUnits, 2.0);

   // Calculate the dark matter density at r
   float dm_scale_density, rho_dm, r_hat;

   r_hat = r / PointSourceGravityCoreRadius;
   dm_scale_density = PointSourceGravityConstant / (4.0 * M_PI * pow(PointSourceGravityCoreRadius, 3.0) * (log(2.0) - 0.5));
   rho_dm =  dm_scale_density * pow(r_hat, -1.0) * pow(1.0 + r_hat, -2.0);

   // Calculate the derivatives
   float* derivs = new float[3];

   derivs[0] = 4.0 * M_PI * pow(r, 2.0) * get_gas_density(r, core_radius, core_density, core_exponent, outer_exponent, redshift);
   derivs[1] = 4.0 * M_PI * pow(r, 2.0) * (get_gas_density(r, core_radius, core_density, core_exponent, outer_exponent, redshift) + rho_dm);
   derivs[2] = - g_code * mass_energy[1] * 4.0 * M_PI * r * get_gas_density(r, core_radius, core_density, core_exponent, outer_exponent, redshift);

   return derivs;
}
