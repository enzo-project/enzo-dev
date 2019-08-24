/***********************************************************************
/
/   GRID CLASS (INITIALIZE THE GRID FOR THE THERMAL INSTABILITY TEST)
/
/   written by: Iryna Butsky and Cameron Hummels
/   date:         June 2018
/   modified1:   
/
/   PURPOSE: Sets up the grid for the ThermalInstability problem type.
/            Your run directory must include the file 'white_noise.in', 
/            which can be generated using the file 'white_noise_generator.py', 
/            which should be in run/Hydro/Hydro-3D/ThermalInstability
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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::ThermalInstabilityInitializeGrid(float TIMeanDensity,
                                       float TIDensityPerturbationAmplitude,
                                       float TIMeanTemperature)
{
   if (ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

   if(debug){
      printf("Entering ThermalInstabilityInitializeGrid\n");
      fflush(stdout);
      }
 
   // Figure out grid quantities and how to access fields 
   int size = 1;

   for (int dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

   int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num, PhiNum;
   int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, MetalNum;

   if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL) {
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

   printf("Getting ready to set up the thermal instability box.\n");

   // Set up the white noise field.
   int pert_size_x;
   int pert_size_y;
   int pert_size_z;

   float*** white_noise_field;

   
   FILE* inf;
   inf = fopen("white_noise.in", "r");

   if (inf == NULL) {
      printf("Could not open file 'white_noise.in', which is needed for the pressure perturbation.\n");
      printf("If it is not here, try running the file 'white_noise_generator.py' in run/ dir.\n");
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
      white_noise_field = new float**[pert_size_z];

      for (int j = 0; j < pert_size_z; j++) {
	white_noise_field[j] = new float*[pert_size_y];

         for (int i = 0; i < pert_size_y; i++){
            white_noise_field[j][i] = new float[pert_size_x];
          }
       }

      for (int k = 0; k < pert_size_z; k++)
         for (int j = 0; j < pert_size_y; j++)
            for (int i = 0; i < pert_size_x; i++)
               white_noise_field[k][j][i] = 0.0;
       

      // Read in the turbulence.
      char* x_s = new char[line_length];
      char* y_s = new char[line_length];
      char* z_s = new char[line_length];

      char* white_noise_s = new char[line_length];

      while (getline(&line, &line_length, inf) != -1) {
         int x, y, z;

         sscanf(line, "%s %s %s %s", x_s, y_s, z_s, white_noise_s);
         x = atoi(x_s);
         y = atoi(y_s);
         z = atoi(z_s);

         white_noise_field[z][y][x] = atof(white_noise_s);
         }
      
      printf("Done reading in white noise.\n");

      // Normalize the turbulent field 
      float ssum = 0.0;

      for (int k = 0; k < pert_size_z; k++)
         for (int j = 0; j < pert_size_y; j++)
            for (int i = 0; i < pert_size_x; i++) 
               ssum += pow(white_noise_field[k][j][i], 2.0);
                      
               
      ssum /= (float)(pert_size_z * pert_size_y * pert_size_x);
      ssum = sqrt(ssum);

      for (int k = 0; k < pert_size_z; k++)
         for (int j = 0; j < pert_size_y; j++)
            for (int i = 0; i < pert_size_x; i++) 
               white_noise_field[k][j][i] *= (TIDensityPerturbationAmplitude / ssum);


      // Free some memory.
      delete [] x_s;
      delete [] y_s;
      delete [] z_s;

      delete [] pert_dim_x_s;
      delete [] pert_dim_y_s;
      delete [] pert_dim_z_s;

      printf("Finished reading in white noise.\n");

   // Set grid values.
   printf("Setting grid values. This might take a minute.\n");

   int cell_index;
   float x, y, z;
   float NormalizedPerturbation, PerturbedDensity, PerturbedTemperature;

   for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
         for (int i = 0; i < GridDimension[0]; i++) {
            cell_index = k * (GridDimension[1] * GridDimension[0]) + j * GridDimension[0] + i;

            // x, y, z, radius in code units
	    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	    y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	    z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

            int pert_index[3];

            pert_index[0] = floor((float)pert_size_x * (x - DomainLeftEdge[0]) / (DomainRightEdge[0] - DomainLeftEdge[0]));
            pert_index[1] = floor((float)pert_size_y * (y - DomainLeftEdge[1]) / (DomainRightEdge[1] - DomainLeftEdge[1]));
            pert_index[2] = floor((float)pert_size_z * (z - DomainLeftEdge[2]) / (DomainRightEdge[2] - DomainLeftEdge[2]));

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

	    NormalizedPerturbation = (1.0 + white_noise_field[pert_index[2]][pert_index[1]][pert_index[0]]);

	    PerturbedDensity = (TIMeanDensity / DensityUnits) * NormalizedPerturbation; 
	    PerturbedTemperature = TIMeanTemperature / TemperatureUnits / NormalizedPerturbation;

            // Set the density                                                                                                                                   
            BaryonField[DensNum][cell_index] = PerturbedDensity;

            // Set the temperature 
            BaryonField[TENum][cell_index] = PerturbedTemperature / ((Gamma - 1.0) * Mu);

            if (DualEnergyFormalism)
	      BaryonField[GENum][cell_index] =  PerturbedTemperature / ((Gamma - 1.0) * Mu);

            // Set the velocity                                                                                                  
            BaryonField[Vel1Num][cell_index] = 0.0;
            BaryonField[Vel2Num][cell_index] = 0.0;
            BaryonField[Vel3Num][cell_index] = 0.0;       
	    
	    if (HydroMethod == MHD_RK){
	      BaryonField[B1Num][cell_index] = 0.0; 
	      BaryonField[B2Num][cell_index] = 0.0; 
	      BaryonField[B3Num][cell_index] = 0.0; 
	      BaryonField[PhiNum][cell_index] = 0.0; 
	    }

            // Set up the chemistry.
            if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE)
               BaryonField[MetalNum][cell_index] = BaryonField[DensNum][cell_index]*TestProblemData.MetallicityField_Fraction;

            /* set multispecies values */

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

   // Free the white noise array
   for (int j = 0; j < pert_size_z; j++) {
      for (int i = 0; i < pert_size_y; i++) {
         delete[] white_noise_field[j][i];
      }
      delete[] white_noise_field[j];
   }
   delete[] white_noise_field;


   printf("All finished initializing this grid.\n");

   // Done with initialization
   if(debug){
      printf("Exiting ThermalInstabilityInitialize\n");
      fflush(stdout);
      }

   return SUCCESS;
}

