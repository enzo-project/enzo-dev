/*----------------------------------------------------------------------------
/
/ INDIVIDUAL STAR PARTICLE PROPERTIES INITIALIZER
/
/ written by: Andrew Emerick
/ date:       March, 2016
/ modified1:
/
/ Purpose:
/
/    Initializes look-up table for individual star ionizing radiation
/    rates based on tables calculated using TLUSTY from Ostar2002
/    (Lanz & Hubeny 2003)
/
----------------------------------------------------------------------------*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "flowdefs.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "communication.h"
#include "CommunicationUtilities.h"



#include "StarParticleData.h"
#include "IndividualStarProperties.h"

int IndividualStarProperties_Initialize(TopGridData &MetaData);
int IndividualStarRadiationProperties_Initialize(void);


int IndividualStarProperties_Initialize(TopGridData &MetaData){

  if (IndividualStarPropertiesData.Teff != NULL && IndividualStarPropertiesData.R != NULL){
    return SUCCESS; // already initialized
  }

  IndividualStarPropertiesData.NumberOfMassBins        = 0;
  IndividualStarPropertiesData.NumberOfMetallicityBins = 0;

  IndividualStarPropertiesData.Zsolar = 0.01524; // see IndividualStarData.h

  if (IndividualStarOutputChemicalTags){
    IndividualStarChemicalTagFilename = new char[80];
    sprintf(IndividualStarChemicalTagFilename, "star_abundances%4.4d.dat",
                 (Eint32) MetaData.DataDumpNumber-1);
  }

  // read in the data to populate the tables
  FILE *fptr = fopen("parsec_zams.in", "r");
  if (fptr == NULL){
    ENZO_FAIL("Error opening stellar properties file, 'parsec_zams.in' \n");
  }

  char line[MAX_LINE_LENGTH];
  IndividualStarPropertiesData.NumberOfMassBins        = 26; // hard code for now
  IndividualStarPropertiesData.NumberOfMetallicityBins = 11;

  /* store bin values (not evenly spaced) */
  IndividualStarPropertiesData.M = new float[IndividualStarPropertiesData.NumberOfMassBins];
  IndividualStarPropertiesData.Z = new float[IndividualStarPropertiesData.NumberOfMetallicityBins];

  /* arrays for data */
  IndividualStarPropertiesData.Teff = new float*[IndividualStarPropertiesData.NumberOfMassBins];
  IndividualStarPropertiesData.R    = new float*[IndividualStarPropertiesData.NumberOfMassBins];
  IndividualStarPropertiesData.L    = new float*[IndividualStarPropertiesData.NumberOfMassBins];
  IndividualStarPropertiesData.lifetime    = new float*[IndividualStarPropertiesData.NumberOfMassBins];
  IndividualStarPropertiesData.agb_start = new float*[IndividualStarPropertiesData.NumberOfMassBins];

  /* fill arrays in each dimension */
  for (int i = 0; i < IndividualStarPropertiesData.NumberOfMassBins; i++){
    IndividualStarPropertiesData.Teff[i] = new float[IndividualStarPropertiesData.NumberOfMetallicityBins];
    IndividualStarPropertiesData.R[i]    = new float[IndividualStarPropertiesData.NumberOfMetallicityBins];
    IndividualStarPropertiesData.L[i]    = new float[IndividualStarPropertiesData.NumberOfMetallicityBins];
    IndividualStarPropertiesData.lifetime[i]    = new float[IndividualStarPropertiesData.NumberOfMetallicityBins];
    IndividualStarPropertiesData.agb_start[i] = new float[IndividualStarPropertiesData.NumberOfMetallicityBins];

    /* need to do NULL things if I move away from hard coding sizes above */
  }

  /* read in the data */
  int i, j, err;
  float L, Teff, R, lifetime, agb_start;
  i = 0; j = 0;
  while( fgets(line, MAX_LINE_LENGTH, fptr) != NULL){

    if(line[0] != '#'){
      err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                         &IndividualStarPropertiesData.M[i],
                         &IndividualStarPropertiesData.Z[j],
                         &L,
                         &Teff,
                         &R,
                         &lifetime,
                         &agb_start);
      // un-log the data
      IndividualStarPropertiesData.L[i][j]           = POW(10.0, L);
      IndividualStarPropertiesData.Teff[i][j]        = POW(10.0, Teff);
      IndividualStarPropertiesData.R[i][j]           = POW(10.0, R);
      IndividualStarPropertiesData.lifetime[i][j]    = lifetime;
      IndividualStarPropertiesData.agb_start[i][j] = agb_start;

      j++;
      if ( j >= IndividualStarPropertiesData.NumberOfMetallicityBins){
        j = 0;
        i++;
      }

    } // end #
  } // end while

  fclose(fptr);

  return SUCCESS;
}

int IndividualStarRadiationProperties_Initialize(void){
  /*
     Load and Initialize the radiative properties of the individual stars. This
     currently uses hard-coded assumptions about the tables that they will
     be interpolating over effective temperature and surface gravity based on
     the OSTAR2002 models of stellar SEDs. Constructing these tables is
     detailed a little more in:

     https://github.com/aemerick/stars/blob/master/radiation/data/extract_sed.py

     and

     https://github.com/aemerick/stars/blob/master/radiation/integrate_SED.py

     NOTE: There was a non-backwards compatable change to this made in April 2020
           with updated yield tables, changing column ordering, changing
           files from storing logged to un-logged values, and adding new rates.
           Input file names were changed to ensure old / new files weren't used
           with the wrong code.
  */

  if (IndividualStarRadData.q0 != NULL && IndividualStarRadData.q1 != NULL){
    return SUCCESS;
  }

  IndividualStarRadData.NumberOfTemperatureBins = 0;
  IndividualStarRadData.NumberOfSGBins          = 0;
  IndividualStarRadData.NumberOfMetallicityBins = 0;

  IndividualStarRadData.Zsolar = 0.017; // see IndividualStarData.h
  // now read in the data to populate the tables

  // read in from file
  FILE *fptrq0 = fopen("q0_photon_rates.in", "r");
  if (fptrq0 == NULL){
    ENZO_FAIL("Error opening q0_photon_rates.in\n");
  }

  char line[MAX_LINE_LENGTH];
/*
  NEED TO PUT CODE IN HERE TO FIND NUMBER OF BINS IN EACH DIMENSION,
  or have user input dimensions by hand (not ideal)...
  For now, just going to hard code the sizes rather than worrying about
  how to do this... unlikely variable sizes will be needed anyway...
*/
  IndividualStarRadData.NumberOfTemperatureBins = INDIVIDUAL_STAR_TEMPERATURE_BINS;
  IndividualStarRadData.NumberOfMetallicityBins = INDIVIDUAL_STAR_METALLICITY_BINS;
  IndividualStarRadData.NumberOfSGBins = INDIVIDUAL_STAR_SG_BINS;

  /* bin values - not evenly spaced, cannot just store spacing and start */
  IndividualStarRadData.T    = new float[IndividualStarRadData.NumberOfTemperatureBins];
  IndividualStarRadData.g    = new float[IndividualStarRadData.NumberOfSGBins];
  IndividualStarRadData.Z    = new float[IndividualStarRadData.NumberOfMetallicityBins];

  /* Code here to read in metallicity bins from table */

  // Metallicities in table are listed from low to high
  // first column in table is z = 0.0, last column is z = 2.0
  IndividualStarRadData.Z[9] = 2.0;
  IndividualStarRadData.Z[8] = 1.0;
  IndividualStarRadData.Z[7] = 0.5;
  IndividualStarRadData.Z[6] = 0.2;
  IndividualStarRadData.Z[5] = 0.1;
  IndividualStarRadData.Z[4] = 1.0/30.0;
  IndividualStarRadData.Z[3] = 0.02;    // 1 / 50.0
  IndividualStarRadData.Z[2] = 0.01;
  IndividualStarRadData.Z[1] = 0.001;
  IndividualStarRadData.Z[0] = 0.0;

  // Load the arrays for ionizing radiation and the other bands
  // q's correspond to HI, HeI, and HeII photon fluxes (count / s / cm^2)
  // while band fluxes are energy fluxes (erg / s / cm^2)
  IndividualStarRadData.q0       = new float**[IndividualStarRadData.NumberOfTemperatureBins];
  IndividualStarRadData.q1       = new float**[IndividualStarRadData.NumberOfTemperatureBins];
  IndividualStarRadData.q2       = new float**[IndividualStarRadData.NumberOfTemperatureBins];
  IndividualStarRadData.IR_flux  = new float**[IndividualStarRadData.NumberOfTemperatureBins];
  IndividualStarRadData.FUV_flux = new float**[IndividualStarRadData.NumberOfTemperatureBins];
  IndividualStarRadData.LW_flux  = new float**[IndividualStarRadData.NumberOfTemperatureBins];

  // fill the arrays in all dimensions
  for(int i = 0; i < IndividualStarRadData.NumberOfTemperatureBins; i++){
    IndividualStarRadData.q0[i] = new float*[IndividualStarRadData.NumberOfSGBins];
    IndividualStarRadData.q1[i] = new float*[IndividualStarRadData.NumberOfSGBins];
    IndividualStarRadData.q2[i] = new float*[IndividualStarRadData.NumberOfSGBins];
    IndividualStarRadData.IR_flux[i]  = new float*[IndividualStarRadData.NumberOfSGBins];
    IndividualStarRadData.FUV_flux[i] = new float*[IndividualStarRadData.NumberOfSGBins];
    IndividualStarRadData.LW_flux[i]  = new float*[IndividualStarRadData.NumberOfSGBins];

    for ( int j = 0; j < IndividualStarRadData.NumberOfSGBins; j ++){
      IndividualStarRadData.q0[i][j] = new float[IndividualStarRadData.NumberOfMetallicityBins];
      IndividualStarRadData.q1[i][j] = new float[IndividualStarRadData.NumberOfMetallicityBins];
      IndividualStarRadData.q2[i][j] = new float[IndividualStarRadData.NumberOfMetallicityBins];
      IndividualStarRadData.IR_flux[i][j]  = new float[IndividualStarRadData.NumberOfMetallicityBins];
      IndividualStarRadData.FUV_flux[i][j] = new float[IndividualStarRadData.NumberOfMetallicityBins];
      IndividualStarRadData.LW_flux[i][j]  = new float[IndividualStarRadData.NumberOfMetallicityBins];

  /* AJE 3/16/16: Need to NULL everything if allow for variable size */
  //    for ( int k = 0; k < INDIVIDUAL_STAR_METALLICITY_BINS; k ++){
  //      IndividualStarq0Data[i][j][k] = NULL;
  //      IndividualStarq1Data[i][j][k] = NULL;
  //    }
    }
  } // loop over to initialize

  /* rewind if bin sizes are read in from file */
  // rewind(fptrq0); // rewind before reading in data

  /* read in ionizing flux q0 */
  int i, j, err;
  i = 0; j = 0;
  while( fgets(line, MAX_LINE_LENGTH, fptrq0) != NULL){
    if(line[0] != '#'){

      err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                   &IndividualStarRadData.T[i],
                   &IndividualStarRadData.g[j],
                   &IndividualStarRadData.q0[i][j][0],
                   &IndividualStarRadData.q0[i][j][1],
                   &IndividualStarRadData.q0[i][j][2],
                   &IndividualStarRadData.q0[i][j][3],
                   &IndividualStarRadData.q0[i][j][4],
                   &IndividualStarRadData.q0[i][j][5],
                   &IndividualStarRadData.q0[i][j][6],
                   &IndividualStarRadData.q0[i][j][7],
                   &IndividualStarRadData.q0[i][j][8],
                   &IndividualStarRadData.q0[i][j][9]);

      // increment i and j
      j++;
      if (j >= IndividualStarRadData.NumberOfSGBins){
        j = 0;
        i++;
      }
    } // check #
  } // end while read

  fclose(fptrq0);

  /* Now do the same for the q1 data */
  FILE *fptrq1 = fopen("q1_photon_rates.in", "r");
  if (fptrq1 == NULL){
    ENZO_FAIL("Error opening q1_photon_rates.in\n");
  }

  i = 0; j = 0;
  while( fgets(line, MAX_LINE_LENGTH, fptrq1) != NULL){
    float temp1, temp2;

    if(line[0] != '#'){

      err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                   &temp1,
                   &temp2,
                   &IndividualStarRadData.q1[i][j][0],
                   &IndividualStarRadData.q1[i][j][1],
                   &IndividualStarRadData.q1[i][j][2],
                   &IndividualStarRadData.q1[i][j][3],
                   &IndividualStarRadData.q1[i][j][4],
                   &IndividualStarRadData.q1[i][j][5],
                   &IndividualStarRadData.q1[i][j][6],
                   &IndividualStarRadData.q1[i][j][7],
                   &IndividualStarRadData.q1[i][j][8],
                   &IndividualStarRadData.q1[i][j][9]);

      j++;
      if (j >= IndividualStarRadData.NumberOfSGBins){
        j = 0;
        i++;
      }
    } // check #
  } // end while read

  fclose(fptrq1);

  /* Now do the same for the q2 data */
  FILE *fptrq2 = fopen("q2_photon_rates.in", "r");
  if (fptrq2 == NULL){
    ENZO_FAIL("Error opening q2_photon_rates.in\n");
  }

  i = 0; j = 0;
  while( fgets(line, MAX_LINE_LENGTH, fptrq2) != NULL){
    float temp1, temp2;

    if(line[0] != '#'){

      err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                   &temp1,
                   &temp2,
                   &IndividualStarRadData.q2[i][j][0],
                   &IndividualStarRadData.q2[i][j][1],
                   &IndividualStarRadData.q2[i][j][2],
                   &IndividualStarRadData.q2[i][j][3],
                   &IndividualStarRadData.q2[i][j][4],
                   &IndividualStarRadData.q2[i][j][5],
                   &IndividualStarRadData.q2[i][j][6],
                   &IndividualStarRadData.q2[i][j][7],
                   &IndividualStarRadData.q2[i][j][8],
                   &IndividualStarRadData.q2[i][j][9]);

      j++;
      if (j >= IndividualStarRadData.NumberOfSGBins){
        j = 0;
        i++;
      }
    } // check #
  } // end while read

  fclose(fptrq2);

  if(IndividualStarIRRadiation){
    FILE *fptrIR = fopen("IR_energy_rates.in","r");
    if (fptrIR == NULL){
      ENZO_FAIL("Error opening IR_energy_rates.in");
    }


    /* read in the IR data*/
    i = 0; j = 0;
    while( fgets(line, MAX_LINE_LENGTH, fptrIR) != NULL){
      float temp1, temp2;

      if(line[0] != '#'){

        err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                     &temp1,
                     &temp2,
                     &IndividualStarRadData.IR_flux[i][j][0],
                     &IndividualStarRadData.IR_flux[i][j][1],
                     &IndividualStarRadData.IR_flux[i][j][2],
                     &IndividualStarRadData.IR_flux[i][j][3],
                     &IndividualStarRadData.IR_flux[i][j][4],
                     &IndividualStarRadData.IR_flux[i][j][5],
                     &IndividualStarRadData.IR_flux[i][j][6],
                     &IndividualStarRadData.IR_flux[i][j][7],
                     &IndividualStarRadData.IR_flux[i][j][8],
                     &IndividualStarRadData.IR_flux[i][j][9]);


        j++;
        if (j >= IndividualStarRadData.NumberOfSGBins){
          j = 0;
          i++;
        }
      }// check #
    }// end while read

    fclose(fptrIR);

    } else{
      for(i = 0; i < IndividualStarRadData.NumberOfTemperatureBins; i++){
        for(j = 0; j < IndividualStarRadData.NumberOfSGBins; j++){
          for(int k = 0; k < IndividualStarRadData.NumberOfMetallicityBins; k++){
            IndividualStarRadData.IR_flux[i][j][k] = 0.0;
          }
        }
      }
    }

  if(IndividualStarFUVHeating){
    FILE *fptrFuv = fopen("FUV_energy_rates.in", "r");
    if (fptrFuv == NULL){
      ENZO_FAIL("Error opening FUV_energy_rates.in");
    }

    /* read in the FUV data*/
    i = 0; j = 0;
    while( fgets(line, MAX_LINE_LENGTH, fptrFuv) != NULL){
      float temp1, temp2;

      if(line[0] != '#'){

        // note metallicities are currently listed in CORRECT order for FUV file

        err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                     &temp1,
                     &temp2,
                     &IndividualStarRadData.FUV_flux[i][j][0],
                     &IndividualStarRadData.FUV_flux[i][j][1],
                     &IndividualStarRadData.FUV_flux[i][j][2],
                     &IndividualStarRadData.FUV_flux[i][j][3],
                     &IndividualStarRadData.FUV_flux[i][j][4],
                     &IndividualStarRadData.FUV_flux[i][j][5],
                     &IndividualStarRadData.FUV_flux[i][j][6],
                     &IndividualStarRadData.FUV_flux[i][j][7],
                     &IndividualStarRadData.FUV_flux[i][j][8],
                     &IndividualStarRadData.FUV_flux[i][j][9]);


        j++;
        if (j >= IndividualStarRadData.NumberOfSGBins){
          j = 0;
          i++;
        }
      } // check #
    } // end while read
    fclose(fptrFuv);

  } else{
    for(i = 0; i < IndividualStarRadData.NumberOfTemperatureBins; i++){
      for(j = 0; j < IndividualStarRadData.NumberOfSGBins; j++){
        for(int k = 0; k < IndividualStarRadData.NumberOfMetallicityBins; k++){
          IndividualStarRadData.FUV_flux[i][j][k] = 0.0;
        }
      }
    }
  }

  if(IndividualStarLWRadiation){
    FILE *fptrLW = fopen("LW_energy_rates.in", "r");
    if (fptrLW == NULL){
      ENZO_FAIL("Error opening LW_energy_rates.in");
    }

    /* read in LW data */
    i = 0; j = 0;
    while( fgets(line, MAX_LINE_LENGTH, fptrLW) != NULL){
      float temp1, temp2;

      if(line[0] != '#'){
        // note metallicities are currently listed in CORRECT order for FUV file

        err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                     &temp1,
                     &temp2,
                     &IndividualStarRadData.LW_flux[i][j][0],
                     &IndividualStarRadData.LW_flux[i][j][1],
                     &IndividualStarRadData.LW_flux[i][j][2],
                     &IndividualStarRadData.LW_flux[i][j][3],
                     &IndividualStarRadData.LW_flux[i][j][4],
                     &IndividualStarRadData.LW_flux[i][j][5],
                     &IndividualStarRadData.LW_flux[i][j][6],
                     &IndividualStarRadData.LW_flux[i][j][7],
                     &IndividualStarRadData.LW_flux[i][j][8],
                     &IndividualStarRadData.LW_flux[i][j][9]);


        j++;
        if (j >= IndividualStarRadData.NumberOfSGBins){
          j = 0;
          i++;
        }
      } // check #


    }

    fclose(fptrLW);
  } else{
    for(i = 0; i < IndividualStarRadData.NumberOfTemperatureBins; i++){
      for(j = 0; j < IndividualStarRadData.NumberOfSGBins; j++){
        for(int k = 0; k < IndividualStarRadData.NumberOfMetallicityBins; k++){
          IndividualStarRadData.LW_flux[i][j][k] = 0.0;
        }
      }
    }
  }

  // convert surface gravity values to linear scale
  for( j = 0; j < IndividualStarRadData.NumberOfSGBins; j++){
    IndividualStarRadData.g[j] = POW(10.0, IndividualStarRadData.g[j]);
  }


  return SUCCESS;
}
