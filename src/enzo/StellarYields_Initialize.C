/*****************************************************************************
/
/ INITIALIZES STELLAR YIELD TABLES
/
/ written by: Andrew Emerick
/ date:       March, 2016
/ modified1:
/
/ PURPOSE: Init stellar yield tables if ON
/
/ RETURNS: SUCCESS or FAIL
/
*****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h> // maybe don't need
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "StarParticleData.h"
#include "StellarYieldsRoutines.h"

/* function prototypes */
void unpack_line_to_yields( char *line, float *dummy);
void initialize_table(StellarYieldsDataType* table);
void fill_table(StellarYieldsDataType *table, FILE *fptr);


int ChemicalSpeciesBaryonFieldNumber(const int &atomic_number);
char* ChemicalSpeciesBaryonFieldLabelByFieldType(const int &field_num);



int InitializeStellarYieldFields(HierarchyEntry &TopGrid,
                                 TopGridData &MetaData,
                                 ExternalBoundary &Exterior,
                                 LevelHierarchyEntry *LevelArray[]){
// Initializes species yields if they do not already exist


  if ( !IndividualStarFollowStellarYields ||
       !(TestProblemData.MultiMetals || MultiMetals)       ||
       !STARMAKE_METHOD(INDIVIDUAL_STAR)){
    return SUCCESS;
  }

  int OldNumberOfBaryonFields = 0, FieldsToAdd = 0;
  int TypesToAdd[MAX_NUMBER_OF_BARYON_FIELDS];
  int ExistingTypes[MAX_NUMBER_OF_BARYON_FIELDS];

  for (int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i ++)
    ExistingTypes[i] = FieldUndefined;

  for(int yield_i = 0; yield_i < StellarYieldsNumberOfSpecies; yield_i++){
    //if(StellarYieldsAtomicNumbers[yield_i] > 2){
      TypesToAdd[FieldsToAdd++] =
                     ChemicalSpeciesBaryonFieldNumber(StellarYieldsAtomicNumbers[yield_i]);
    //}
  } // loop over tracer fields to add

  for (int i = FieldsToAdd; i < MAX_NUMBER_OF_BARYON_FIELDS; i++){
    TypesToAdd[i] = FieldUndefined;
  }

  /* Check if the fields already exist */
  OldNumberOfBaryonFields = LevelArray[0]->GridData->
    ReturnNumberOfBaryonFields();
  LevelArray[0]->GridData->ReturnFieldType(ExistingTypes);

  for (int i = 0; i < FieldsToAdd; i++){
    for (int j = 0; j < OldNumberOfBaryonFields; j++){
      if(TypesToAdd[i] == ExistingTypes[j]) {

        for (int k = i; k < FieldsToAdd; k++){
          TypesToAdd[k] = TypesToAdd[k+1];
        }
        i--;

        break;
      } // endif
    } // end oldnumberofbaryonfields loop
  } // end fields to add loop

  FieldsToAdd = 0;
  while (TypesToAdd[FieldsToAdd] != FieldUndefined)
    FieldsToAdd++;

  // Add the fields
  if (FieldsToAdd > 0 && debug)
    fprintf(stdout, "InitializeStellarYieldsFields: Increasing baryon fields "
             "from %"ISYM" to %"ISYM"\n", OldNumberOfBaryonFields,
              OldNumberOfBaryonFields + FieldsToAdd);

  // Add an extra one?? (copied over from RT, but do I actually need the +1?)
  if (OldNumberOfBaryonFields+FieldsToAdd+1 > MAX_NUMBER_OF_BARYON_FIELDS)
    ENZO_FAIL("Exceeds MAX_NUMBER_OF_BARYON_FIELDS. Please increase and re-compile.");

  LevelHierarchyEntry *Temp;

  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++){
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel){
      Temp->GridData->AddFields(TypesToAdd, FieldsToAdd);
    }
  }

  // Add external boundaries
  for (int i = 0; i < FieldsToAdd; i++){
    Exterior.AddField(TypesToAdd[i]);
  }

  for (int i = 0; i < FieldsToAdd; i ++){
   //if(StellarYieldsAtomicNumbers[i] > 2){
     DataLabel[OldNumberOfBaryonFields+i] =\
          ChemicalSpeciesBaryonFieldLabelByFieldType(TypesToAdd[i]);
   //}
  }

  return SUCCESS;
}

int InitializeStellarYields(const float &time){
  /* ------------------------------------------------------
   * InitializeStellarYields
   * -------------------------------------------------------
   * A. Emerick - April 2016
   *
   * Initializes stellar yields lookup tables if requested
   * -------------------------------------------------------*/

  // Make sure we have the right parameters set to do this
  //
  // Requires:
  //           MultiMetals == 2
  //           ChemicalEjecta ON
  // Useless unless: (as of May 2016)
  //           IndividualStar SF method
  if( !IndividualStarFollowStellarYields ||
      !TestProblemData.MultiMetals       ||
      !STARMAKE_METHOD(INDIVIDUAL_STAR)) {
    return SUCCESS;

  } else if (IndividualStarFollowStellarYields && !TestProblemData.MultiMetals){
    printf("Failure in InitializeStellarYields. MultiMetals must be enabled to follow yields\n");
    return FAIL;
  }

  if (StellarYieldsSNData.M != NULL && StellarYieldsSNData.Z != NULL){
    return SUCCESS; // already initialized
  }

  // AJE: Hard code he number of bins for now
  //     - I want to fix this but this is not priority -
  //     - unfixed as of April 2016 -
  StellarYieldsSNData.NumberOfMassBins        = 12;
  StellarYieldsSNData.NumberOfMetallicityBins =  5;
  StellarYieldsSNData.NumberOfYields          = StellarYieldsNumberOfSpecies;

  StellarYieldsWindData.NumberOfMassBins        = 12;
  StellarYieldsWindData.NumberOfMetallicityBins =  5;
  StellarYieldsWindData.NumberOfYields          = StellarYieldsNumberOfSpecies;

  StellarYieldsMassiveStarData.NumberOfMassBins        = 30;
  StellarYieldsMassiveStarData.NumberOfMetallicityBins = 12;
  StellarYieldsMassiveStarData.NumberOfYields       = StellarYieldsNumberOfSpecies;

  StellarYieldsPopIIIData.NumberOfMassBins = 7 + 14 ;  //  7 from Nomoto+2016 for Type II (13 < M < 40)
                                                       // 14 from Heger+Woosley2002 for PISN (60 < M < 130)
  StellarYieldsPopIIIData.NumberOfMetallicityBins = 1;
  StellarYieldsPopIIIData.NumberOfYields = StellarYieldsNumberOfSpecies;

  // read in data from files - one table for each yield type:
  //   1) core collapse supernova
  //   2) stellar winds
  //   3) Massive star tables
  //
  FILE *fptr_sn = fopen("stellar_yields_sn.in", "r");
  if (fptr_sn == NULL){
    ENZO_FAIL("Error opening stellar yields SN file, 'stellar_yields_sn.in");
  }
  FILE *fptr_wind = fopen("stellar_yields_wind.in", "r");
  if (fptr_wind == NULL){
    ENZO_FAIL("Error opening stellar yields wind file, 'stellar_yields_wind.in'");
  }

  FILE *fptr_mstar = fopen("stellar_yields_massive_star.in", "r");
  if (fptr_mstar == NULL){
    ENZO_FAIL("Error opening stellar yields massive stars, 'stellar_yields_massive_star.in'");
  }

  /* Initialize tables with empty pointers */
  initialize_table(&StellarYieldsSNData);
  initialize_table(&StellarYieldsWindData);
  initialize_table(&StellarYieldsMassiveStarData);

  /* Now fill the tables with data from respective files */
  fill_table(&StellarYieldsSNData, fptr_sn);
  fill_table(&StellarYieldsWindData, fptr_wind);
  fill_table(&StellarYieldsMassiveStarData, fptr_mstar);

  /* close files */
  fclose(fptr_sn);
  fclose(fptr_wind);
  fclose(fptr_mstar);

  if (IndividualStarPopIIIFormation){

    FILE *fptr_popIII = fopen("popIII_yields.in", "r");

    if (fptr_popIII == NULL){
      ENZO_FAIL("Error opening stellar yields for pop III stars, 'popIII_yields.in'");
    }

    initialize_table(&StellarYieldsPopIIIData);
    fill_table(&StellarYieldsPopIIIData, fptr_popIII);
    fclose(fptr_popIII);
  }

  /* If we are doing artificial injection events */
  if (MetalMixingExperiment) {

   MixingExperimentData.NumberOfEvents = 0;

   MixingExperimentData.time = new float[MAX_TIME_ACTIONS];

   MixingExperimentData.xpos = new float[MAX_TIME_ACTIONS];
   MixingExperimentData.ypos = new float[MAX_TIME_ACTIONS];
   MixingExperimentData.zpos = new float[MAX_TIME_ACTIONS];

   MixingExperimentData.M_ej = new float[MAX_TIME_ACTIONS];
   MixingExperimentData.E_ej = new float[MAX_TIME_ACTIONS];

   MixingExperimentData.anums = new int[StellarYieldsNumberOfSpecies];
   MixingExperimentData.yield = new float*[MAX_TIME_ACTIONS];

   /* Zero everything */
   for (int i = 0; i < MAX_TIME_ACTIONS; i++){

     MixingExperimentData.time[i] = -1.0;

     MixingExperimentData.xpos[i] = -1.0;
     MixingExperimentData.ypos[i] = -1.0;
     MixingExperimentData.zpos[i] = -1.0;

     MixingExperimentData.M_ej[i] = 0.0;
     MixingExperimentData.E_ej[i] = 0.0;

     MixingExperimentData.yield[i] = new float[StellarYieldsNumberOfSpecies];
     for (int j = 0; j < StellarYieldsNumberOfSpecies; j ++){
       MixingExperimentData.yield[i][j] = 0.0; // MASS (not mass fraction of event)
     }
   }

   for (int j = 0; j < StellarYieldsNumberOfSpecies; j ++){
     MixingExperimentData.anums[j] = -1;
   }

    FILE *fptr_mix = fopen("mixing_events.in", "r");
    if (fptr_mix == NULL){
      ENZO_FAIL("Error opening metal mixing experiment events file, 'mixing_events.in'");
    }

    const int max_column_number = 87; /* bad to hard code this */
    float *dummy = new float[max_column_number];

    char line[MAX_LINE_LENGTH];

    int i = 0, d = 0;
    while ( fgets(line, MAX_LINE_LENGTH, fptr_mix) != NULL){
      if (line[0] != '#'){

        // just to be sure, reset dumyy variable every time
        for (d = 0; d < max_column_number; d++){
          dummy[d] = -1.0;
        }

        unpack_line_to_yields(line, dummy);

        MixingExperimentData.time[i] = dummy[0];     // time of event in code units

        // set time to negative if current time > event time (otherwise they will happen again on restars
        TimeActionTime[i]            = (time >= MixingExperimentData.time[i]) ?
                                       -MixingExperimentData.time[i] : MixingExperimentData.time[i];
        TimeActionType[i]            = 4;                            // hard coded for this experiment
        TimeActionParameter[i]       = 0.0;     // does not actually need to be set
        // TimeActionRedshift[i]     = 0.0;     // not currently used

        /* Positions are in code units */
        MixingExperimentData.xpos[i] = dummy[1];
        MixingExperimentData.ypos[i] = dummy[2];
        MixingExperimentData.zpos[i] = dummy[3];

        MixingExperimentData.M_ej[i] = dummy[4];    // Assumed to be in solar masses
        MixingExperimentData.E_ej[i] = dummy[5];    // Assumed to be in erg


        // for remaining fields, alternate between atomic number and
        // corresponding
        d = 6;
        while(dummy[d] > 0){

          int anum = dummy[d];

          // match index of stellar abundance with index of
          // species in global atomic numbers list just to keep things easier
          // when using current feedback injection machinery. All other
          // species followed but not listed in events table file will
          // have ejection masses of zero
          int index = -1;
          for (int j = 0; j < StellarYieldsNumberOfSpecies; j++){
            if (anum == StellarYieldsAtomicNumbers[j]) index = j;
          }

          // save all atomic numbers used for this experiment
          for (int count = 0; count < StellarYieldsNumberOfSpecies; count++){
            if (anum == MixingExperimentData.anums[count]){
              break;
            } else if (MixingExperimentData.anums[count] < 0){
              MixingExperimentData.anums[count] = anum;
            }
          }

          if (index < 0) ENZO_FAIL("Error initializing MetalMixingExperiment. Yield does not exist\n");

          MixingExperimentData.yield[i][index] = dummy[d+1]; // yield - assumed to be in Msun

          d += 2;
        }

        MixingExperimentData.NumberOfEvents++;
        i++;
      }
    }

    if (debug){
      fprintf(stdout,"Succesfully initialized Metal mixing experiment with %"ISYM" events\n", MixingExperimentData.NumberOfEvents);
      for(int i = 0; i < MixingExperimentData.NumberOfEvents; i ++){
        for(int j = 0; j < StellarYieldsNumberOfSpecies; j++){
          fprintf(stdout, " %"ESYM, MixingExperimentData.yield[i][j]);
        }
      fprintf(stdout,"\n");
      }
    }


    fclose(fptr_mix);
  }

  return SUCCESS;
}


void unpack_line_to_yields( char *line, float *dummy){
/* -----------------------------------------------------------
 * unpack_line_to_yields
 * -----------------------------------------------------------
 * Lets do something gross, :)
 * ----------------------------------------------------------- */

    int err;

    /* I'm sorry - really, I am*/
    err = sscanf(line, " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                   &dummy[ 0], &dummy[ 1], &dummy[ 2], &dummy[ 3], &dummy[ 4], &dummy[ 5], &dummy[ 6], &dummy[ 7], &dummy[ 8],
                   &dummy[ 9], &dummy[10], &dummy[11], &dummy[12], &dummy[13], &dummy[14], &dummy[15], &dummy[16],
                   &dummy[17], &dummy[18], &dummy[19], &dummy[20], &dummy[21], &dummy[22], &dummy[23], &dummy[24],
                   &dummy[25], &dummy[26], &dummy[27], &dummy[28], &dummy[29], &dummy[30], &dummy[31], &dummy[32],
                   &dummy[33], &dummy[34], &dummy[35], &dummy[36], &dummy[37], &dummy[38], &dummy[39], &dummy[40],
                   &dummy[41], &dummy[42], &dummy[43], &dummy[44], &dummy[45], &dummy[46], &dummy[47], &dummy[48],
                   &dummy[49], &dummy[50], &dummy[51], &dummy[52], &dummy[53], &dummy[54], &dummy[55], &dummy[56],
                   &dummy[57], &dummy[58], &dummy[59], &dummy[60], &dummy[61], &dummy[62], &dummy[63], &dummy[64],
                   &dummy[65], &dummy[66], &dummy[67], &dummy[68], &dummy[69], &dummy[70], &dummy[71], &dummy[72],
                   &dummy[73], &dummy[74], &dummy[75], &dummy[76], &dummy[77], &dummy[78], &dummy[79], &dummy[80],
                   &dummy[81], &dummy[82], &dummy[83], &dummy[84], &dummy[85], &dummy[86]);
}

void initialize_table(StellarYieldsDataType* table){
  /* -----------------------------------------------
   * fill_table
   * -----------------------------------------------
   */

  table->M    = new float[table->NumberOfMassBins];
  table->Z    = new float[table->NumberOfMetallicityBins];
  table->Mtot = new float*[table->NumberOfMassBins];
  table->Metal_Mtot = new float*[table->NumberOfMassBins];
  table->Yields = new float**[table->NumberOfMassBins];

  for (int i = 0; i < table->NumberOfMassBins; i++){
    table->Yields[i] = new float*[table->NumberOfMetallicityBins];

    table->Mtot[i] = new float[table->NumberOfMetallicityBins];
    table->Metal_Mtot[i] = new float[table->NumberOfMetallicityBins];

    for (int j = 0; j < table->NumberOfMetallicityBins; j++){
      table->Yields[i][j] = new float [table->NumberOfYields];
    }
  }

  return;
}

void fill_table(StellarYieldsDataType *table, FILE *fptr){

  const int max_column_number = 87;
  float *dummy = new float[max_column_number];

  char line[MAX_LINE_LENGTH];

  int i,j;
  i = 0; j = 0;
  while ( fgets(line, MAX_LINE_LENGTH, fptr) != NULL){
    if (line[0] != '#'){

      // just to be sure, reset dumyy variable every time
      for (int d = 0; d < max_column_number; d++){
        dummy[d] = 0.0;
      }

      unpack_line_to_yields(line, dummy);

      table->M[i] = dummy[0];
      table->Z[j] = dummy[1];
      table->Mtot[i][j]       = dummy[2];
      table->Metal_Mtot[i][j] = dummy[3];


      // file column numbers are atomic numbers + 1,
      // if first column is 0. Loop over number of yields
      // and pick only the ones we want
      for (int k = 0; k < table->NumberOfYields; k++){
        table->Yields[i][j][k] = dummy[3 + *(StellarYieldsAtomicNumbers+k)];
      }

      // iterate counters and reset if needed
      j++;
      if( j >= table->NumberOfMetallicityBins){
        j=0;
        i++;
      }
    } // end if
  }

  return;
}
