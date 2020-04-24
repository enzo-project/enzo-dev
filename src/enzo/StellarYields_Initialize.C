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
//#include <hdf5.h>
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

#ifdef NEWYIELDTABLES
int fill_table(hid_t file_id,
               StellarYieldsDataType *table,
               std::string filename,
               std::string dname,
               int null_if_no_group = FALSE);
#else
void fill_table(StellarYieldsDataType *table, FILE *fptr);
#endif

int ChemicalSpeciesBaryonFieldNumber(const int &atomic_number);
char* ChemicalSpeciesBaryonFieldLabelByFieldType(const int &field_num);


int InitializeStellarYieldFields(HierarchyEntry &TopGrid,
                                 TopGridData &MetaData,
                                 ExternalBoundary &Exterior,
                                 LevelHierarchyEntry *LevelArray[]){
// Initializes species yields if they do not already exist


  if ( !IndividualStarFollowStellarYields ||
       !(MultiMetals)       ||
       !STARMAKE_METHOD(INDIVIDUAL_STAR)){
    return SUCCESS;
  }

  int OldNumberOfBaryonFields = 0, FieldsToAdd = 0;
  int TypesToAdd[MAX_NUMBER_OF_BARYON_FIELDS];
  int ExistingTypes[MAX_NUMBER_OF_BARYON_FIELDS];

  for (int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i ++)
    ExistingTypes[i] = FieldUndefined;

  for(int yield_i = 0; yield_i < StellarYieldsNumberOfSpecies; yield_i++){
    if(StellarYieldsAtomicNumbers[yield_i] > 2){
      TypesToAdd[FieldsToAdd++] =
                     ChemicalSpeciesBaryonFieldNumber(StellarYieldsAtomicNumbers[yield_i]);
    }
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
  if (FieldsToAdd > 0 && debug){
    fprintf(stdout, "InitializeStellarYieldsFields: Increasing baryon fields "
             "from %"ISYM" to %"ISYM"\n", OldNumberOfBaryonFields,
              OldNumberOfBaryonFields + FieldsToAdd);
    fprintf(stdout, "Adding:   \n");
    for (int k = 0; k < FieldsToAdd; k ++){
      fprintf(stdout, "Field Number  %"ISYM"\n", TypesToAdd[k]);
    }
  }

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
//   if(StellarYieldsAtomicNumbers[i] > 2){
     DataLabel[OldNumberOfBaryonFields+i] =\
          ChemicalSpeciesBaryonFieldLabelByFieldType(TypesToAdd[i]);
//   }
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
      !MultiMetals       ||
      !STARMAKE_METHOD(INDIVIDUAL_STAR)) {
    return SUCCESS;

  } else if (IndividualStarFollowStellarYields && !MultiMetals){
    printf("Failure in InitializeStellarYields. MultiMetals must be enabled to follow yields\n");
    return FAIL;
  }

  if (StellarYieldsSNData.M != NULL && StellarYieldsSNData.Z != NULL){
    return SUCCESS; // already initialized
  }

#ifdef NEWYIELDTABLES

  std::string filename = "IndividualStarYields.h5";

  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  fill_table(file_id, &StellarYieldsSNData,           filename, "SN");
  fill_table(file_id, &StellarYieldsWindData,         filename, "Wind");
  fill_table(file_id, &StellarYieldsPopIIIData,       filename, "PopIII");

// treat these two a little differently for now so this is backwards compatabile
// with original model.
  int agb_check     = fill_table(file_id, &StellarYieldsAGBData,          filename, "AGB", TRUE);
  int massive_check = fill_table(file_id, &StellarYieldsMassiveStarData,  filename, "Massive_star", TRUE);

  if (agb_check && massive_check){
    ENZO_FAIL("Yield table has both AGB and Massive Star data. Do not know how to handle this.");
  } else if ( !(agb_check) && !(massive_check)){
    ENZO_FAIL("Yield table does not have either AGB or Massive_star data. Do not know how to handle this.");
  }

  // if agb_check is true, AGB yields are present and we are using the new methods
  // if massive_check is true, massive star stellar wind yields are available and we
  //     are using the old methods. In this case, the Wind yields contain both
  //     AGB stars and stars up to 25 Msun. Massive_star contains winds for stars above 25 Msun.

  herr_t status = H5Fclose (file_id);
  herr_t h5_error = -1;
  if (status == h5_error){
    ENZO_VFAIL("Error closing %s \n", filename.c_str());
  }

#else


  // AJE: Hard code he number of bins for now
  //     - I want to fix this but this is not priority -
  //     - unfixed as of April 2016 -
  StellarYieldsSNData.Nm        = 12;
  StellarYieldsSNData.Nz =  5;
  StellarYieldsSNData.Ny          = StellarYieldsNumberOfSpecies;

  StellarYieldsWindData.Nm        = 12;
  StellarYieldsWindData.Nz =  5;
  StellarYieldsWindData.Ny          = StellarYieldsNumberOfSpecies;

  StellarYieldsMassiveStarData.Nm        = 30;
  StellarYieldsMassiveStarData.Nz = 12;
  StellarYieldsMassiveStarData.Ny       = StellarYieldsNumberOfSpecies;

  StellarYieldsPopIIIData.Nm = 120 + 14 ;  // 120 from Heger+Woosley2010 for Type II (10 < M < 100)
                                                         //  14 from Heger+Woosley2002 for PISN    (140 < M < 260)
  StellarYieldsPopIIIData.Nz = 1;
  StellarYieldsPopIIIData.Ny = StellarYieldsNumberOfSpecies;

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


#endif


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


#ifdef NEWYIELDTABLES
void initialize_table(StellarYieldsDataType* table){

  /* fill table in 1D - makes lookup faster*/

  const int Nm = table->Nm;
  const int Nz = table->Nz;
  const int Ny = table->Ny;

  table->M    = new float[Nm];
  table->Z    = new float[Nz];
  table->Mtot = new float[Nm*Nz];
  table->Metal_Mtot = new float [Nm*Nz];
  table->Yields = new float [Nm*Nz*Ny];

  for (int i = 0; i < Nm; i++){
    table->M[i] = 0.0;

    for (int j = 0; j < Nz; j++){
      table->Z[j] = 0.0;

      table->Mtot[i + j*Nm] = 0.0;
      table->Metal_Mtot[i + j*Nm] = 0.0;

      for (int k =0; k < Ny; k++){
        table->Yields[i + (j + k*Nz)*Nm] = 0.0;
      }
    }
  }

  return;
}

#else

void initialize_table(StellarYieldsDataType* table){
  /* -----------------------------------------------
   * Initialize table
   * -----------------------------------------------
   */

  table->M    = new float[table->Nm];
  table->Z    = new float[table->Nz];
  table->Mtot = new float*[table->Nm];
  table->Metal_Mtot = new float*[table->Nm];
  table->Yields = new float**[table->Nm];

  for (int i = 0; i < table->Nm; i++){
    table->Yields[i] = new float*[table->Nz];

    table->Mtot[i] = new float[table->Nz];
    table->Metal_Mtot[i] = new float[table->Nz];

    for (int j = 0; j < table->Nz; j++){
      table->Yields[i][j] = new float [table->Ny];
    }
  }

  return;
}

#endif


#ifdef NEWYIELDTABLES

int fill_table(hid_t file_id,
               StellarYieldsDataType *table,
               std::string filename,
               std::string dname,
               int null_if_no_group /*Default FALSE */
               ){

  /*
    Generic routine to read in stellar yields from an
    HDF5 file. Yields are assumed to be the mass yield (in Msun)
    for grid points in stellar mass (Msun) and metallicity (fraction),
    requiring 2D interpolation.

    If null_if_no_group is TRUE, checks to see if the group exists first
    and fills the struct pointers with NULL if it does not, returning FAIL.

    If FALSE (default) calls ENZO_FAIL if it cannot find group.

  */


  /* Read an HDF5 table */

  hid_t    dset_id, dspace_id;
  herr_t   status;
  herr_t   h5_error = -1;


  // Read Info dataset
  //  (not necessary)
/*
  dset_id = H5Dopen(file_id, "/Info");
  if (dset_id == h5error){
    ENZO_VFAIL("Can't open 'Info' in dataset in %s.\n",
               filename.c_str());
  }

  int strlen = (int)(H5Dget_storage_size(dset_id));
  char
*/

  // Set yields to load (don't load all available elements!)

  table->Ny = StellarYieldsNumberOfSpecies;

  // Check if desired group exists

  status = H5Gget_objinfo(file_id, ("/"+dname).c_str(), 0, NULL);
  if (status == h5_error){

    if (null_if_no_group){
      // Allow code to keep working, but set table pointers to NULL
      printf("Group name %s not found in %s. Continuing anyway\n",dname.c_str(),filename.c_str());

      table->Nm = -1; table->Nz = -1; table->Ny = -1;

      table->dm = 0; table->dy = 0; table->dz = 0;

      table->M  = NULL;
      table->Z  = NULL;
      table->Mtot = NULL;
      table->Metal_Mtot = NULL;
      table->Yields = NULL;

      return FAIL;

    } else {
      ENZO_VFAIL("No group %s found in %s\n", dname.c_str(), filename.c_str());
    }
  }

  if (null_if_no_group){

  }
  // Find mass bins in desired dataset

  dset_id = H5Dopen(file_id, ("/"+dname+"/M").c_str());
  if (dset_id == h5_error){
    ENZO_VFAIL("Can't open 'M' in %s in file %s.\n",dname.c_str(),
                                                    filename.c_str());
  }

  dspace_id = H5Dget_space(dset_id);
  if (dspace_id == h5_error){
    ENZO_VFAIL("Can't open 'M' dataspace in %s in file %s.\n",
               dname.c_str(), filename.c_str());
  }

  table->Nm = H5Sget_simple_extent_npoints(dspace_id);
  if (table->Nm <= 0) {
    ENZO_VFAIL("Cannot propertly read mass bins ('M') in %s in %s.\n",
               dname.c_str(), filename.c_str());
  }

  H5Sclose(dspace_id);
  H5Dclose(dset_id);

  // Find metallicity bins in desired dataset

  dset_id = H5Dopen(file_id, ("/"+dname+"/Z").c_str());
  if (dset_id == h5_error){
    ENZO_VFAIL("Can't open 'Z' in %s in file %s.\n",dname.c_str(),
                                                    filename.c_str());
  }

  dspace_id = H5Dget_space(dset_id);
  if (dspace_id == h5_error){
    ENZO_VFAIL("Can't open 'Z' dataspace in %s in file %s.\n",
               dname.c_str(), filename.c_str());
  }

  table->Nz = H5Sget_simple_extent_npoints(dspace_id);
  if (table->Nz <= 0) {
    ENZO_VFAIL("Cannot propertly read Z bins ('Z') in %s in %s.\n",
               dname.c_str(), filename.c_str());
  }

  H5Sclose(dspace_id);
  H5Dclose(dset_id);

  // Find total number of available yields in the dataset

  dset_id = H5Dopen(file_id, ("/"+dname+"/atomic_numbers").c_str());
  if (dset_id == h5_error){
    ENZO_VFAIL("Can't open 'atomic_numbers' in %s in file %s.\n",dname.c_str(),
                                                    filename.c_str());
  }

  dspace_id = H5Dget_space(dset_id);
  if (dspace_id == h5_error){
    ENZO_VFAIL("Can't open 'atomic_numbers' dataspace in %s in file %s.\n",
               dname.c_str(), filename.c_str());
  }

  int Nyields = H5Sget_simple_extent_npoints(dspace_id);
  if (Nyields <= 0) {
    ENZO_VFAIL("Cannot propertly read number of yields from 'atomic_numbers' in %s in %s.\n",
               dname.c_str(), filename.c_str());
  }

  H5Sclose(dspace_id);
  H5Dclose(dset_id);

  table->size = table->Nm * table->Nz * table->Ny;

  // allocate space for table

  initialize_table(table);

  // Now read in yields

  if(! read_dataset(file_id, ("/"+dname+"/M").c_str(),
                    table->M)){
    ENZO_VFAIL("Error reading dataset 'M' in %s in %s.\n",
               dname.c_str(), filename.c_str());
  }

  if( !read_dataset(file_id, ("/"+dname+"/Z").c_str(),
                    table->Z)){
    ENZO_VFAIL("Error reading dataset 'Z' in %s in %s.\n",
               dname.c_str(), filename.c_str());
  }

  // Find list of atomic numbers in the dataset
  //   list *should* start with -1 and 0 since these
  //   are the codes for total mass (-1) and metal mass (0)
  //   which should be first two columns in yields

  int * temp_anum = new int [Nyields];
  for (int i = 0; i < Nyields; i++) temp_anum[i] = -1;

  dset_id = H5Dopen(file_id, ("/"+dname+"/atomic_numbers").c_str());
  if (dset_id == h5_error){
    ENZO_VFAIL("Error opening atomic_numbers for %s in %s\n",
               dname.c_str(), filename.c_str());
  }

  status = H5Dread(dset_id, HDF5_I8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_anum);
  if (status == h5_error){
    ENZO_VFAIL("Error reading in atomic_numbers for %s in %s\n",
               dname.c_str(), filename.c_str());
  }

  // now, allocate a temporary array to load in ALL yields
  // from the table

  int temp_size = table->Nm * table->Nz *
                   Nyields;

  float *temp_yields = new float [temp_size];
  for(int i = 0 ; i < temp_size; i++) temp_yields[i] = 0.0;

  dset_id = H5Dopen(file_id, ("/"+dname+"/yields").c_str());
  if (dset_id == h5_error){
    ENZO_VFAIL("Error opening yields for %s in %s\n",
               dname.c_str(), filename.c_str());
  }

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_yields);
  if (status == h5_error){
    ENZO_VFAIL("Error reading yields dataset for %s in %s\n",
               dname.c_str(), filename.c_str());
  }

  // first column of yields is the total mass ejected
  // second column is the total metals
  /* We want the table to be constructed such that adjacent
     values are adjacent in mass at fixed Z and element. They
     are read-in in reversed order, however */

  for(int j =0; j < table->Nz; j++){
    for(int i =0; i < table->Nm; i++){

      int index = YIELD_INDEX(i,j,0,table->Nm,table->Nz);

      table->Mtot[index] =
         temp_yields[0 + (j + i*table->Nz)*Nyields];

      table->Metal_Mtot[index] =
         temp_yields[1 + (j + i*table->Nz)*Nyields];
    }
  }
  // now loop through and grab the elements we need
  int temp_k = 0;
  for (int k = 0; k < table->Ny; k++){

    // atomic number we want
    int anum = StellarYieldsAtomicNumbers[k];

    // corresponding yield column in array
    // assume yield table is in atomic number order
    int kk;
    for(kk = temp_k; kk < Nyields; kk++){
      if (anum == temp_anum[kk]){
        temp_k = kk;
        break;
      }
    }

    if (kk >= Nyields){

      // this means that the yield is not present
      // print a warning and set value to zero below

      printf("Yield set %s in yield table %s does not have anum = %i. Setting to 0 \n",
             dname.c_str(), filename.c_str(), anum);

      for (int j =0; j < table->Nz; j++){
        for(int i = 0; i < table->Nm; i++){
          table->Yields[YIELD_INDEX(i,j,k,table->Nm,table->Nz)] = 0.0;
        }
      }

      temp_k = 0; // reset

    } else {

      for (int j = 0; j < table->Nz; j++){
        for(int i = 0; i < table->Nm; i++){

          // get index in saved table yield set
          int index = YIELD_INDEX(i,j,k,table->Nm,table->Nz);

          // get index in temporary yield set
          int itemp  = (temp_k) + (j + i*table->Nz)*Nyields;

          table->Yields[index] = temp_yields[itemp];
        }
      }
    } // else if available
  }

  /* Save index offsets for next item in each dimension for convenience */
  table->dm = 1;                     // next mass
  table->dz = table->Nm;             // next metallicity
  table->dy = table->Nm * table->Nz; // next yield

  status = H5Dclose(dset_id);
  if (status == h5_error){
    ENZO_VFAIL("Error closing yields dataset in %s in %s\n",
               dname.c_str(), filename.c_str());
  }

  /* Delete temporary yields and atomic numbers */
  delete [] temp_yields;
  delete [] temp_anum;

  return SUCCESS;
}

#else

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
      for (int k = 0; k < table->Ny; k++){
        table->Yields[i][j][k] = dummy[3 + *(StellarYieldsAtomicNumbers+k)];
      }

      // iterate counters and reset if needed
      j++;
      if( j >= table->Nz){
        j=0;
        i++;
      }
    } // end if
  }

  return;
}

#endif

#ifdef NEWYIELDTABLES
int read_dataset(hid_t file_id, const char *dset_name, double *buffer) {
  hid_t dset_id;
  herr_t status;
  herr_t h5_error = -1;

  dset_id =  H5Dopen(file_id, dset_name);
  if (dset_id == h5_error) {
    fprintf(stderr, "Failed to open dataset 'z'.\n");
    return FAIL;
  }

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  if (status == h5_error) {
    fprintf(stderr, "Failed to read dataset 'z'.\n");
    return FAIL;
  }

  H5Dclose(dset_id);

  return SUCCESS;
}
#endif
