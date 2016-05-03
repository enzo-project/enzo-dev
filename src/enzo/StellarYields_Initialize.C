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
#include "StarParticleData.h"
#include "StellarYieldsRoutines.h"

/* function prototypes */
void unpack_line_to_yields( char *line, float *dummy);



int InitializeStellarYields(void){
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
  if( !IndividualStarFollowStellarYields &&
      !TestProblemData.MultiMetals       &&
      !STARMAKE_METHOD(INDIVIDUAL_STAR)) {
    printf("Returning.... did not meet initialize conditions");
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

  // read in data from files - one table for each yield type:
  //   1) core collapse supernova
  //   2) stellar winds
  //
  FILE *fptr_sn = fopen("stellar_yields_sn.in", "r");
  if (fptr_sn == NULL){
    ENZO_FAIL("Error opening stellar yields SN file, 'stellar_yields_sn.in");
  }
  FILE *fptr_wind = fopen("stellar_yields_wind.in", "r");
  if (fptr_wind == NULL){
    ENZO_FAIL("Error opening stellar yields wind file, 'stellar_yields_wind.in'");
  }

  char line[MAX_LINE_LENGTH];

  /* bins are not evenly spaced - save bin values */
  StellarYieldsSNData.M = new float[StellarYieldsSNData.NumberOfMassBins];
  StellarYieldsSNData.Z = new float[StellarYieldsSNData.NumberOfMetallicityBins];

  StellarYieldsWindData.M = new float[StellarYieldsWindData.NumberOfMassBins];
  StellarYieldsWindData.Z = new float[StellarYieldsWindData.NumberOfMetallicityBins];

  /* total mass for each M and Z */
  StellarYieldsSNData.Mtot   = new float*[StellarYieldsSNData.NumberOfMassBins];
  StellarYieldsSNData.Metal_Mtot = new float*[StellarYieldsSNData.NumberOfMassBins];

  StellarYieldsWindData.Mtot = new float*[StellarYieldsWindData.NumberOfMassBins];
  StellarYieldsWindData.Metal_Mtot = new float*[StellarYieldsWindData.NumberOfMassBins];

  /* start making arrays for the data */
  StellarYieldsSNData.Yields   = new float**[StellarYieldsSNData.NumberOfMassBins];
  StellarYieldsWindData.Yields = new float**[StellarYieldsWindData.NumberOfMassBins];

  for (int i = 0; i < StellarYieldsSNData.NumberOfMassBins; i++){
    StellarYieldsSNData.Yields[i] = new float*[StellarYieldsSNData.NumberOfMetallicityBins];

    StellarYieldsSNData.Mtot[i]   = new float[StellarYieldsSNData.NumberOfMetallicityBins];
    StellarYieldsSNData.Metal_Mtot[i] = new float[StellarYieldsSNData.NumberOfMetallicityBins];

    for (int j = 0; j < StellarYieldsSNData.NumberOfMetallicityBins; j++){ 
      StellarYieldsSNData.Yields[i][j] = new float[ StellarYieldsSNData.NumberOfYields ]; /* AJE: Make parameter */
    }
  }

  for (int i = 0; i < StellarYieldsWindData.NumberOfMassBins; i++){
    StellarYieldsWindData.Yields[i] = new float*[StellarYieldsWindData.NumberOfMetallicityBins];

    StellarYieldsWindData.Mtot[i]   = new float[StellarYieldsWindData.NumberOfMetallicityBins];
    StellarYieldsWindData.Metal_Mtot[i] = new float[StellarYieldsWindData.NumberOfMetallicityBins];

    for (int j = 0; j < StellarYieldsWindData.NumberOfMetallicityBins; j++){ 
      StellarYieldsWindData.Yields[i][j] = new float[ StellarYieldsWindData.NumberOfYields ]; /* AJE: Make parameter */
    }
  }

  /* now read in the data from file */
  int i, j, err;
  i = 0; j = 0;

  float *dummy = new float[87];

  /* We are doing the unnatractive thing of unpacking a 80+ column file when
     we only want a few columns. Hence the function call to unpack the line 
     to a dummy variable */
  while( fgets(line, MAX_LINE_LENGTH, fptr_sn) != NULL){
    if(line[0] != '#'){

      unpack_line_to_yields(line, dummy);

      StellarYieldsSNData.M[i]                = dummy[0];
      StellarYieldsSNData.Z[j]                = dummy[1];
      StellarYieldsSNData.Mtot[i][j]          = dummy[2];
      StellarYieldsSNData.Metal_Mtot[i][j]    = dummy[3];

      /* file column numbers are atomic numbers + 1, if first column is 0
         loop over number of yields provided by user and pick only the ones we
         want  2 in below is number of non-element columns - 1 */
      for (int k = 0; k < StellarYieldsSNData.NumberOfYields; k++){
        StellarYieldsSNData.Yields[i][j][k] = dummy[3 + *(StellarYieldsAtomicNumbers+k)];
      }

      j++;
      if ( j >= StellarYieldsSNData.NumberOfMetallicityBins){
        j = 0;
        i++;
      }
    }
  } // end sn yields read in

  i = 0; j = 0;

  while( fgets(line, MAX_LINE_LENGTH, fptr_wind) != NULL){
    if(line[0] != '#'){

      unpack_line_to_yields(line, dummy);

      StellarYieldsWindData.M[i]                = dummy[0];
      StellarYieldsWindData.Z[j]                = dummy[1];
      StellarYieldsWindData.Mtot[i][j]          = dummy[2];
      StellarYieldsWindData.Metal_Mtot[i][j]    = dummy[3];

      /* see above comment in previous while loop */
      for (int k = 0; k < StellarYieldsWindData.NumberOfYields; k++){
        StellarYieldsWindData.Yields[i][j][k] = dummy[3 + *(StellarYieldsAtomicNumbers+k)];
      }

      j++;
      if ( j >= StellarYieldsWindData.NumberOfMetallicityBins){
        j = 0;
        i++;
      }
    }
  }

  delete [] dummy;


  /* close files */
  fclose(fptr_sn);
  fclose(fptr_wind);

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

