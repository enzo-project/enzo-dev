/***********************************************************************
/
/  INITIALIZE THE TABULATED HM12 Photo-ionization/heating rates
/
/  written by: Gabriel Altay
/  date:       April, 2013
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
  
int InitializeHM12Photorates()
{


  /* Open input file for data. */
 
  FILE *fptr = fopen("hm12_photorates.dat", "r");
  if (fptr == NULL) {
    ENZO_FAIL("Error opening hm12_photorates.dat\n");
  }
 
  /* Read rate data, skipping over comments (count first). 
   Note that photoheating rates are per ion in eV/s */
 
  int index = 0;
  char line[MAX_LINE_LENGTH];
  double ev2erg;

  ev2erg = 1.60217657e-12;

  RateData.HM12NumberOfRedshiftBins = 0;

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    if (line[0] != '#') RateData.HM12NumberOfRedshiftBins++;
  }

  RateData.HM12Redshifts = new float[RateData.HM12NumberOfRedshiftBins];
  RateData.HM12GH1 = new float[RateData.HM12NumberOfRedshiftBins];
  RateData.HM12GHe1 = new float[RateData.HM12NumberOfRedshiftBins];
  RateData.HM12GHe2 = new float[RateData.HM12NumberOfRedshiftBins];
  RateData.HM12GhH1 = new float[RateData.HM12NumberOfRedshiftBins];
  RateData.HM12GhHe1 = new float[RateData.HM12NumberOfRedshiftBins];
  RateData.HM12GhHe2 = new float[RateData.HM12NumberOfRedshiftBins];  
  RateData.HM12Compton = new float[RateData.HM12NumberOfRedshiftBins];  

  rewind(fptr);
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    if (line[0] != '#')
      if (sscanf(line, "%"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM,
		 &RateData.HM12Redshifts[index],
		 &RateData.HM12GH1[index],
		 &RateData.HM12GhH1[index],
		 &RateData.HM12GHe1[index],
		 &RateData.HM12GhHe1[index],
		 &RateData.HM12GHe2[index],
		 &RateData.HM12GhHe2[index], 
		 &RateData.HM12Compton[index]) == 8) {
	
	RateData.HM12GH1[index] = log10( RateData.HM12GH1[index] ); 
	RateData.HM12GhH1[index] = log10( RateData.HM12GhH1[index]*ev2erg ); 

	RateData.HM12GHe1[index] = log10( RateData.HM12GHe1[index] );
	RateData.HM12GhHe1[index] = log10( RateData.HM12GhHe1[index]*ev2erg ); 

	RateData.HM12GHe2[index]  = log10( RateData.HM12GHe2[index] ); 
	RateData.HM12GhHe2[index] = log10( RateData.HM12GhHe2[index]*ev2erg );  

	RateData.HM12Compton[index] = log10( RateData.HM12Compton[index]*ev2erg ); 
	
	index++;
	
      }
  }

  fclose(fptr);
 
  RateData.HM12RedshiftLo = RateData.HM12Redshifts[0];
  RateData.HM12RedshiftHi = RateData.HM12Redshifts[RateData.HM12NumberOfRedshiftBins-1];

  if (debug) {
    printf("InitializeHM12Photorates: HM12NumberOfRedshiftBins = %"ISYM"\n",
	   RateData.HM12NumberOfRedshiftBins);
    printf("InitializeHM12Photorates: HM12RedshiftLo = %"ESYM"\n",
	   RateData.HM12RedshiftLo);
    printf("InitializeHM12Photorates: HM12RedshiftHi = %"ESYM"\n",
	   RateData.HM12RedshiftHi);
  }

 
  return SUCCESS;
}
