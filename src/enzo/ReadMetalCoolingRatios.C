/***********************************************************************
/
/  READ METAL-COOLING RATIOS
/
/  written by: John Wise
/  date:       May, 2008
/  modified1:
/
/ PURPOSE: For analysis purposes, read a table with ratios of a single
/          metal cooling line and total metal cooling tabulated by 
/          temperature and electron fraction.
/
/  RETURNS: ENZO_SUCCESS or FAIL
/
************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

#define NUMBER_OF_COOLANTS 11

int ReadMetalCoolingRatios(char *filename)
{

  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int i, nbins;

  if ((fptr = fopen(filename, "r")) == NULL) {
    ENZO_VFAIL("Error opening metal cooling table %s\n", filename)
  }

  // The second and third lines have the number of bins and temp/x_e
  // ranges, respectively.

  fgets(line, MAX_LINE_LENGTH, fptr);
  fgets(line, MAX_LINE_LENGTH, fptr);
  if ((sscanf(line, "# %"ISYM" %"ISYM, &CoolData.MR_NumberOfTemperatureBins, 
	      &CoolData.MR_NumberOfElectronFracBins)) != 2) {
    ENZO_FAIL("Error reading number of bins (line 2)\n");
  }

  fgets(line, MAX_LINE_LENGTH, fptr);
  if ((sscanf(line, "# %"GSYM" %"GSYM" %"GSYM" %"GSYM, 
	      &CoolData.MR_TemperatureStart, &CoolData.MR_TemperatureEnd,
	      &CoolData.MR_ElectronFracStart, &CoolData.MR_ElectronFracEnd)) != 4) {
    ENZO_FAIL("Error reading number of ranges (line 3)\n");
  }

  int prev_pos;
  int icool, ixe, itemp, index;
  float dummy;

  nbins = NUMBER_OF_COOLANTS * CoolData.MR_NumberOfTemperatureBins * 
    CoolData.MR_NumberOfElectronFracBins;
  CoolData.metal_ratios = new float[nbins];

  // Pass comments
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    if (line[0] != '#') {
      fseek(fptr, prev_pos, SEEK_SET);
      break;
    } // ENDIF not comment
    prev_pos = ftell(fptr);
  } // ENDWHILE file

    // Read the rate table
  index = 0;
  for (ixe = 0; ixe < CoolData.MR_NumberOfElectronFracBins; ixe++) {

    // Throw away line that indicates electron fraction bin and value
    fscanf(fptr, "%"GSYM, &dummy);
    fscanf(fptr, "%"GSYM, &dummy);
    //fgets(line, MAX_LINE_LENGTH, fptr);

    for (itemp = 0; itemp < CoolData.MR_NumberOfTemperatureBins; itemp++) {

      // Throw away first column (temperature)
      fscanf(fptr, "%"GSYM, &dummy);

      // Read ratios
      for (icool = 0; icool < NUMBER_OF_COOLANTS; icool++, index++) {
	if ((fscanf(fptr, "%"GSYM, &CoolData.metal_ratios[index])) == EOF) {
	  ENZO_VFAIL("EOF reached at itemp = %"ISYM", ixe = %"ISYM", icool = %"ISYM"\n", 
		  itemp, ixe, icool)

	}
      } // ENDFOR icool
    } // ENDFOR temperature
  } // ENDFOR electron fraction

  fclose(fptr);

  return SUCCESS;

}
