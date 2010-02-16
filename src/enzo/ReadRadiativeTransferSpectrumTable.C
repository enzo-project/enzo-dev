/***********************************************************************
/
/  READ THE SPECTRUM TABLE
/
/  written by: Ji-hoon Kim
/  date:       February, 2010
/  modified1:
/
/  PURPOSE:
/    For RadiaitveTransferTraceSpectrum = TRUE, initialize the spectrum
/    table.  The table must be tabulated against column density (rows),
/    and should contain (1) the fraction of photons left in the ray, and
/    (2) the mean energy of the spectrum at a given column density.
/
/  RETURNS: SUCCESS or FAIL
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

int ReadRadiativeTransferSpectrumTable(float TemperatureUnits, float LengthUnits, 
				       float aUnits, float DensityUnits, float TimeUnits)
{

  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int i, nbins;
  float TemperatureRange[2];

  // open the file

  if ((fptr = fopen(RadiativeTransferTraceSpectrumTable, "r")) == NULL) {
    fprintf(stderr, "Error opening spectrum table %s\n", 
	    RadiativeTransferTraceSpectrumTable);
    ENZO_FAIL("");
  }

  // The second line gives the number of lines in this file

  fgets(line, MAX_LINE_LENGTH, fptr); // pass the first line
  fgets(line, MAX_LINE_LENGTH, fptr);
  if ((sscanf(line, "# %"ISYM, &nbins)) != 1) {
    fprintf(stderr, "Error reading number of bins (line 2)\n");
    ENZO_FAIL("");
  }

  // Initialize

  RadiativeTransferSpectrumTable.NumberOfColumnDensityBins = nbins;
  RadiativeTransferSpectrumTable.columndensity_table = new float[nbins];
  RadiativeTransferSpectrumTable.fractionphotons_table[0] = new float[nbins];
  RadiativeTransferSpectrumTable.fractionphotons_table[1] = new float[nbins];
  RadiativeTransferSpectrumTable.fractionphotons_table[2] = new float[nbins];
  RadiativeTransferSpectrumTable.meanenergy_table = new float[nbins];

  // Read the spectrum table

  for (i = 0; i < nbins; i++) {
      if (fscanf(fptr, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %f",
		 &RadiativeTransferSpectrumTable.columndensity_table[i],
		 &RadiativeTransferSpectrumTable.fractionphotons_table[0][i],
		 &RadiativeTransferSpectrumTable.fractionphotons_table[1][i],
		 &RadiativeTransferSpectrumTable.fractionphotons_table[2][i],
		 &RadiativeTransferSpectrumTable.meanenergy_table[i])
	  != 5) {
	fprintf(stderr, "Error reading RadiationData line %"ISYM"\n", i);
	ENZO_FAIL("");
      }
  }

  fclose(fptr);


#ifdef OUTPUT_SPECTRUM_FOR_CHECK
  // Write out the spectrum table for check

  FILE *fptr2;

  fptr2 = fopen("spectrum_table.out", "w");
  for (i = 0; i < RadiativeTransferSpectrumTable.NumberOfColumnDensityBins; i++) 
    fprintf(fptr2, "%g  %"FSYM"  %"FSYM"  %"FSYM"  %f\n", 
	    RadiativeTransferSpectrumTable.columndensity_table[i],
	    RadiativeTransferSpectrumTable.fractionphotons_table[0][i],
	    RadiativeTransferSpectrumTable.fractionphotons_table[1][i],
	    RadiativeTransferSpectrumTable.fractionphotons_table[2][i],
	    RadiativeTransferSpectrumTable.meanenergy_table[i]);

  fclose(fptr2);
#endif

  return SUCCESS;

}
