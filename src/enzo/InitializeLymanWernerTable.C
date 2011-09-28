/***********************************************************************
/
/  INITIALIZE THE TABULATED LW BACKGROUND
/
/  written by: Michele Trenti and Britton Smith
/  date:       April, 2009
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
  
int InitializeLymanWernerTable()
{
  
  /* Open input file for data. */
 
  FILE *fptr = fopen("LW_J21.in", "r");
  if (fptr == NULL) {
    ENZO_FAIL("Error opening LW_J21.in\n");
  }
 
  /* Read rate data, skipping over comments (count first). */
 
  int index = 0;
  char line[MAX_LINE_LENGTH];
  RadiationData.NumberOfLWRedshiftBins = 0;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    if (line[0] != '#') RadiationData.NumberOfLWRedshiftBins++;
  RadiationData.LymanWerner_J21 = new float[RadiationData.NumberOfLWRedshiftBins];
  RadiationData.LymanWerner_redshift = new float[RadiationData.NumberOfLWRedshiftBins];
  rewind(fptr);
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    if (line[0] != '#')
      if (sscanf(line, "%"ESYM" %"ESYM,
		 &RadiationData.LymanWerner_redshift[index],
		 &RadiationData.LymanWerner_J21[index]) == 2) {
	index++;
      }
  RadiationData.NumberOfLWRedshiftBins = index;
  fclose(fptr);
 
  if (debug) {
    printf("InitializeLymanWernerTable: NumberOfLWRedshiftBins = %"ISYM"\n",
	   RadiationData.NumberOfLWRedshiftBins);
    printf("InitializeLymanWernerTable: RedshiftStart = %"ESYM"\n",
	   RadiationData.LymanWerner_redshift[0]);
    printf("InitializeLymanWernerTable: RedshiftEnd = %"ESYM"\n",
	   RadiationData.LymanWerner_redshift[RadiationData.NumberOfLWRedshiftBins-1]);
  }
 
  return SUCCESS;
}
