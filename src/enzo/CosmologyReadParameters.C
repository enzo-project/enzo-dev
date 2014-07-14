/***********************************************************************
/
/  READS COSMOLOGY PARAMETERS FROM INPUT FILE
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
 
int CosmologyComputeTimeFromRedshift(FLOAT Redshift, FLOAT *TimeCodeUnits);
 
int CosmologyReadParameters(FILE *fptr, FLOAT *StopTime, FLOAT *InitTime)
{
 
  int i, OutputNumber;
  FLOAT CurrentRedshift;
  char line[MAX_LINE_LENGTH], *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
 
  /* Set defaults. */
 
  HubbleConstantNow    = 0.701;
  OmegaMatterNow       = 0.279;
  OmegaLambdaNow       = 0.721;
  ComovingBoxSize      = 64;
  MaxExpansionRate     = 0.01;
  InitialRedshift      = 20;
  FinalRedshift        = 0;
 
  for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++) {
    CosmologyOutputRedshift[i]     = -1;  // Never!!
    CosmologyOutputRedshiftName[i] = NULL;
  }
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    int ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "CosmologyHubbleConstantNow = %"FSYM,
		  &HubbleConstantNow);
    ret += sscanf(line, "CosmologyOmegaMatterNow = %"FSYM, &OmegaMatterNow);
    ret += sscanf(line, "CosmologyOmegaLambdaNow = %"FSYM, &OmegaLambdaNow);
    ret += sscanf(line, "CosmologyComovingBoxSize = %"FSYM, &ComovingBoxSize);
    ret += sscanf(line, "CosmologyMaxExpansionRate = %"FSYM,
		  &MaxExpansionRate);
    ret += sscanf(line, "CosmologyInitialRedshift = %"PSYM, &InitialRedshift);
    ret += sscanf(line, "CosmologyFinalRedshift = %"PSYM, &FinalRedshift);
    ret += sscanf(line, "CosmologyCurrentRedshift = %"PSYM, &CurrentRedshift);
 
    if (sscanf(line, "CosmologyOutputRedshift[%"ISYM"] =", &OutputNumber) == 1)
      ret += sscanf(line, "CosmologyOutputRedshift[%"ISYM"] = %"PSYM,
		    &OutputNumber, &CosmologyOutputRedshift[OutputNumber]);
    if (sscanf(line, "CosmologyOutputRedshiftName[%"ISYM"] = %s",
	       &OutputNumber, dummy) == 2)
      CosmologyOutputRedshiftName[OutputNumber] = dummy;
 
    /* If the dummy char space was used, then make another. */
 
    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#' &&
	strstr(line, "Cosmology") && !strstr(line, "CosmologySimulation"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
  /* Initialize by finding the time at the initial redshift. */
 
  if (CosmologyComputeTimeFromRedshift(InitialRedshift,
				       &InitialTimeInCodeUnits) == FAIL) {
    ENZO_FAIL("Error in ComputeTimeFromRedshift.\n");
  }
  if (*InitTime == 0.0)
    *InitTime = InitialTimeInCodeUnits;
 
  /* Now find the time at the end of the simulation. */
 
  if (CosmologyComputeTimeFromRedshift(FinalRedshift, StopTime) == FAIL) {
    ENZO_FAIL("Error in ComputeTimeFromRedshift.\n");
  }
 
  /* Convert the output redshifts into time, for later convenience. */
 
  for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
    if (CosmologyOutputRedshift[i] != -1)
      CosmologyComputeTimeFromRedshift(CosmologyOutputRedshift[i],
				       &CosmologyOutputRedshiftTime[i]);
 
  /* Convert the time action redshift into time. */
 
  for (i = 0; i < MAX_TIME_ACTIONS; i++)
    if (TimeActionRedshift[i] != -1)

      CosmologyComputeTimeFromRedshift(TimeActionRedshift[i],
				       &TimeActionTime[i]);

  delete [] dummy;
 
  return SUCCESS;
}
