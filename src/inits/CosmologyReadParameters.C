/***********************************************************************
/
/  READS COSMOLOGY PARAMETERS FROM INPUT FILE
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:  Robert Harkness
/  date:       November, 2003
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
 
#include "macros_and_parameters.h"
#include "CosmologyParameters.h"
 
int CosmologyComputeTimeFromRedshift(FLOAT Redshift, FLOAT *TimeCodeUnits);
 
int CosmologyReadParameters(FILE *fptr)
{
 
  char line[MAX_LINE_LENGTH], *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
 
  // Set defaults
 
  HubbleConstantNow    = 0.701;
  OmegaMatterNow       = 0.279;
  OmegaLambdaNow       = 0.721;
  OmegaWDMNow          = 0.0;
  OmegaHDMNow          = 0;
  OmegaBaryonNow       = 0.0462;
  ComovingBoxSize      = 64;
  InitialRedshift      = 20;
 
  // Read input from file
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    int ret = 0;
 
    // Read parameters
 
    ret += sscanf(line, "CosmologyHubbleConstantNow = %"FSYM, &HubbleConstantNow);
    ret += sscanf(line, "CosmologyOmegaMatterNow = %"FSYM, &OmegaMatterNow);
    ret += sscanf(line, "CosmologyOmegaLambdaNow = %"FSYM, &OmegaLambdaNow);
    ret += sscanf(line, "CosmologyOmegaWDMNow = %"FSYM, &OmegaWDMNow);
    ret += sscanf(line, "CosmologyOmegaHDMNow = %"FSYM, &OmegaHDMNow);
    ret += sscanf(line, "CosmologyOmegaBaryonNow = %"FSYM, &OmegaBaryonNow);
    ret += sscanf(line, "CosmologyComovingBoxSize = %"FSYM, &ComovingBoxSize);
    ret += sscanf(line, "CosmologyInitialRedshift = %"FSYM, &InitialRedshift);
 
    // If the dummy char space was used, then make another
 
    if (*dummy != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
    // If the line is suspicious, issue a warning
 
    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#' && strstr(line, "Cosmology"))
      fprintf(stderr, "Warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }
 
/*
    fprintf(stderr, "CosmologyHubbleConstantNow = %16.8e\n", HubbleConstantNow);
    fprintf(stderr, "CosmologyOmegaMatterNow    = %16.8e\n", OmegaMatterNow);
    fprintf(stderr, "CosmologyOmegaLambdaNow    = %16.8e\n", OmegaLambdaNow);
    fprintf(stderr, "CosmologyOmegaWDMNow       = %16.8e\n", OmegaWDMNow);
    fprintf(stderr, "CosmologyOmegaHDMNow       = %16.8e\n", OmegaHDMNow);
    fprintf(stderr, "CosmologyOmegaBaryonNow    = %16.8e\n", OmegaBaryonNow);
    fprintf(stderr, "CosmologyComovingBoxSize   = %16.8e\n", ComovingBoxSize);
    fprintf(stderr, "CosmologyInitialRedshift   = %16.8e\n", InitialRedshift);
*/
 
  return SUCCESS;
}
