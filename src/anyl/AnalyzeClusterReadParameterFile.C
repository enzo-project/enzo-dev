/***********************************************************************
/
/  READS A PARAMETER FILE FOR ANALYZE CLUSTER
/
/  written by: Greg Bryan
/  date:       August, 1997
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../enzo/macros_and_parameters.h"
#include "../enzo/typedefs.h"
#include "../enzo/global_data.h"
#include "../enzo/CosmologyParameters.h"
#include "AnalyzeClusters.h"



/* function prototypes */

int AnalyzeClusterReadParameterFile(char *filename, int &NumberOfCenters,
				    FLOAT *CenterList[],
				    AnalyzeClusterParameters *parm)

{

  int dim, ret, j;
  char line[MAX_LINE_LENGTH], *char_dummy = new char[MAX_LINE_LENGTH];
  FLOAT center[MAX_DIMENSION], float_dummy;

  /* Set default vaules. */

  if (ComovingCoordinates != 1) {  
    InitialRedshift = 0; 
    //FinalRedshift = 0;
    HubbleConstantNow = 0.7; 
    OmegaMatterNow = 0.3;
    OmegaLambdaNow = 0.7;
    //    float ComovingBoxSize = 1;
    //    float MaxExpansionRate = 1;
  }  

  float BoxSize = 1;

  if (ComovingCoordinates)
    BoxSize = ComovingBoxSize/HubbleConstantNow;

  char *CenterListName = NULL;
  center[0]   = FLOAT_UNDEFINED;

  parm->rinner      = 0.0001*BoxSize; 
  parm->router      = 0.1*BoxSize;
  parm->npoints     = 16;
  parm->virial_dens = FLOAT_UNDEFINED;
  parm->MeanVelocityVirialFraction = 1.0;
  parm->ColdTemperatureCutoff = 15000;    // in K
  parm->ColdTemperatureCutoffVirialFraction = FLOAT_UNDEFINED;
  parm->VirialTemperatureNormalization = 1.0;  // M-T-z normalization
  parm->LowerDensityCutoff     = 1.0e14;  /* in solar masses/Mpc^3 */
  parm->UpperDensityCutoff     = 1.0e35;  /* a big number */
  parm->ComputeDiskInformation = FALSE;
  parm->DiskImageSize          = 100;
  parm->DiskRadius             = 0.2;  // as a fraction of virial radius
  parm->XrayLowerCutoffkeV     = 0.5;
  parm->XrayUpperCutoffkeV     = 2.5;
  parm->XrayTableFileName      = NULL;
  parm->ComputeClumpingFactor  = FALSE;
  parm->MetaData               = NULL;
  parm->DiskRadiusCutoff       = BoxSize;
  parm->LinearProfileRadiusForVertical = TRUE; 
  parm->PrintGlobalProfileValues = FALSE; 

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    center[dim] = FLOAT_UNDEFINED;

  /* Open file. */

  FILE *fptr;
  if ((fptr = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "error opening file %s\n", filename);
    exit(EXIT_FAILURE);
  }

  /* Read file. */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    ret += sscanf(line, "Rinner = %"PSYM, &parm->rinner);
    ret += sscanf(line, "Router = %"PSYM, &parm->router);
    ret += sscanf(line, "CenterPosition = %"PSYM" %"PSYM" %"PSYM, 
		  &center[0], &center[1], &center[2]);
    ret += sscanf(line, "NumberOfPoints = %"ISYM, &parm->npoints);
    ret += sscanf(line, "VirialDensity = %f", &parm->virial_dens);
    ret += sscanf(line, "MeanVelocityVirialFraction = %f",
		  &parm->MeanVelocityVirialFraction);
    ret += sscanf(line, "ColdTemperatureCutoff = %f",
		  &parm->ColdTemperatureCutoff);
    ret += sscanf(line, "ColdTemperatureCutoffVirialFraction = %f",
		  &parm->ColdTemperatureCutoffVirialFraction);
    ret += sscanf(line, "VirialTemperatureNormalization = %f",
		  &parm->VirialTemperatureNormalization);
    ret += sscanf(line, "LowerDensityCutoff = %f", 
		  &parm->LowerDensityCutoff);
    ret += sscanf(line, "UpperDensityCutoff = %f", 
		  &parm->UpperDensityCutoff);
    ret += sscanf(line, "ComputeDiskInformation = %d", 
		  &parm->ComputeDiskInformation);
    ret += sscanf(line, "DiskImageSize = %d", &parm->DiskImageSize);
    ret += sscanf(line, "DiskRadius = %f", &parm->DiskRadius);
    ret += sscanf(line, "XrayLowerCutoffkeV = %f", 
		  &parm->XrayLowerCutoffkeV);
    ret += sscanf(line, "XrayUpperCutoffkeV = %f", 
		  &parm->XrayUpperCutoffkeV);
    ret += sscanf(line, "ComputeClumpingFactor = %d", 
		  &parm->ComputeClumpingFactor);
    ret += sscanf(line, "DiskRadiusCutoff = %f", &parm->DiskRadiusCutoff);
    ret += sscanf(line, "LinearProfileRadiusForVertical = %d", 
		  &parm->LinearProfileRadiusForVertical);
    ret += sscanf(line, "PrintGlobalProfileValues = %d", 
		  &parm->PrintGlobalProfileValues);

    if (sscanf(line, "CenterListName = %s", char_dummy) == 1) 
      CenterListName = char_dummy;

    if (sscanf(line, "XrayTableFileName = %s", char_dummy) == 1) 
      parm->XrayTableFileName = char_dummy;

    if (sscanf(line, "MetaData = %s", char_dummy) == 1) 
      parm->MetaData = char_dummy;

    /* If the dummy char space was used, then make another. */

    if (*char_dummy != 0) {
      char_dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") != NULL && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s", line);

  }

  /* Error check. */ 

    if (parm->rinner/BoxSize > 1 || parm->router/BoxSize > 1 || (parm->router-parm->rinner <= 0.)) {    
    fprintf(stderr, "Rinner or Router > BoxSize.\n  %g %g", parm->rinner, parm->router);
    exit(EXIT_FAILURE);
  }  


  /* Close file. */

  fclose(fptr);

  /* If a CenterListName was specified, read this file, otherwise copy
     Center into CenterList. */

  if (CenterListName == NULL) {

    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CenterList[dim] = new FLOAT[1];
      CenterList[dim][0] = center[dim];
      //fprintf(stderr, "CenterList[dim][0] = %g, center[dim] = %g", CenterList[dim][0], center[dim]);
    }
    NumberOfCenters = 1;

  } else {

    /* Open CenterListName. */

    if ((fptr = fopen(CenterListName, "r")) == NULL) {
      fprintf(stderr, "error opening CenterListName %s\n", CenterListName);
      exit(EXIT_FAILURE);
    }

    /* Count lines and allocate space. */
	
    NumberOfCenters = 0;
    while (fscanf(fptr, "%"PSYM, &float_dummy) > 0) {
      NumberOfCenters++;
      fgets(char_dummy, MAX_LINE_LENGTH, fptr);
    }
    if (debug) printf("NumberOfCenters = %d\n", NumberOfCenters);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      CenterList[dim] = new FLOAT[NumberOfCenters];
    rewind(fptr);

    /* Read data. */

    j = 0;
    while (fscanf(fptr, "%"PSYM" %"PSYM" %"PSYM, CenterList[0]+j, 
		  CenterList[1]+j, CenterList[2]+j) == 3) {
      fgets(char_dummy, MAX_LINE_LENGTH, fptr);  // get rid of rest of line
      j++;
    }
    if (j != NumberOfCenters) {
      fprintf(stderr, "Counting error (%d/%d) in %s\n", j, 
	      NumberOfCenters, CenterListName);
      exit(EXIT_FAILURE);
    }
    
  } // end: if (CenterListName == NULL)

  return SUCCESS;
}
