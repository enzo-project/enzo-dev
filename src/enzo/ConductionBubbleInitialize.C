////////////////////////////////////////////////////////////////////////////////
//
//  Conduction Test Problem
//
//  written by: David A. Ventimiglia, Brian O'Shea
//  date:       June 2009
//  modified:  
//
//  PURPOSE: 
//
//  RETURNS: SUCCESS or FAIL
//
////////////////////////////////////////////////////////////////////////////////
 
#include <string.h>
#include <stdio.h>
#include <math.h>
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

// Problem Initializer
int ConductionBubbleInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData){

  if(debug){
    printf("Entering ConductionBubbleInitialize\n");
    fflush(stdout);
  }

  char line[MAX_LINE_LENGTH];
  float LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  int i, j, dim, ret;

  float ConductionBubbleDensity = 1.0;
  float ConductionBubbleTotalEnergy = 1.0;
  float ConductionBubbleVelocity[3] = {0.0,0.0,0.0};
  float ConductionBubbleInitialUniformBField[MAX_DIMENSION];  // in Gauss

  FLOAT ConductionBubbleRadiusOfBubble = 0.1;  // units of box size
  int   ConductionBubblePulseType = 1;  // pulse type
  float ConductionBubbleDeltaEntropy = 0.1;    // multiple of entropy
  float ConductionBubbleMidpointEntropy = 20.0;  // entropy in kev*cm^2
  float ConductionBubbleEntropyGradient = 1.0;   // kev*cm^2 / kpc
  float ConductionBubbleMidpointTemperature = 5.0e+7;  // Kelvin
  FLOAT ConductionBubbleCenter[MAX_DIMENSION] = {0.5,0.5,0.5};


  // Read parameters

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;
    ret += sscanf(line, "ConductionBubbleRadiusOfBubble = %"PSYM, &ConductionBubbleRadiusOfBubble);
    ret += sscanf(line, "ConductionBubblePulseType = %"ISYM, &ConductionBubblePulseType);
    ret += sscanf(line, "ConductionBubbleDeltaEntropy = %"FSYM, &ConductionBubbleDeltaEntropy);
    ret += sscanf(line, "ConductionBubbleMidpointEntropy = %"FSYM, &ConductionBubbleMidpointEntropy);
    ret += sscanf(line, "ConductionBubbleEntropyGradient = %"FSYM, &ConductionBubbleEntropyGradient);
    ret += sscanf(line, "ConductionBubbleMidpointTemperature = %"FSYM, &ConductionBubbleMidpointTemperature);

    ret += sscanf(line, "ConductionBubbleCenter = %"PSYM" %"PSYM" %"PSYM, &ConductionBubbleCenter[0],
		  &ConductionBubbleCenter[1],&ConductionBubbleCenter[2]);
    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    if (ret == 0 && 
	strstr(line, "=") && strstr(line, "ConductionBubble") &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "*** warning: the following parameter line was not interpreted:\n%s\n", line);
    }
  }

  // Create a uniform grid
  if (TopGrid.GridData->InitializeUniformGrid(ConductionBubbleDensity,
					      ConductionBubbleTotalEnergy,
					      ConductionBubbleTotalEnergy,
					      ConductionBubbleVelocity,
					      ConductionBubbleInitialUniformBField) == FAIL) {
    ENZO_FAIL("Error in InitializeUniformGrid.");
  }
  
  // Then perturb it
  if (TopGrid.GridData->ConductionBubbleInitialize(ConductionBubbleRadiusOfBubble, 
						   ConductionBubblePulseType,
						   ConductionBubbleDeltaEntropy, 
						   ConductionBubbleMidpointEntropy, 
						   ConductionBubbleEntropyGradient, 
						   ConductionBubbleMidpointTemperature,
						   ConductionBubbleCenter) == FAIL) {
    ENZO_FAIL("Error in ConductionBubbleInitialize.");
  }

  // set up field names and units
  i = 0;
  DataLabel[i++] = "Density";
  DataLabel[i++] = "Total_Energy";
  if (DualEnergyFormalism) {DataLabel[i++] = "Gas_Energy";}
  if (MetaData.TopGridRank > 0) {DataLabel[i++] = "x-velocity";}
  if (MetaData.TopGridRank > 1) {DataLabel[i++] = "y-velocity";}
  if (MetaData.TopGridRank > 2) {DataLabel[i++] = "z-velocity";}

  if (TestProblemData.UseMetallicityField)
    DataLabel[i++] = "Metal_Density";

  for (j=0; j < i; j++) 
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "ConductionBubbleRadiusOfBubble = %"PSYM"\n", ConductionBubbleRadiusOfBubble);
    fprintf(Outfptr, "ConductionBubblePulseType = %"ISYM"\n", ConductionBubblePulseType);
    fprintf(Outfptr, "ConductionBubbleDeltaEntropy = %"FSYM"\n", ConductionBubbleDeltaEntropy);
    fprintf(Outfptr, "ConductionBubbleMidpointEntropy = %"FSYM"\n", ConductionBubbleMidpointEntropy);
    fprintf(Outfptr, "ConductionBubbleEntropyGradient = %"FSYM"\n", ConductionBubbleEntropyGradient);
    fprintf(Outfptr, "ConductionBubbleMidpointTemperature = %"FSYM"\n", ConductionBubbleMidpointTemperature);
    fprintf(Outfptr, "ConductionBubbleCenter = %"PSYM" %"PSYM" %"PSYM"\n", ConductionBubbleCenter,
		  ConductionBubbleCenter+1,ConductionBubbleCenter+2);
    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);
  }

  if(debug){
    printf("Exiting ConductionBubbleInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;

}
