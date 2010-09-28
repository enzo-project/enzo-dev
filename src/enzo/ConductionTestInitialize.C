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
int ConductionTestInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData){

  if(debug){
    printf("Entering ConductionTestInitialize\n");
    fflush(stdout);
  }

  char line[MAX_LINE_LENGTH];
  float LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  int i, j, dim, ret;

  float ConductionTestDensity = 1.0;
  float ConductionTestTotalEnergy = 1.0;
  float ConductionTestVelocity[3] = {0.0,0.0,0.0};
  float ConductionTestInitialUniformBField[MAX_DIMENSION];  // in Gauss

  for (int i = 0; i<MetaData.TopGridRank; i++) {ConductionTestVelocity[i] = 0.0;}
 
  // Read parameters
  float ConductionTestPulseHeight;
  FLOAT ConductionTestPulseWidth;
  int ConductionTestPulseType = 0;

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;
    ret += sscanf(line, "ConductionTestPulseHeight = %"FSYM, &ConductionTestPulseHeight);
    ret += sscanf(line, "ConductionTestPulseWidth = %"PSYM, &ConductionTestPulseWidth);
    ret += sscanf(line, "ConductionTestPulseType = %"ISYM, &ConductionTestPulseType);
    if (ret == 0 && 
	strstr(line, "=") && strstr(line, "ConductionTest") &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "*** warning: the following parameter line was not interpreted:\n%s\n", line);
    }
  }

  // Create a uniform grid
  if (TopGrid.GridData->InitializeUniformGrid(ConductionTestDensity,
					      ConductionTestTotalEnergy,
					      ConductionTestTotalEnergy,
					      ConductionTestVelocity,
					      ConductionTestInitialUniformBField) == FAIL) {
    ENZO_FAIL("Error in InitializeUniformGrid.");
  }
  
  // Then perturb it
  if (TopGrid.GridData->ConductionTestInitialize(ConductionTestPulseHeight, 
						 ConductionTestPulseWidth, 
						 ConductionTestPulseType) == FAIL) {
    ENZO_FAIL("Error in ConductionTestInitialize.");
  }

  // set up field names and units
  i = 0;
  DataLabel[i++] = "Density";
  DataLabel[i++] = "Total_Energy";
  if (DualEnergyFormalism) {DataLabel[i++] = "Gas_Energy";}
  if (MetaData.TopGridRank > 0) {DataLabel[i++] = "x-velocity";}
  if (MetaData.TopGridRank > 1) {DataLabel[i++] = "y-velocity";}
  if (MetaData.TopGridRank > 2) {DataLabel[i++] = "z-velocity";}
  for (j=0; j < i; j++) 
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "ConductionTestPulseHeight = %"FSYM"\n", ConductionTestPulseHeight);
    fprintf(Outfptr, "ConductionTestPulseWidth = %"PSYM"\n", ConductionTestPulseWidth);
    fprintf(Outfptr, "ConductionTestPulseType = %"ISYM"\n", ConductionTestPulseType);
  }

  if(debug){
    printf("Exiting ConductionTestInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;

}
