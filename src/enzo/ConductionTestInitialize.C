////////////////////////////////////////////////////////////////////////////////
//
//  Conduction Test Problem
//
//  written by: David A. Ventimiglia, Brian O'Shea
//  date:       June 2009
//  modified:   February 2011 by BWO
//
//  PURPOSE: Initializes a temperature pulse in a uniform medium.  Depending
//     on which version of the problem type is chosen, the hydro solver will be
//     turned on or off.
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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

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
  float ConductionTestTemperature = 1.0;
  float ConductionTestTotalEnergy = 1.0;
  float ConductionTestGasEnergy = 1.0;
  float ConductionTestVelocity[3] = {0.0,0.0,0.0};
  FLOAT ConductionTestPulseCenter[3] = {0.5,0.5,0.5};
  float ConductionTestInitialUniformBField[3] = {0.0,0.0,0.0};  // in Gauss

  for (int i = 0; i<MetaData.TopGridRank; i++) {ConductionTestVelocity[i] = 0.0;}
 
  float ConductionTestPulseHeight;
  FLOAT ConductionTestPulseWidth;
  int ConductionTestPulseType = 0, ConductionTestFieldGeometry=0;

  // Read parameters
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;
    ret += sscanf(line, "ConductionTestPulseHeight = %"FSYM, &ConductionTestPulseHeight);
    ret += sscanf(line, "ConductionTestTotalEnergy = %"FSYM, &ConductionTestTotalEnergy);
    ret += sscanf(line, "ConductionTestTemperature = %"FSYM, &ConductionTestTemperature);
    ret += sscanf(line, "ConductionTestDensity = %"FSYM, &ConductionTestDensity);
    ret += sscanf(line, "ConductionTestPulseWidth = %"PSYM, &ConductionTestPulseWidth);
    ret += sscanf(line, "ConductionTestPulseCenter = %"PSYM" %"PSYM" %"PSYM, &ConductionTestPulseCenter[0],
		  &ConductionTestPulseCenter[1], &ConductionTestPulseCenter[2]);
    ret += sscanf(line, "ConductionTestPulseType = %"ISYM, &ConductionTestPulseType);
    ret += sscanf(line, "ConductionTestFieldGeometry = %"ISYM, &ConductionTestFieldGeometry);
    ret += sscanf(line, "ConductionTestPulseBFieldX = %"FSYM,&ConductionTestInitialUniformBField[0]);
    ret += sscanf(line, "ConductionTestPulseBFieldY = %"FSYM,&ConductionTestInitialUniformBField[1]);
    ret += sscanf(line, "ConductionTestPulseBFieldZ = %"FSYM,&ConductionTestInitialUniformBField[2]);
    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    if (ret == 0 && 
	strstr(line, "=") && strstr(line, "ConductionTest") &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "*** warning: the following parameter line was not interpreted:\n%s\n", line);
    }
  }

  float BFieldVal=0.0;

  if(ConductionTestFieldGeometry >= 1){
    BFieldVal = fabs(ConductionTestInitialUniformBField[0]);

    ConductionTestInitialUniformBField[0] =
      ConductionTestInitialUniformBField[1] =
      ConductionTestInitialUniformBField[2]=0.0;
    
  }

  float DensityUnits=1.0, LengthUnits=1.0, TemperatureUnits=1.0, TimeUnits=1.0,
    VelocityUnits=1.0;
  double MassUnits=1.0;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, 0.0) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  float Boltzmann = 1.38e-16, mu = 0.6, mh=1.67e-24;

  // User can set 
  if(ConductionTestTemperature > 1.0){

    ConductionTestTotalEnergy = (Boltzmann*ConductionTestTemperature)/((Gamma - 1.0)*mu*mh);
    ConductionTestTotalEnergy /= (VelocityUnits*VelocityUnits);
    ConductionTestGasEnergy = ConductionTestTotalEnergy;
    printf("ConductionTestTotalEnergy is %e and ConductionTestTemperature is %e\n\n",ConductionTestTotalEnergy, ConductionTestTemperature);
    fflush(stdout);
  }

  if (HydroMethod == MHD_RK){
    float MagneticUnits = sqrt(DensityUnits*4.0*M_PI)*VelocityUnits;
    ConductionTestInitialUniformBField[0] /= MagneticUnits;
    ConductionTestInitialUniformBField[1] /= MagneticUnits;
    ConductionTestInitialUniformBField[2] /= MagneticUnits;

    ConductionTestTotalEnergy += 0.5*(ConductionTestInitialUniformBField[0]*ConductionTestInitialUniformBField[0] + 
				      ConductionTestInitialUniformBField[1]*ConductionTestInitialUniformBField[1] + 
				      ConductionTestInitialUniformBField[2]*ConductionTestInitialUniformBField[2])/ConductionTestDensity;
  }

  // Create a uniform grid
  if (TopGrid.GridData->InitializeUniformGrid(ConductionTestDensity,
					      ConductionTestTotalEnergy,
					      ConductionTestGasEnergy,
					      ConductionTestVelocity,
					      ConductionTestInitialUniformBField) == FAIL) {
    ENZO_FAIL("Error in InitializeUniformGrid.");
  }
  
  // Then perturb it
  if (TopGrid.GridData->ConductionTestInitialize(ConductionTestPulseHeight, 
						 ConductionTestPulseWidth, 
						 ConductionTestPulseType,
						 ConductionTestPulseCenter,
						 ConductionTestFieldGeometry, 
						 BFieldVal) == FAIL) {
    ENZO_FAIL("Error in ConductionTestInitialize.");
  }

  // set up field names and units
  i = 0;
  DataLabel[i++] = "Density";
  DataLabel[i++] = "TotalEnergy";
  if (DualEnergyFormalism) {DataLabel[i++] = "GasEnergy";}
  if (MetaData.TopGridRank > 0) {DataLabel[i++] = "x-velocity";}
  if (MetaData.TopGridRank > 1 || HydroMethod > 2) {DataLabel[i++] = "y-velocity";}
  if (MetaData.TopGridRank > 2 || HydroMethod > 2) {DataLabel[i++] = "z-velocity";}
  if (HydroMethod == MHD_RK) {
    DataLabel[i++] = "Bx";
    DataLabel[i++] = "By";
    DataLabel[i++] = "Bz";
    DataLabel[i++] = "Phi";
    if(UseDivergenceCleaning){
      DataLabel[i++] = "Phip";
    }
  }

  if (TestProblemData.UseMetallicityField)
    DataLabel[i++] = "Metal_Density";

  for (j=0; j < i; j++) 
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "ConductionTestPulseHeight = %"FSYM"\n", ConductionTestPulseHeight);
    fprintf(Outfptr, "ConductionTestPulseWidth = %"PSYM"\n", ConductionTestPulseWidth);
    fprintf(Outfptr, "ConductionTestPulseType = %"ISYM"\n", ConductionTestPulseType);

    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

  }

  if(debug){
    printf("Exiting ConductionTestInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;

}
