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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

// Problem Initializer
int ConductionCloudInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid, TopGridData &MetaData){

  if(debug){
    printf("Entering ConductionCloudInitialize\n");
    fflush(stdout);
  }

  char line[MAX_LINE_LENGTH];
  float LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  int i, j, dim, ret;

  float ConductionCloudDensity = 1.0;
  float ConductionCloudTemperature = 1.0;
  float ConductionCloudTotalEnergy = 1.0;
  float ConductionCloudVelocity[3] = {0.0,0.0,0.0};
  float ConductionCloudInitialUniformBField[3] = {0.0,0.0,0.0};  // in Gauss

  for (int i = 0; i<MetaData.TopGridRank; i++) {ConductionCloudVelocity[i] = 0.0;}
 
  float ConductionCloudPulseHeight;
  FLOAT ConductionCloudPulseWidth;
  int ConductionCloudPulseType = 0;

  TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)

  // Read parameters
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;
    ret += sscanf(line, "ConductionCloudPulseHeight = %"FSYM, &ConductionCloudPulseHeight);
    ret += sscanf(line, "ConductionCloudTotalEnergy = %"FSYM, &ConductionCloudTotalEnergy);
    ret += sscanf(line, "ConductionCloudTemperature = %"FSYM, &ConductionCloudTemperature);
    ret += sscanf(line, "ConductionCloudDensity = %"FSYM, &ConductionCloudDensity);
    ret += sscanf(line, "ConductionCloudPulseWidth = %"PSYM, &ConductionCloudPulseWidth);
    ret += sscanf(line, "ConductionCloudPulseType = %"ISYM, &ConductionCloudPulseType);

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    /* read in more general test parameters to set species, turn on color fields, etc. */
    ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);
    ret += sscanf(line, "TestProblemDeuteriumToHydrogenRatio = %"FSYM, &TestProblemData.DeuteriumToHydrogenRatio);

    ret += sscanf(line, "TestProblemInitialHIFractionInner  = %"FSYM, &TestProblemData.HI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHIIFractionInner  = %"FSYM, &TestProblemData.HII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHeIFractionInner  = %"FSYM, &TestProblemData.HeI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHeIIFractionInner  = %"FSYM, &TestProblemData.HeII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHeIIIFractionInner  = %"FSYM, &TestProblemData.HeIII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHMFractionInner  = %"FSYM, &TestProblemData.HM_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialH2IFractionInner  = %"FSYM, &TestProblemData.H2I_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialH2IIFractionInner  = %"FSYM, &TestProblemData.H2II_Fraction_Inner);

    ret += sscanf(line, "TestProblemInitialDIFractionInner  = %"FSYM, &TestProblemData.DI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialDIIFractionInner  = %"FSYM, &TestProblemData.DII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHDIFractionInner  = %"FSYM, &TestProblemData.HDI_Fraction_Inner);

    ret += sscanf(line, "TestProblemInitialHIFraction  = %"FSYM, &TestProblemData.HI_Fraction);
    ret += sscanf(line, "TestProblemInitialHIIFraction  = %"FSYM, &TestProblemData.HII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIFraction  = %"FSYM, &TestProblemData.HeI_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIFraction  = %"FSYM, &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIIFraction  = %"FSYM, &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "TestProblemInitialHMFraction  = %"FSYM, &TestProblemData.HM_Fraction);
    ret += sscanf(line, "TestProblemInitialH2IFraction  = %"FSYM, &TestProblemData.H2I_Fraction);
    ret += sscanf(line, "TestProblemInitialH2IIFraction  = %"FSYM, &TestProblemData.H2II_Fraction);

    ret += sscanf(line, "TestProblemInitialDIFraction  = %"FSYM, &TestProblemData.DI_Fraction);
    ret += sscanf(line, "TestProblemInitialDIIFraction  = %"FSYM, &TestProblemData.DII_Fraction);
    ret += sscanf(line, "TestProblemInitialHDIFraction  = %"FSYM, &TestProblemData.HDI_Fraction);

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    ret += sscanf(line, "TestProblemMultiMetals  = %"ISYM, &TestProblemData.MultiMetals);
    ret += sscanf(line, "TestProblemInitialMultiMetalsField1Fraction  = %"FSYM, &TestProblemData.MultiMetalsField1_Fraction);
    ret += sscanf(line, "TestProblemInitialMultiMetalsField2Fraction  = %"FSYM, &TestProblemData.MultiMetalsField2_Fraction);

    if (ret == 0 && 
	strstr(line, "=") && strstr(line, "ConductionCloud") &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "*** warning: the following parameter line was not interpreted:\n%s\n", line);
    }
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

  if(ConductionCloudTemperature > 1.0){

    ConductionCloudTotalEnergy = (Boltzmann*ConductionCloudTemperature)/((Gamma - 1.0)*mu*mh);
    ConductionCloudTotalEnergy /= (VelocityUnits*VelocityUnits);
    printf("ConductionCloudTotalEnergy is %e and ConductionCloudTemperature is %e\n\n",ConductionCloudTotalEnergy, ConductionCloudTemperature);
    fflush(stdout);
  }

  // Create a uniform grid
  if (TopGrid.GridData->InitializeUniformGrid(ConductionCloudDensity,
					      ConductionCloudTotalEnergy,
					      ConductionCloudTotalEnergy,
					      ConductionCloudVelocity,
					      ConductionCloudInitialUniformBField) == FAIL) {
    ENZO_FAIL("Error in InitializeUniformGrid.");
  }
  
  // Then perturb it
  if (TopGrid.GridData->ConductionCloudInitialize(ConductionCloudPulseHeight, 
						 ConductionCloudPulseWidth, 
						 ConductionCloudPulseType) == FAIL) {
    ENZO_FAIL("Error in ConductionCloudInitialize.");
  }

  // set up field names and units
  i = 0;
  DataLabel[i++] = "Density";
  DataLabel[i++] = "TotalEnergy";
  if (DualEnergyFormalism) {DataLabel[i++] = "GasEnergy";}
  if (MetaData.TopGridRank > 0) {DataLabel[i++] = "x-velocity";}
  if (MetaData.TopGridRank > 1) {DataLabel[i++] = "y-velocity";}
  if (MetaData.TopGridRank > 2) {DataLabel[i++] = "z-velocity";}

  if (TestProblemData.MultiSpecies) {
    DataLabel[i++] = "Electron_Density";
    DataLabel[i++] = "HI_Density";
    DataLabel[i++] = "HII_Density";
    DataLabel[i++] = "HeI_Density";
    DataLabel[i++] = "HeII_Density";
    DataLabel[i++] = "HeIII_Density";
    if (TestProblemData.MultiSpecies > 1) {
      DataLabel[i++] = "HM_Density";
      DataLabel[i++] = "H2I_Density";
      DataLabel[i++] = "H2II_Density";
    }
    if (TestProblemData.MultiSpecies > 2) {
      DataLabel[i++] = "DI_Density";
      DataLabel[i++] = "DII_Density";
      DataLabel[i++] = "HDI_Density";
    }
  }

  if (TestProblemData.UseMetallicityField)
    DataLabel[i++] = "Metal_Density";

  for (j=0; j < i; j++) 
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "ConductionCloudPulseHeight = %"FSYM"\n", ConductionCloudPulseHeight);
    fprintf(Outfptr, "ConductionCloudPulseWidth = %"PSYM"\n", ConductionCloudPulseWidth);
    fprintf(Outfptr, "ConductionCloudPulseType = %"ISYM"\n", ConductionCloudPulseType);

    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);


    fprintf(Outfptr, "TestProblemHydrogenFractionByMass = %"FSYM"\n",   TestProblemData.HydrogenFractionByMass);
    fprintf(Outfptr, "TestProblemDeuteriumToHydrogenRatio = %"FSYM"\n", TestProblemData.DeuteriumToHydrogenRatio);

    fprintf(Outfptr, "TestProblemInitialHIFractionInner  = %"FSYM"\n", TestProblemData.HI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHIIFractionInner  = %"FSYM"\n", TestProblemData.HII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHeIFractionInner  = %"FSYM"\n", TestProblemData.HeI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHeIIFractionInner  = %"FSYM"\n", TestProblemData.HeII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHeIIIFractionInner  = %"FSYM"\n", TestProblemData.HeIII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHMFractionInner  = %"FSYM"\n", TestProblemData.HM_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialH2IFractionInner  = %"FSYM"\n", TestProblemData.H2I_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialH2IIFractionInner  = %"FSYM"\n", TestProblemData.H2II_Fraction_Inner);

    fprintf(Outfptr, "TestProblemInitialDIFractionInner  = %"FSYM"\n", TestProblemData.DI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialDIIFractionInner  = %"FSYM"\n", TestProblemData.DII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHDIFractionInner  = %"FSYM"\n", TestProblemData.HDI_Fraction_Inner);

    fprintf(Outfptr, "TestProblemInitialHIFraction  = %"FSYM"\n", TestProblemData.HI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHIIFraction  = %"FSYM"\n", TestProblemData.HII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIFraction  = %"FSYM"\n", TestProblemData.HeI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIFraction  = %"FSYM"\n", TestProblemData.HeII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIIFraction  = %"FSYM"\n", TestProblemData.HeIII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHMFraction  = %"FSYM"\n", TestProblemData.HM_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IFraction  = %"FSYM"\n", TestProblemData.H2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IIFraction  = %"FSYM"\n", TestProblemData.H2II_Fraction);

    fprintf(Outfptr, "TestProblemInitialDIFraction  = %"FSYM"\n", TestProblemData.DI_Fraction);
    fprintf(Outfptr, "TestProblemInitialDIIFraction  = %"FSYM"\n", TestProblemData.DII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHDIFraction  = %"FSYM"\n", TestProblemData.HDI_Fraction);

    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

    fprintf(Outfptr, "TestProblemMultiMetals  = %"ISYM"\n", TestProblemData.MultiMetals);
    fprintf(Outfptr, "TestProblemInitialMultiMetalsField1Fraction  = %"FSYM"\n", TestProblemData.MultiMetalsField1_Fraction);
    fprintf(Outfptr, "TestProblemInitialMultiMetalsField2Fraction  = %"FSYM"\n", TestProblemData.MultiMetalsField2_Fraction);

  }

  if(debug){
    printf("Exiting ConductionCloudInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;

}
