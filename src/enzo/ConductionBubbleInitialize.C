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
  float ConductionBubbleInitialUniformBField[3] = {0.0,0.0,0.0};  // in Gauss

  FLOAT ConductionBubbleRadiusOfBubble = 0.1;  // units of box size
  int   ConductionBubblePulseType = 1;  // pulse type
  float ConductionBubbleDeltaEntropy = 0.1;    // multiple of entropy
  float ConductionBubbleMidpointEntropy = 20.0;  // entropy in kev*cm^2
  float ConductionBubbleEntropyGradient = 1.0;   // kev*cm^2 / kpc
  float ConductionBubbleMidpointTemperature = 5.0e+7;  // Kelvin
  FLOAT ConductionBubbleCenter[MAX_DIMENSION] = {0.5,0.5,0.5};

  TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)
 
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
    printf("Exiting ConductionBubbleInitialize\n");
    fflush(stdout);
  }

  return SUCCESS;

}
