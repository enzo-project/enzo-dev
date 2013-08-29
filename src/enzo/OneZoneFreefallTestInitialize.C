/***********************************************************************
/
/  INITIALIZE ONE-ZONE FREE-FALL TEST PROBLEM
/
/  written by: Britton Smith
/  date:       JANUARY 2009
/  modified1:  
/
/  PURPOSE: Test chemistry and cooling in free-fall collapse.
/
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
 
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
	     float *VelocityUnits, FLOAT Time);
 
int OneZoneFreefallTestInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
				  TopGridData &MetaData)
{

  fprintf(stderr,"Initializing one-zone free-fall test.\n");

  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";
  char *MetalName = "Metal_Density";
  char *ExtraNames[2] = {"Z_Field1", "Z_Field2"};

   /* parameter declarations */
 
  FLOAT ConstantDensitySubgridLeft, ConstantDensitySubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  float ConstantDensityVelocity[3]   = {0.0, 0.0, 0.0};

  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];

  MaximumRefinementLevel = 0;

  float dx = (DomainRightEdge[0] - DomainLeftEdge[0])/
     MetaData.TopGridDims[0];

  float OneZoneFreefallTestInitialDensity = 1.0;
  float OneZoneFreefallTestMinimumEnergy = 10.0;
  float OneZoneFreefallTestMaximumEnergy = 1000.0;
  float OneZoneFreefallTestMinimumMetallicity = 1e-6;
  float OneZoneFreefallTestMaximumMetallicity = 1e-2;
  TestProblemData.OneZoneFreefallTimestepFraction = 1e-3;
  TestProblemData.OneZoneFreefallUseEffectiveGamma = FALSE;

//   /* set no subgrids by default. */
 
  ConstantDensitySubgridLeft         = 0.0;    // start of subgrid(s)
  ConstantDensitySubgridRight        = 0.0;    // end of subgrid(s)

  TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)

  /* read input from file */

  int comment_count = 0;
 
  while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL) 
      && (comment_count < 2)) {
 
    ret = 0;
 
    /* read parameters specifically for constant density problem */

    /* read in more general test parameters to set species, turn on color fields, etc. */
    ret += sscanf(line, "OneZoneFreefallTestInitialDensity = %"FSYM, &OneZoneFreefallTestInitialDensity);
    ret += sscanf(line, "OneZoneFreefallTestMinimumEnergy = %"FSYM, &OneZoneFreefallTestMinimumEnergy);
    ret += sscanf(line, "OneZoneFreefallTestMaximumEnergy = %"FSYM, &OneZoneFreefallTestMaximumEnergy);
    ret += sscanf(line, "OneZoneFreefallTestMinimumMetallicity = %"FSYM, &OneZoneFreefallTestMinimumMetallicity);
    ret += sscanf(line, "OneZoneFreefallTestMaximumMetallicity = %"FSYM, &OneZoneFreefallTestMaximumMetallicity);
    ret += sscanf(line, "OneZoneFreefallTimestepFraction = %"FSYM, 
		  &TestProblemData.OneZoneFreefallTimestepFraction);
    ret += sscanf(line, "OneZoneFreefallUseEffectiveGamma = %"ISYM,
                  &TestProblemData.OneZoneFreefallUseEffectiveGamma);

    ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);
    ret += sscanf(line, "TestProblemDeuteriumToHydrogenRatio = %"FSYM, &TestProblemData.DeuteriumToHydrogenRatio);
    ret += sscanf(line, "TestProblemInitialHIFraction  = %"FSYM, &TestProblemData.HI_Fraction);
    ret += sscanf(line, "TestProblemInitialHIIFraction  = %"FSYM, &TestProblemData.HII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIFraction  = %"FSYM, &TestProblemData.HeI_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIFraction  = %"FSYM, &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIIIFraction  = %"FSYM, &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "TestProblemInitialHMFraction  = %"FSYM, &TestProblemData.HM_Fraction);
    ret += sscanf(line, "TestProblemInitialH2IFraction  = %"FSYM, &TestProblemData.H2I_Fraction);
    ret += sscanf(line, "TestProblemInitialH2IIFraction  = %"FSYM, &TestProblemData.H2II_Fraction);
    ret += sscanf(line, "TestProblemInitialDIFraction  = %"FSYM, &TestProblemData.DI_Fraction);
    ret += sscanf(line, "TestProblemInitialDIIFraction  = %"FSYM, &TestProblemData.DII_Fraction);
    ret += sscanf(line, "TestProblemInitialHDIFraction  = %"FSYM, &TestProblemData.HDI_Fraction);
    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);

    if (strstr(line, "\"\"\"")              ) comment_count++;

    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && (strstr(line, "CoolingDensity") || strstr(line, "TestProblem")) &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,
	      "*** warning: the following parameter line was not interpreted:\n%s\n",
	      line);
 
  } // end input from parameter file

  /* Set constant for analytical free-fall collapse. */
  TestProblemData.OneZoneFreefallConstant = pow(OneZoneFreefallTestInitialDensity, -0.5);

  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, 0.0) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* set the periodic boundaries */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }
 
 //  /* set up grid */
 
  if (TopGrid.GridData->OneZoneFreefallTestInitializeGrid(OneZoneFreefallTestInitialDensity,
							  OneZoneFreefallTestMinimumEnergy,
							  OneZoneFreefallTestMaximumEnergy,
							  OneZoneFreefallTestMinimumMetallicity,
							  OneZoneFreefallTestMaximumMetallicity) == FAIL) {
    ENZO_FAIL("Error in OneZoneFreefallTestInitializeGrid.");
  }
 
  /* set up field names and units -- NOTE: these absolutely MUST be in 
     the same order that they are in Grid_InitializeUniformGrids.C, or 
     else you'll find out that data gets written into incorrectly-named
     fields.  Just FYI. */

  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if(DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;

  if(MetaData.TopGridRank > 1)
    DataLabel[i++] = Vel2Name;

  if(MetaData.TopGridRank > 2)
    DataLabel[i++] = Vel3Name;

  if (TestProblemData.MultiSpecies) {
    DataLabel[i++] = ElectronName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = HIIName;
    DataLabel[i++] = HeIName;
    DataLabel[i++] = HeIIName;
    DataLabel[i++] = HeIIIName;
    if (TestProblemData.MultiSpecies > 1) {
      DataLabel[i++] = HMName;
      DataLabel[i++] = H2IName;
      DataLabel[i++] = H2IIName;
    }
    if (TestProblemData.MultiSpecies > 2) {
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
    }
  }
 
  if (TestProblemData.UseMetallicityField) {
    DataLabel[i++] = MetalName;

    if(TestProblemData.MultiMetals){
      DataLabel[i++] = ExtraNames[0];
      DataLabel[i++] = ExtraNames[1];
    }
  }
 
  for(j=0; j < i; j++)
    DataUnits[j] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "OneZoneFreefallTestInitialDensity = %"FSYM"\n", OneZoneFreefallTestInitialDensity);
    fprintf(Outfptr, "OneZoneFreefallTestMinimumEnergy = %"FSYM"\n", OneZoneFreefallTestMinimumEnergy);
    fprintf(Outfptr, "OneZoneFreefallTestMaximumEnergy = %"FSYM"\n", OneZoneFreefallTestMaximumEnergy);
    fprintf(Outfptr, "OneZoneFreefallTestMinimumMetallicity = %"FSYM"\n", OneZoneFreefallTestMinimumMetallicity);
    fprintf(Outfptr, "OneZoneFreefallTestMaximumMetallicity = %"FSYM"\n", OneZoneFreefallTestMaximumMetallicity);
    fprintf(Outfptr, "OneZoneFreefallTimestepFraction = %"FSYM"\n", TestProblemData.OneZoneFreefallTimestepFraction);
    fprintf(Outfptr, "OneZoneFreefallUseEffectiveGamma = %"ISYM"\n", TestProblemData.OneZoneFreefallUseEffectiveGamma);

    fprintf(Outfptr, "TestProblemHydrogenFractionByMass = %"FSYM"\n",   TestProblemData.HydrogenFractionByMass);
    fprintf(Outfptr, "TestProblemDeuteriumToHydrogenRatio = %"FSYM"\n", TestProblemData.DeuteriumToHydrogenRatio);
    fprintf(Outfptr, "TestProblemInitialHIFraction  = %"FSYM"\n", TestProblemData.HI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHIIFraction  = %"FSYM"\n", TestProblemData.HII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIFraction  = %"FSYM"\n", TestProblemData.HeI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIFraction  = %"FSYM"\n", TestProblemData.HeII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIIIFraction  = %"FSYM"\n", TestProblemData.HeIII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHMFraction  = %"FSYM"\n", TestProblemData.HM_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IFraction  = %"FSYM"\n", TestProblemData.H2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IIFraction  = %"FSYM"\n", TestProblemData.H2II_Fraction);
    fprintf(Outfptr, "TestProblemInitialDIFraction  = %"FSYM"\n", TestProblemData.DI_Fraction);
    fprintf(Outfptr, "TestProblemInitialDIIFraction  = %"FSYM"\n", TestProblemData.DII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHDIFraction  = %"FSYM"\n", TestProblemData.HDI_Fraction);
    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);

  } //   if (MyProcessorNumber == ROOT_PROCESSOR) 
 
  return SUCCESS;
 
}
