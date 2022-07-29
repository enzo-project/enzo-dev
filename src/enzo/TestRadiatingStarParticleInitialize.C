/***********************************************************************
/
/  INITIALIZE A RADIATING STAR PARTICLE TEST
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1: Danielle Skinner
/  date:      July, 2022
/
/  PURPOSE:
/    Initialize a radiating Pop III star particle test in a uniform medium. 
/    Based on the Single Star Particle Test and Photon Test.
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

static float TestStarParticleInitialFractionHII   = 1.2e-5;
static float TestStarParticleInitialFractionHeII  = 1.0e-14;
static float TestStarParticleInitialFractionHeIII = 1.0e-17;
static float TestStarParticleInitialFractionHM    = 2.0e-9;
static float TestStarParticleInitialFractionH2I   = 2.0e-20;
static float TestStarParticleInitialFractionH2II  = 3.0e-14;

int TestRadiatingStarParticleInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			       TopGridData &MetaData,float *Initialdt)
{
  const char *DensName = "Density";   
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *MetalName = "Metal_Density";
  const char *ElectronName = "Electron_Density";
  const char *HIName    = "HI_Density";
  const char *HIIName   = "HII_Density";
  const char *HeIName   = "HeI_Density";
  const char *HeIIName  = "HeII_Density";
  const char *HeIIIName = "HeIII_Density";
  const char *HMName    = "HM_Density";
  const char *H2IName   = "H2I_Density";
  const char *H2IIName  = "H2II_Density";
  const char *DIName    = "DI_Density";
  const char *DIIName   = "DII_Density";
  const char *HDIName   = "HDI_Density";
  const char *kphHIName    = "HI_kph";
  const char *gammaName  = "PhotoGamma";
  const char *kphHeIName   = "HeI_kph";   
  const char *kphHeIIName  = "HeII_kph";
  const char *kphHMName   = "HM_kph";   
  const char *kdissH2IIName  = "H2II_kdiss";
  const char *kdissH2IName = "H2I_kdiss"; 
  const char *RadAccel1Name = "x-RadPressure";
  const char *RadAccel2Name = "y-RadPressure";
  const char *RadAccel3Name = "z-RadPressure";


  /* declarations */

  char  line[MAX_LINE_LENGTH];
  char *dummy = new char[MAX_LINE_LENGTH];
  int   dim, ret, i, source;
  dummy[0] = 0;

  /* Error check. */


  /* set default parameters */

  float TestStarParticleDensity     = 1.0;
  float TestStarParticleEnergy      = 1.0;
  float TestStarParticleVelocity[3] = {0.0, 0.0, 0.0};
  FLOAT TestStarParticleStarVelocity[3] = {0.0, 0.0, 0.0};
  FLOAT TestStarParticleStarPosition[3] = {0.5, 0.5, 0.5};
   float TestStarParticleBField[3]   = {0.0, 0.0, 0.0};
  float TestStarParticleStarMass    = 100.0;
  int TestProblemUseMetallicityField = 1;
  float TestProblemInitialMetallicityFraction = 2e-3; // 0.1 Zsun




  TestProblemData.MultiSpecies = MultiSpecies;
  TestProblemData.UseMetallicityField = TestProblemUseMetallicityField;
  TestProblemData.MetallicityField_Fraction = TestProblemInitialMetallicityFraction;

  /* read input from file */

  rewind(fptr);
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "TestStarParticleDensity = %"FSYM,
		  &TestStarParticleDensity);
    ret += sscanf(line, "TestStarParticleEnergy = %"FSYM,
		  &TestStarParticleEnergy);
    ret += sscanf(line, "TestStarParticleStarMass = %"FSYM,
		  &TestStarParticleStarMass);
    ret += sscanf(line,"TestStarParticleStarVelocity = %"PSYM" %"PSYM" %"PSYM, 
		  &TestStarParticleStarVelocity[0],
		  &TestStarParticleStarVelocity[1],
		  &TestStarParticleStarVelocity[2]);
    ret += sscanf(line,"TestStarParticleStarPosition = %"PSYM" %"PSYM" %"PSYM, 
		  &TestStarParticleStarPosition[0],
		  &TestStarParticleStarPosition[1],
		  &TestStarParticleStarPosition[2]);
    

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction); 

    //ret += sscanf(line, "TestProblemInitialHIFraction  = %"FSYM, &TestProblemData.HI_Fraction);
    //ret += sscanf(line, "TestProblemInitialHIIFraction  = %"FSYM, &TestProblemData.HII_Fraction);
    //ret += sscanf(line, "TestProblemInitialHeIFraction  = %"FSYM, &TestProblemData.HeI_Fraction);
    ////ret += sscanf(line, "TestProblemInitialHeIIFraction  = %"FSYM, &TestProblemData.HeII_Fraction);
    //ret += sscanf(line, "TestProblemInitialHeIIIIFraction  = %"FSYM, &TestProblemData.HeIII_Fraction);
    //ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);


    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "TestStarParticle") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up uniform grid as of before explosion */

  if (TopGrid.GridData->
      TestRadiatingStarParticleInitializeGrid(TestStarParticleStarMass,
				     Initialdt, 
				     TestStarParticleStarVelocity,
				     TestStarParticleStarPosition,
             TestStarParticleDensity,
             TestStarParticleEnergy, 
             TestStarParticleVelocity) == FAIL)
  ENZO_FAIL("Error in TestRadiatingStarParticleInitializeGrid.\n");

  /* set up field names and units */
  
  int count = 0;
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = (char*) GEName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;
  if (MultiSpecies){
    DataLabel[count++] = (char*) ElectronName;
    DataLabel[count++] = (char*) HIName;
    DataLabel[count++] = (char*) HIIName;
    DataLabel[count++] = (char*) HeIName;
    DataLabel[count++] = (char*) HeIIName;
    DataLabel[count++] = (char*) HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[count++] = (char*) HMName;
      DataLabel[count++] = (char*) H2IName;
      DataLabel[count++] = (char*) H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[count++] = (char*) DIName;
      DataLabel[count++] = (char*) DIIName;
      DataLabel[count++] = (char*) HDIName;
    }
  } // if Multispecies
  if (TestProblemData.UseMetallicityField)
    DataLabel[count++] = (char*) MetalName;

  if (RadiativeTransfer)
    if (MultiSpecies) {
      DataLabel[count++]  = (char*) kphHIName;
      DataLabel[count++]  = (char*) gammaName;
      DataLabel[count++]  = (char*) kphHeIName;
      DataLabel[count++]  = (char*) kphHeIIName;
      if (MultiSpecies > 1) {
	DataLabel[count++]= (char*) kdissH2IName;
	DataLabel[count++]= (char*) kdissH2IIName;
	DataLabel[count++]= (char*) kphHMName;
      }
    } // if RadiativeTransfer

  if (RadiationPressure) {
    DataLabel[count++]  = (char*) RadAccel1Name;
    DataLabel[count++]  = (char*) RadAccel2Name;
    DataLabel[count++]  = (char*) RadAccel3Name;
  }


  int j;
  for(j=0; j < count; j++)
    DataUnits[j] = NULL;

   
  /* Write parameters to parameter output file */
  
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "TestStarParticleDensity = %"FSYM"\n",
	    TestStarParticleDensity);
    fprintf(Outfptr, "TestStarParticleEnergy = %"FSYM"\n",
	    TestStarParticleEnergy);
    fprintf(Outfptr, "MetallicityField_Fraction = %"FSYM"\n",
            TestProblemData.MetallicityField_Fraction);
  }

  fprintf(stderr, "TestStarParticleDensity = %"FSYM"\n",
	  TestStarParticleDensity);
  fprintf(stderr, "TestStarParticleEnergy = %"FSYM"\n",
	  TestStarParticleEnergy);
  fprintf(stderr, "MetallicityField_Fraction = %"FSYM"\n",
	  TestProblemData.MetallicityField_Fraction);



  return SUCCESS;

}

