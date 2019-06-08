
/***********************************************************************
/
/  INITIALIZE A COLLAPSE TEST
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:
/
/  PURPOSE:
/    Set up a number of spherical objects
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include <stdlib.h>
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
#include "LevelHierarchy.h"
#include "TopGridData.h"

int FDMCollapseInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  const char *DensName = "Density";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *ColourName = "SN_Colour";
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
  const char *MetalName = "Metal_Density";
  const char *RePsiName = "Re_Psi"; 
  const char *ImPsiName = "Im_Psi"; 
  const char *FDMDensName = "FDMDensity"; 

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set up grid */

  if (TopGrid.GridData->FDMCollapseInitializeGrid() == FAIL) {
    ENZO_FAIL("Error in CollapseTestInitializeGrid.");
  }


  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) TEName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;

  if (QuantumPressure) {
    DataLabel[count++] = (char*) RePsiName;
    DataLabel[count++] = (char*) ImPsiName;
    DataLabel[count++] = (char*) FDMDensName;
  }

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  /*if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "CollapseTestNumberOfSpheres    = %"ISYM"\n",
	    CollapseTestNumberOfSpheres);
    fprintf(Outfptr, "CollapseTestRefineAtStart      = %"ISYM"\n",
	    CollapseTestRefineAtStart);
    fprintf(Outfptr, "CollapseTestUseParticles       = %"ISYM"\n",
	    CollapseTestUseParticles);
    fprintf(Outfptr, "CollapseTestUseColour          = %"ISYM"\n",
	    CollapseTestUseColour);
    fprintf(Outfptr, "CollapseTestUseMetals          = %"ISYM"\n",
	    CollapseTestUseMetals);
    fprintf(Outfptr, "CollapseTestInitialTemperature = %"FSYM"\n",
	    CollapseTestInitialTemperature);
    fprintf(Outfptr, "CollapseTestInitialDensity     = %"FSYM"\n",
	    CollapseTestInitialDensity);
    fprintf(Outfptr, "CollapseTestUniformVelocity    = %"FSYM" %"FSYM" %"FSYM"\n",
	    CollapseTestUniformVelocity[0], CollapseTestUniformVelocity[1],
	    CollapseTestUniformVelocity[2]);
    for (sphere = 0; sphere < CollapseTestNumberOfSpheres; sphere++) {
      fprintf(Outfptr, "CollapseTestSphereType[%"ISYM"] = %"ISYM"\n", sphere,
	      CollapseTestSphereType[sphere]);
      fprintf(Outfptr, "CollapseTestSphereRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereRadius[sphere]);
      fprintf(Outfptr, "CollapseTestSphereCoreRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereCoreRadius[sphere]);
      fprintf(Outfptr, "CollapseTestSphereDensity[%"ISYM"] = %"FSYM"\n", sphere,
	      CollapseTestSphereDensity[sphere]);
      fprintf(Outfptr, "CollapseTestSphereTemperature[%"ISYM"] = %"FSYM"\n", sphere,
	      CollapseTestSphereTemperature[sphere]);
      fprintf(Outfptr, "CollapseTestSphereMetallicity[%"ISYM"] = %"FSYM"\n", sphere,
	      CollapseTestSphereMetallicity[sphere]);
      fprintf(Outfptr, "CollapseTestSpherePosition[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CollapseTestSpherePosition[sphere]);
      fprintf(Outfptr, "CollapseTestSphereVelocity[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CollapseTestSphereVelocity[sphere]);
      fprintf(Outfptr, "CollapseTestFracKeplerianRot[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestFracKeplerianRot[sphere]);
      fprintf(Outfptr, "CollapseTestSphereTurbulence[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestSphereTurbulence[sphere]);
      fprintf(Outfptr, "CollapseTestSphereCutOff[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestSphereCutOff[sphere]);
      fprintf(Outfptr, "CollapseTestSphereAng1[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestSphereAng1[sphere]);
      fprintf(Outfptr, "CollapseTestSphereAng2[%"ISYM"] = %"GOUTSYM"\n", sphere,
              CollapseTestSphereAng2[sphere]);
      fprintf(Outfptr, "CollapseTestSphereNumShells[%"ISYM"] = %"ISYM"\n", sphere,
              CollapseTestSphereNumShells[sphere]);
      fprintf(Outfptr, "CollapseTestSphereConstantPressure[%"ISYM"] = %"ISYM"\n", sphere,
	      CollapseTestSphereConstantPressure[sphere]);
      fprintf(Outfptr, "CollapseTestSphereSmoothSurface[%"ISYM"] = %"ISYM"\n", sphere,
	      CollapseTestSphereSmoothSurface[sphere]);
      fprintf(Outfptr, "CollapseTestSphereSmoothRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereSmoothRadius[sphere]);
      fprintf(Outfptr, "CollapseTestSphereHIIFraction[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereHIIFraction[sphere]);
      fprintf(Outfptr, "CollapseTestSphereHeIIFraction[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereHeIIFraction[sphere]);
      fprintf(Outfptr, "CollapseTestSphereHeIIIFraction[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereHeIIIFraction[sphere]);
      fprintf(Outfptr, "CollapseTestSphereH2IFraction[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      CollapseTestSphereH2IFraction[sphere]);
    }
  }*/

  return SUCCESS;

}
