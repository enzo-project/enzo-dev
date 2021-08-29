
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
  const char *GravPotName = "GravPotential";
  int CollapseTestUseParticles    = FALSE;
  float CollapseTestParticleMeanDensity = FLOAT_UNDEFINED;

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set up grid */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;

    /* read parameters */

    ret += sscanf(line, "CollapseTestUseParticles = %"ISYM, 
		  &CollapseTestUseParticles);
    ret += sscanf(line, "CollapseTestParticleMeanDensity = %"FSYM,
		  &CollapseTestParticleMeanDensity);
      
    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "CollapseTest") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } 

  if (TopGrid.GridData->FDMCollapseInitializeGrid(CollapseTestUseParticles, CollapseTestParticleMeanDensity) == FAIL) {
    ENZO_FAIL("Error in FDMCollapseInitializeGrid.");
  }

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) TEName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;

  DataLabel[count++] = (char*) RePsiName;
  DataLabel[count++] = (char*) ImPsiName;
  DataLabel[count++] = (char*) FDMDensName;
  
  if(WritePotential)
	DataLabel[count++] = (char*) GravPotName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  return SUCCESS;

}
