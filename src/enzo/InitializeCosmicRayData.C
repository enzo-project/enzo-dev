/***********************************************************************
/
/  INITIALIZE THE COSMIC RAY ACCELERATION EFFICIENCIES
/
/  written by: Sam Skillman
/  date:       June, 2008
/  modified1:
/
/  PURPOSE: Initializes the cosmic ray efficiency table that is used to
/  interpolate the pre-shock CR component and Mach number
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
int InitializeCosmicRayData(void)
{
  /* Open input file for data. */
 
  FILE *fptr = fopen("cosmic_ray.dat", "r");
  if (fptr == NULL) {
    fprintf(stderr, "Error opening cosmic_ray.dat\n");
    return FAIL;
  }
 
  /* Read CR efficiency data, skipping over comments. */
  //Pre-CR is in normal, Mach in log.
  int index = 0;
  int i, j;
  float precr, mach, creff;
  char line[MAX_LINE_LENGTH];
  CosmicRayData.CRNumberPrePopValues = 2; //Pre-CR
  CosmicRayData.CRNumberMachValues = 512; //Mach

  CosmicRayData.CREfficiency =     
    new float*[CosmicRayData.CRNumberPrePopValues];
  for(i=0; i < CosmicRayData.CRNumberPrePopValues; i++)
    CosmicRayData.CREfficiency[i] = new float[CosmicRayData.CRNumberMachValues];

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
    if (line[0] != '#')
      if (sscanf(line, "%"FSYM" %"FSYM" %"FSYM, &precr, &mach,
		 &creff) == 3) {
	
	i = index / CosmicRayData.CRNumberMachValues;
	j = index - (i * CosmicRayData.CRNumberMachValues);
	
	CosmicRayData.CREfficiency[i][j] = creff;

	if (index == 0){
	  CosmicRayData.CRMinPrePop = precr;
	  CosmicRayData.CRMinMach   = log10(mach);
	}
	if (index == (CosmicRayData.CRNumberMachValues -1))
	  CosmicRayData.CRMaxMach   = log10(mach);
	if (index == (CosmicRayData.CRNumberMachValues*
		      CosmicRayData.CRNumberPrePopValues -1))
	  CosmicRayData.CRMaxPrePop = precr;
	
	index++;
      }
  fclose(fptr);

  if (index < (CosmicRayData.CRNumberPrePopValues * 
	       CosmicRayData.CRNumberMachValues - 1)){
    printf("Number of entries smaller than expected: %i\n",index);
    return FAIL;
  }

  if (debug) {
    printf("InitializeCosmicRayData: MinPrePop = %"GSYM"\n",
	   CosmicRayData.CRMinPrePop);
    printf("InitializeCosmicRayData: MaxPrePop = %"GSYM"\n",
	   CosmicRayData.CRMaxPrePop);
    printf("InitializeCosmicRayData: MinMach = %"GSYM"\n",
	   CosmicRayData.CRMinMach);
    printf("InitializeCosmicRayData: MaxMach = %"GSYM"\n",
	   CosmicRayData.CRMaxMach);
    printf("InitializeCosmicRayData: Lines Read = %"ISYM"\n",
	   index);
  }
 
  return SUCCESS;
}
