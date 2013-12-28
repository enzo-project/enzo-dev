/***********************************************************************
/
/  WRITE OUT GRID DATA AT THE LOCATION OF ANY TRACER PARTICLES
/
/  written by: Greg Bryan
/  date:       March, 2004
/  modified1:  BWO, Dec. 2013: tracer particles now work for any
/              dimensionality.  also converted velocity output to
/              a runtime parameter.
/
/  PURPOSE: This function walks through the grids and writes out,
/     for each grid, the position, density and temperature at the
/     location of tracer particles which are advected at the grid
/     velocity.  These follow stream-lines through the gas.
/     Each processor creates its own file to eliminate communication.
/
************************************************************************/
 
 
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
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
 
static char SummarySuffix[] = ".summary";
 
int WriteTracerParticleData(char *basename, int dumpnumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime)
{
  /* Exit if tracer particles not turned on. */
 
  if (TracerParticleOn == FALSE)
    return SUCCESS;
 
  /* declarations */
 
  char id[7], name[MAX_LINE_LENGTH];
  FILE *fptr = NULL;
  int level, Zero = 0;
 
  /* Compute redshift and units. */
 
  FLOAT a = 1, dadt, Redshift = 0;
  float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
    TemperatureUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, WriteTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(WriteTime, &a, &dadt);
    Redshift = (1.0+InitialRedshift)/a - 1;
  }
 
  /* Create output filename */
 
  strcpy(name, basename);
  sprintf(id, "%-d", MyProcessorNumber);  /* create processor # */
  strcat(name, id);
 
  if (debug)
    printf("WriteTracerParticleData: writing file %s.\n", name);
 
  /* Open file. */
 
  if ((fptr = fopen(name, "ab")) == NULL) {
    ENZO_VFAIL("Error opening output file %s\n", name)
  }
 
  /* Output header information to the file :
        0) time in problem units
        1) redshift
        2) time conversion factor
        3) length conversion factor
        4) density conversion factor
        5) velocity conversion factor
        6) Number of values per tracer particle
           This is currently either 2 + simulation dimensionality (density and temperature
	   plus position), or 2 + 2*simulation dimensionality (as before, but with velocity).  */
 
  float float_temp = float(WriteTime);  // convert from FLOAT to float
  fwrite((void*) &float_temp, sizeof(float), 1, fptr);
  float_temp = float(Redshift);
  fwrite((void*) &float_temp, sizeof(float), 1, fptr);
  fwrite((void*) &TimeUnits, sizeof(float), 1, fptr);
  fwrite((void*) &LengthUnits, sizeof(float), 1, fptr);
  fwrite((void*) &DensityUnits, sizeof(float), 1, fptr);
  fwrite((void*) &VelocityUnits, sizeof(float), 1, fptr);

  int int_temp;
  

  int_temp = 2 + MetaData->TopGridRank;  // 2 fields we're tracking plus TopGridRank position fields

  if(TracerParticleOutputVelocity)
    int_temp += MetaData->TopGridRank;  // TopGridRank velocity fields


  // printf("tracer particle: writing %d fields and sizeof(int) is %d\n",int_temp,sizeof(int));

  fwrite((void*) &int_temp, sizeof(int), 1, fptr);  // write number of fields
 
  /* --------------------------------------------------------------- */
  /* Loop over grids and write grid data to files. */
 
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    /* Loop over all the grids. */
 
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
 
      /* Write out grid info (also deletes the under subgrid field). */
 
      if (Temp->GridData->TracerParticleOutputData(fptr, WriteTime) == FAIL) {
	ENZO_FAIL("Error in grid->OutputTracerParticleData.\n");

      }
 
      /* Next grid on this level. */
 
      Temp = Temp->NextGridThisLevel;
 
    } // end loop over grids
 
  } // end loop over levels
 
  /* Write a 0 to indicate no more data at this time. */
 
  fwrite((void*) &Zero, sizeof(int), 1, fptr);
 
  /* close file. */
 
  fclose(fptr);
 
  return SUCCESS;
}
