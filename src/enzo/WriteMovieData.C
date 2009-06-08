/***********************************************************************
/
/  WRITE OUT THE SO-CALLED MOVIE DATA
/
/  written by: Greg Bryan
/  date:       March, 2000
/  modified1:
/
/  PURPOSE: This function walks through the grids and writes out selected
/     information (positions and densities) of grid and particle data
/     within a selected region.  Each processor creates its own file.
/
************************************************************************/
 
 
#include <string.h>
#include <stdio.h>


#include <math.h>
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
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
 
 
static char DMSuffix[]      = ".dm";
static char GridSuffix[]    = ".grid";
static char StarSuffix[]    = ".star";
static char SummarySuffix[] = ".summary";
 
int WriteMovieData(char *basename, int filenumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime)
{
 
  /* declarations */
 
  char id[7], name[MAX_LINE_LENGTH], gridname[MAX_LINE_LENGTH],
       dmname[MAX_LINE_LENGTH], starname[MAX_LINE_LENGTH],
       summaryname[MAX_LINE_LENGTH];
  FILE *Gridfptr = NULL, *DMfptr = NULL, *Starfptr = NULL, *Summaryfptr;
  int level, type, value;
 
  /* Create main name */
 
  strcpy(name, basename);
  sprintf(id, "%6.6"ISYM, filenumber);   /* create output # */
  strcat(name, id);
  sprintf(id, "%-d", MyProcessorNumber);  /* create processor # */
  //  if (NumberOfProcessors == 1)
  //    strcpy(id, "");
 
  if (debug)
    printf("WriteMovieData: writing file %s.\n", name);
 
  /* Set the Cell width of the root grid. */
 
  float BaseRadius = (DomainRightEdge[0] - DomainLeftEdge[0])/
      FLOAT(MetaData->TopGridDims[0]);
 
  /* Set data file name. */
 
  strcpy(gridname, name);
  strcat(gridname, GridSuffix);
  strcat(gridname, id);
 
  strcpy(dmname, name);
  strcat(dmname, DMSuffix);
  strcat(dmname, id);
 
  strcpy(starname, name);
  strcat(starname, StarSuffix);
  strcat(starname, id);
 
  strcpy(summaryname, name);
  strcat(summaryname, SummarySuffix);
  strcat(summaryname, id);
 
  /* Open files. */
 
  if ((Gridfptr = fopen(gridname, "w")) == NULL) {
    fprintf(stderr, "Error opening grid output file %s\n", gridname);
    return FAIL;
  }
 
  if ((DMfptr = fopen(dmname, "w")) == NULL) {
    fprintf(stderr, "Error opening dm output file %s\n", dmname);
    return FAIL;
  }
 
  if (StarParticleCreation > 0)
    if ((Starfptr = fopen(starname, "w")) == NULL) {
      fprintf(stderr, "Error opening star output file %s\n", starname);
      return FAIL;
    }
 
  /* --------------------------------------------------------------- */
  /* Loop over grids and write grid data to files. */
 
  int NumberOfPoints[3] ={0,0,0}, NumberOfValuesPerPoint[3] = {0,0,0};
  char *PointValueNames[3][20];  /* should make 20 into parameter. */
  for (type = 0; type < 3; type++)
    for (value = 0; value < 20; value++)
      PointValueNames[type][value] = NULL;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    /* Loop over all the grids. */
 
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
 
      /* Initialize the UNDER_SUBGRID_FIELD for this grid. */
 
      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
 
      /* Zero the solution (on this grid) which is underneath any subgrid
	 (so we get only the high resolution solution from the subgrid). */
 
      LevelHierarchyEntry *Temp2 = LevelArray[level+1];
      while (Temp2 != NULL) {
	Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData,
						 ZERO_UNDER_SUBGRID_FIELD);
	Temp2 = Temp2->NextGridThisLevel;
      }
 
      /* Write out grid info (also deletes the under subgrid field). */
 
      if (Temp->GridData->OutputGridMovieData(Gridfptr, DMfptr, Starfptr,
				    MetaData->MovieRegionLeftEdge,
				    MetaData->MovieRegionRightEdge,
				    WriteTime, NumberOfPoints,
				    NumberOfValuesPerPoint,
				    PointValueNames, BaseRadius) == FAIL) {
	fprintf(stderr, "Error in grid->OutputGridMovieData.\n");
	return FAIL;
      }
 
      /* Next grid on this level. */
 
      Temp = Temp->NextGridThisLevel;
 
    } // end loop over grids
 
  } // end loop over levels
 
  /* close files. */
 
  if (Gridfptr != NULL)
    fclose(Gridfptr);
  if (DMfptr != NULL)
    fclose(DMfptr);
  if (Starfptr != NULL)
    fclose(Starfptr);
 
  /* Compute redshift. */
 
  FLOAT a = 1, dadt, Redshift = 0;
  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1, MassUnits=1;
  double Mpc = 3.086e24, SolarMass = 1.989e33;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, WriteTime) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;    
  }
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(WriteTime, &a, &dadt);
    Redshift = (1.0+InitialRedshift)/a - 1;
  }
 
  if ((Summaryfptr = fopen(summaryname, "w")) == NULL) {
    fprintf(stderr, "Error opening summary output file %s\n", summaryname);
    return FAIL;
  }
 
  /* Output summary information:
     0-2 - Number of grid point/dm/star points (3 integers)
     3-5 - Number of values for grid/dm/star points (3 integers)
     6   - time (code units)  (1 float)
     7   - redshift. (1 float)
     8   - dimension of problem (1 integer)
     9   - number of processors (1 integer)
     10-12 - left corner of region used for extraction (3 floats)
     13-15 - right corner of region used for extraction (3 floats) */
 
  fprintf(Summaryfptr, "NumberOfParticles          = ");
  WriteListOfInts(Summaryfptr, 3, NumberOfPoints);
  fprintf(Summaryfptr, "NumberOfValuesPerParticle  = ");
  WriteListOfInts(Summaryfptr, 3, NumberOfValuesPerPoint);
  fprintf(Summaryfptr, "Time                       = %"GOUTSYM"\n", WriteTime);
  fprintf(Summaryfptr, "Redshift                   = %"GOUTSYM"\n", Redshift);
  fprintf(Summaryfptr, "Rank                       = %"ISYM"\n",
	  MetaData->TopGridRank);
  fprintf(Summaryfptr, "NumberOfProcessors         = %"ISYM"\n",
	  NumberOfProcessors);
  fprintf(Summaryfptr, "RegionLeftEdge             = ");
  WriteListOfFloats(Summaryfptr, 3, MetaData->MovieRegionLeftEdge);
  fprintf(Summaryfptr, "RegionRightEdge            = ");
  WriteListOfFloats(Summaryfptr, 3, MetaData->MovieRegionRightEdge);
  fprintf(Summaryfptr, "#LengthCGSConversionFactor     = %"GSYM"\n", LengthUnits);
  fprintf(Summaryfptr, "#MassCGSConversionFactor       = %"GSYM"\n",
	  double(DensityUnits)*POW(double(LengthUnits), 3));
  fprintf(Summaryfptr, "#LengthMpcConversionFactor     = %"GSYM"\n",
	  LengthUnits/Mpc);
  fprintf(Summaryfptr, "#MassSolarMassConversionFactor = %"GSYM"\n",
	  double(DensityUnits)*POW(double(LengthUnits), 3)/SolarMass);
  fprintf(Summaryfptr, "#VelocityCGSConversionFactor   = %"GSYM"\n", VelocityUnits);
  fprintf(Summaryfptr, "#VelocityPositionFactor        = %"GOUTSYM"\n", a);
 
  for (type = 0; type < 3; type++)
    if (NumberOfValuesPerPoint[type] > 0) {
      fprintf(Summaryfptr, "\n#Quantities for particle type %"ISYM" ", type);
      if (type == 0) fprintf(Summaryfptr, "(grid):\n");
      if (type == 1) fprintf(Summaryfptr, "(dark-matter):\n");
      if (type == 2) fprintf(Summaryfptr, "(star):\n");
      for (value = 0; value < NumberOfValuesPerPoint[type]; value++)
	fprintf(Summaryfptr, "# %s\n", PointValueNames[type][value]);
    }
 
  fclose(Summaryfptr);
 
  return SUCCESS;
}
