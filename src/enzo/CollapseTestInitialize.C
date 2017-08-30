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

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, double *MassUnits, FLOAT Time);

static float CollapseTestInitialFractionHII   = 1.2e-5;
static float CollapseTestInitialFractionHeII  = 1.0e-14;
static float CollapseTestInitialFractionHeIII = 1.0e-17;
static float CollapseTestInitialFractionHM    = 2.0e-9;
static float CollapseTestInitialFractionH2I   = 2.0e-20;
static float CollapseTestInitialFractionH2II  = 3.0e-14;

int CollapseTestInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior)
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

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int CollapseTestNumberOfSpheres = 1;
  int CollapseTestRefineAtStart   = TRUE;
  int CollapseTestUseParticles    = FALSE;
  float CollapseTestParticleMeanDensity = FLOAT_UNDEFINED;
  int CollapseTestUseColour       = FALSE;
  int CollapseTestUseMetals       = FALSE;
  int CollapseTestWind            = FALSE;
  float CollapseTestInitialTemperature = 1000;
  float CollapseTestInitialDensity     = 1.0;
  float CollapseTestSphereDensity[MAX_SPHERES],
    CollapseTestSphereTemperature[MAX_SPHERES],
    CollapseTestSphereVelocity[MAX_SPHERES][MAX_DIMENSION],
    CollapseTestUniformVelocity[MAX_DIMENSION],
    CollapseTestFracKeplerianRot[MAX_SPHERES],
    CollapseTestSphereTurbulence[MAX_SPHERES],
    CollapseTestSphereDispersion[MAX_SPHERES],
    CollapseTestSphereCutOff[MAX_SPHERES],
    CollapseTestSphereAng1[MAX_SPHERES],
    CollapseTestSphereAng2[MAX_SPHERES],
    CollapseTestSphereMetallicity[MAX_SPHERES],
    CollapseTestSphereSmoothRadius[MAX_SPHERES],
    CollapseTestSphereHIIFraction[MAX_SPHERES],
    CollapseTestSphereHeIIFraction[MAX_SPHERES],
    CollapseTestSphereHeIIIFraction[MAX_SPHERES],
    CollapseTestSphereH2IFraction[MAX_SPHERES],
    CollapseTestWindVelocity[MAX_DIMENSION];
  int CollapseTestSphereNumShells[MAX_SPHERES],
    CollapseTestSphereInitialLevel[MAX_SPHERES],
    CollapseTestSphereType[MAX_SPHERES],
    CollapseTestSphereConstantPressure[MAX_SPHERES],
    CollapseTestSphereSmoothSurface[MAX_SPHERES];
  FLOAT CollapseTestSphereRadius[MAX_SPHERES],
    CollapseTestSphereCoreRadius[MAX_SPHERES],
    CollapseTestSpherePosition[MAX_SPHERES][MAX_DIMENSION];

  for (sphere = 0; sphere < MAX_SPHERES; sphere++) {
    CollapseTestSphereRadius[sphere]     = 1.0;
    CollapseTestSphereCoreRadius[sphere] = 0.1;
    CollapseTestSphereDensity[sphere]    = 1.0;
    CollapseTestSphereTemperature[sphere] = 1.0;
    CollapseTestFracKeplerianRot[sphere] = 0.0;
    CollapseTestSphereTurbulence[sphere] = 0.0;
    CollapseTestSphereDispersion[sphere] = 0.0;
    CollapseTestSphereCutOff[sphere] = 6.5;
    CollapseTestSphereAng1[sphere] = 0;
    CollapseTestSphereAng2[sphere] = 0;
    CollapseTestSphereNumShells[sphere] = 1;
    CollapseTestSphereSmoothRadius[sphere] = 1.2;
    CollapseTestSphereMetallicity[sphere] = tiny_number;
    CollapseTestSphereInitialLevel[sphere] = 0;
    CollapseTestSphereHIIFraction[sphere] = CollapseTestInitialFractionHII;
    CollapseTestSphereHeIIFraction[sphere] = CollapseTestInitialFractionHeII;
    CollapseTestSphereHeIIIFraction[sphere] = CollapseTestInitialFractionHeIII;
    CollapseTestSphereH2IFraction[sphere] = CollapseTestInitialFractionH2I;

    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CollapseTestSpherePosition[sphere][dim] = 0.5*(DomainLeftEdge[dim] +
						     DomainRightEdge[dim]);
      CollapseTestSphereVelocity[sphere][dim] = 0;
    }
    CollapseTestSphereType[sphere]       = 0;
    CollapseTestSphereConstantPressure[sphere] = FALSE;
    CollapseTestSphereSmoothSurface[sphere] = FALSE;
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    CollapseTestUniformVelocity[dim] = 0;
    CollapseTestWindVelocity[dim] = 0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "CollapseTestNumberOfSpheres = %"ISYM,
		  &CollapseTestNumberOfSpheres);
    ret += sscanf(line, "CollapseTestRefineAtStart = %"ISYM, 
		  &CollapseTestRefineAtStart);
    ret += sscanf(line, "CollapseTestUseParticles = %"ISYM, 
		  &CollapseTestUseParticles);
    ret += sscanf(line, "CollapseTestParticleMeanDensity = %"FSYM,
		  &CollapseTestParticleMeanDensity);
    ret += sscanf(line, "CollapseTestUseColour = %"ISYM, 
		  &CollapseTestUseColour);
    ret += sscanf(line, "CollapseTestUseMetals = %"ISYM, 
		  &CollapseTestUseMetals);
    ret += sscanf(line, "CollapseTestWind = %"ISYM, 
                  &CollapseTestWind);
    ret += sscanf(line, "CollapseTestInitialTemperature = %"FSYM, 
		  &CollapseTestInitialTemperature);
    ret += sscanf(line, "CollapseTestInitialDensity = %"FSYM,
		  &CollapseTestInitialDensity);
    ret += sscanf(line, "CollapseTestUniformVelocity = %"FSYM" %"FSYM" %"FSYM, 
		  CollapseTestUniformVelocity, CollapseTestUniformVelocity+1,
		  CollapseTestUniformVelocity+2);
    if (sscanf(line, "CollapseTestSphereType[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereType[%"ISYM"] = %"ISYM, &sphere,
		    &CollapseTestSphereType[sphere]);
    if (sscanf(line, "CollapseTestSphereRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereRadius[%"ISYM"] = %"PSYM, &sphere,
		    &CollapseTestSphereRadius[sphere]);
    if (sscanf(line, "CollapseTestSphereCoreRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereCoreRadius[%"ISYM"] = %"PSYM, &sphere,
		    &CollapseTestSphereCoreRadius[sphere]);
    if (sscanf(line, "CollapseTestSphereDensity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereDensity[%"ISYM"] = %"FSYM, &sphere,
		    &CollapseTestSphereDensity[sphere]);
    if (sscanf(line, "CollapseTestSphereTemperature[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereTemperature[%"ISYM"] = %"FSYM, &sphere,
		    &CollapseTestSphereTemperature[sphere]);
    if (sscanf(line, "CollapseTestSphereMetallicity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereMetallicity[%"ISYM"] = %"FSYM, &sphere,
		    &CollapseTestSphereMetallicity[sphere]);
    if (sscanf(line, "CollapseTestSpherePosition[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSpherePosition[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM, 
		    &sphere, &CollapseTestSpherePosition[sphere][0],
		    &CollapseTestSpherePosition[sphere][1],
		    &CollapseTestSpherePosition[sphere][2]);
    if (sscanf(line, "CollapseTestSphereVelocity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereVelocity[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
		    &sphere, &CollapseTestSphereVelocity[sphere][0],
		    &CollapseTestSphereVelocity[sphere][1],
		    &CollapseTestSphereVelocity[sphere][2]);
    if (sscanf(line, "CollapseTestFracKeplerianRot[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestFracKeplerianRot[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestFracKeplerianRot[sphere]);
    if (sscanf(line, "CollapseTestSphereTurbulence[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereTurbulence[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestSphereTurbulence[sphere]);
    if (sscanf(line, "CollapseTestSphereDispersion[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereDispersion[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestSphereDispersion[sphere]);
    if (sscanf(line, "CollapseTestSphereCutOff[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereCutOff[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestSphereCutOff[sphere]);
    if (sscanf(line, "CollapseTestSphereAng1[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereAng1[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestSphereAng1[sphere]);
    if (sscanf(line, "CollapseTestSphereAng2[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereAng2[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestSphereAng2[sphere]);
    if (sscanf(line, "CollapseTestSphereNumShells[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereNumShells[%"ISYM"] = %"ISYM, &sphere,
                    &CollapseTestSphereNumShells[sphere]);
    if (sscanf(line, "CollapseTestSphereInitialLevel[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereInitialLevel[%"ISYM"] = %"ISYM, &sphere,
                    &CollapseTestSphereInitialLevel[sphere]);
    if (sscanf(line, "CollapseTestSphereConstantPressure[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereConstantPressure[%"ISYM"] = %"ISYM, &sphere,
		    &CollapseTestSphereConstantPressure[sphere]);
    if (sscanf(line, "CollapseTestSphereSmoothSurface[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereSmoothSurface[%"ISYM"] = %"ISYM, &sphere,
		    &CollapseTestSphereSmoothSurface[sphere]);
    if (sscanf(line, "CollapseTestSphereSmoothRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereSmoothRadius[%"FSYM"] = %"FSYM, &sphere,
		    &CollapseTestSphereSmoothRadius[sphere]);
    if (sscanf(line, "CollapseTestSphereHIIFraction[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereHIIFraction[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestSphereHIIFraction[sphere]);
    if (sscanf(line, "CollapseTestSphereHeIIFraction[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereHeIIFraction[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestSphereHeIIFraction[sphere]);
    if (sscanf(line, "CollapseTestSphereHeIIIFraction[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereHeIIIFraction[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestSphereHeIIIFraction[sphere]);
    if (sscanf(line, "CollapseTestSphereH2IFraction[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "CollapseTestSphereH2IFraction[%"ISYM"] = %"FSYM, &sphere,
                    &CollapseTestSphereH2IFraction[sphere]);

    ret += sscanf(line, "CollapseTestInitialFractionHII = %"FSYM,
		  &CollapseTestInitialFractionHII);
    ret += sscanf(line, "CollapseTestInitialFractionHeII = %"FSYM,
		  &CollapseTestInitialFractionHeII);
    ret += sscanf(line, "CollapseTestInitialFractionHeIII = %"FSYM,
		  &CollapseTestInitialFractionHeIII);
    ret += sscanf(line, "CollapseTestInitialFractionHM = %"FSYM,
		  &CollapseTestInitialFractionHM);
    ret += sscanf(line, "CollapseTestInitialFractionH2I = %"FSYM,
		  &CollapseTestInitialFractionH2I);
    ret += sscanf(line, "CollapseTestInitialFractionH2II = %"FSYM,
		  &CollapseTestInitialFractionH2II);

    ret += sscanf(line, "CollapseTestWindVelocity = %"FSYM" %"FSYM" %"FSYM, 
                  &CollapseTestWindVelocity[0],&CollapseTestWindVelocity[1],&CollapseTestWindVelocity[2]);


    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "CollapseTest") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up grid */

  if (TopGrid.GridData->CollapseTestInitializeGrid(
	     CollapseTestNumberOfSpheres, CollapseTestSphereRadius,
	     CollapseTestSphereCoreRadius, CollapseTestSphereDensity,
	     CollapseTestSphereTemperature, CollapseTestSphereMetallicity,
	     CollapseTestSpherePosition, CollapseTestSphereVelocity,
             CollapseTestFracKeplerianRot, CollapseTestSphereTurbulence,
	     CollapseTestSphereDispersion,
             CollapseTestSphereCutOff, CollapseTestSphereAng1,
             CollapseTestSphereAng2, CollapseTestSphereNumShells,
	     CollapseTestSphereType, CollapseTestSphereConstantPressure,
	     CollapseTestSphereSmoothSurface, CollapseTestSphereSmoothRadius, 
	     CollapseTestSphereHIIFraction, CollapseTestSphereHeIIFraction,
	     CollapseTestSphereHeIIIFraction, CollapseTestSphereH2IFraction,
	     CollapseTestUseParticles, CollapseTestParticleMeanDensity,
             CollapseTestUniformVelocity, CollapseTestUseColour,
	     CollapseTestUseMetals, 
             CollapseTestInitialTemperature, CollapseTestInitialDensity,
	     0,
	     CollapseTestInitialFractionHII, CollapseTestInitialFractionHeII,
	     CollapseTestInitialFractionHeIII, CollapseTestInitialFractionHM,
	     CollapseTestInitialFractionH2I, CollapseTestInitialFractionH2II) == FAIL) {
    ENZO_FAIL("Error in CollapseTestInitializeGrid.");
  }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */
  for (int count = 0; count < MAX_FLAGGING_METHODS; count++)
    if (MinimumMassForRefinement[count] == FLOAT_UNDEFINED) {
      MinimumMassForRefinement[count] = MinimumOverDensityForRefinement[count];
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
        MinimumMassForRefinement[count] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	     float(MetaData.TopGridDims[dim]);
    }


  /* If requested and there are no manual settings of the refinement
     of spheres, refine the grid to the desired level. */

  int MaxInitialLevel = 0;
  for (sphere = 0; sphere < CollapseTestNumberOfSpheres; sphere++)
    MaxInitialLevel = max(MaxInitialLevel, CollapseTestSphereInitialLevel[sphere]);

  if (CollapseTestRefineAtStart) {

    /* If the user specified an initial refinement level for a sphere,
       then manually create the hierarchy first. */

    if (MaxInitialLevel > 0) {

      int lev, max_level;
      float dx;
      HierarchyEntry **Subgrid;
      int NumberOfSubgridDims[MAX_DIMENSION];
      FLOAT ThisLeftEdge[MAX_DIMENSION], ThisRightEdge[MAX_DIMENSION];

      for (sphere = 0; sphere < CollapseTestNumberOfSpheres; sphere++) {
	
	max_level = CollapseTestSphereInitialLevel[sphere];
	if (max_level > 0) {

	  Subgrid = new HierarchyEntry*[max_level];
	  for (lev = 0; lev < max_level; lev++)
	    Subgrid[lev] = new HierarchyEntry;

	  for (lev = 0; lev < max_level; lev++) {
	    
	    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	      dx = 1.0 / float(MetaData.TopGridDims[dim]) / POW(RefineBy, lev);
	      ThisLeftEdge[dim] = CollapseTestSpherePosition[sphere][dim] -
		0.5 * CollapseTestSphereRadius[sphere] - 2*dx;  // plus some buffer
	      ThisLeftEdge[dim] = nint(ThisLeftEdge[dim] / dx) * dx;
	      ThisRightEdge[dim] = CollapseTestSpherePosition[sphere][dim] +
		0.5 * CollapseTestSphereRadius[sphere] + 2*dx;
	      ThisRightEdge[dim] = nint(ThisRightEdge[dim] / dx) * dx;
	      NumberOfSubgridDims[dim] = 
		nint((ThisRightEdge[dim] - ThisLeftEdge[dim]) / 
		     (DomainRightEdge[dim] - DomainLeftEdge[dim]) / dx);		
	    } // ENDFOR dims

	    if (debug)
	      printf("CollapseTest:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n",
		     lev+1, NumberOfSubgridDims[0]);
	    
	    if (NumberOfSubgridDims[0] > 0) {

	      // Insert into AMR hierarchy
	      if (lev == 0) {
		Subgrid[lev]->NextGridThisLevel = TopGrid.NextGridNextLevel;
		TopGrid.NextGridNextLevel = Subgrid[lev];
		Subgrid[lev]->ParentGrid = &TopGrid;
	      } else {
		Subgrid[lev]->NextGridThisLevel = NULL;
		Subgrid[lev]->ParentGrid = Subgrid[lev-1];
	      }
	      if (lev == max_level-1)
		Subgrid[lev]->NextGridNextLevel = NULL;
	      else
		Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];

	      // Create grid
	      for (dim = 0; dim < MetaData.TopGridRank; dim++)
		NumberOfSubgridDims[dim] += 2*NumberOfGhostZones;
	      Subgrid[lev]->GridData = new grid;
	      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
	      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, 
						  NumberOfSubgridDims,
						  ThisLeftEdge,
						  ThisRightEdge, 0);


	      if (Subgrid[lev]->GridData->CollapseTestInitializeGrid(
	          CollapseTestNumberOfSpheres, CollapseTestSphereRadius,
		  CollapseTestSphereCoreRadius, CollapseTestSphereDensity,
		  CollapseTestSphereTemperature, CollapseTestSphereMetallicity,
		  CollapseTestSpherePosition, CollapseTestSphereVelocity,
		  CollapseTestFracKeplerianRot, CollapseTestSphereTurbulence,
		  CollapseTestSphereDispersion,
		  CollapseTestSphereCutOff, CollapseTestSphereAng1,
		  CollapseTestSphereAng2, CollapseTestSphereNumShells,
		  CollapseTestSphereType, CollapseTestSphereConstantPressure,
		  CollapseTestSphereSmoothSurface, CollapseTestSphereSmoothRadius, 
		  CollapseTestSphereHIIFraction, CollapseTestSphereHeIIFraction,
		  CollapseTestSphereHeIIIFraction, CollapseTestSphereH2IFraction,
		  CollapseTestUseParticles, CollapseTestParticleMeanDensity,
		  CollapseTestUniformVelocity, CollapseTestUseColour,
		  CollapseTestUseMetals,
		  CollapseTestInitialTemperature, CollapseTestInitialDensity,
		  lev-1,
		  CollapseTestInitialFractionHII, CollapseTestInitialFractionHeII,
		  CollapseTestInitialFractionHeIII, CollapseTestInitialFractionHM,
		  CollapseTestInitialFractionH2I, CollapseTestInitialFractionH2II) == FAIL) {
		ENZO_FAIL("Error in CollapseTestInitializeGrid.");
	      }
	      
	    } // ENDIF zones exist
	  } // ENDFOR levels
	} // ENDIF max_level > 0
      } // ENDFOR spheres
    } // ENDIF MaxInitialLevel > 0

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    if (MaxInitialLevel == 0) {
      for (level = 0; level < MaximumRefinementLevel; level++) {
	if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	  ENZO_FAIL("Error in RebuildHierarchy.");
	}
	if (LevelArray[level+1] == NULL)
	  break;
	LevelHierarchyEntry *Temp = LevelArray[level+1];
	while (Temp != NULL) {
	  if (Temp->GridData->CollapseTestInitializeGrid(
		 CollapseTestNumberOfSpheres, CollapseTestSphereRadius,
		 CollapseTestSphereCoreRadius, CollapseTestSphereDensity,
		 CollapseTestSphereTemperature, CollapseTestSphereMetallicity,
		 CollapseTestSpherePosition, CollapseTestSphereVelocity,
		 CollapseTestFracKeplerianRot, CollapseTestSphereTurbulence,
		 CollapseTestSphereDispersion,
		 CollapseTestSphereCutOff, CollapseTestSphereAng1,
		 CollapseTestSphereAng2, CollapseTestSphereNumShells,
		 CollapseTestSphereType, CollapseTestSphereConstantPressure,
		 CollapseTestSphereSmoothSurface, CollapseTestSphereSmoothRadius, 
		 CollapseTestSphereHIIFraction, CollapseTestSphereHeIIFraction,
		 CollapseTestSphereHeIIIFraction, CollapseTestSphereH2IFraction,
		 CollapseTestUseParticles, CollapseTestParticleMeanDensity,
		 CollapseTestUniformVelocity, CollapseTestUseColour,
		 CollapseTestUseMetals,
		 CollapseTestInitialTemperature, CollapseTestInitialDensity,
		 level+1,
		 CollapseTestInitialFractionHII, CollapseTestInitialFractionHeII,
		 CollapseTestInitialFractionHeIII, CollapseTestInitialFractionHM,
		 CollapseTestInitialFractionH2I, CollapseTestInitialFractionH2II) == FAIL) {
	    ENZO_FAIL("Error in CollapseTestInitializeGrid.");
	  }
	  Temp = Temp->NextGridThisLevel;
	}
      } // end: loop over levels
    } // ENDELSE manually set refinement levels

      /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
			      *LevelArray[level-1]->GridData) == FAIL) {
	  ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.");
	}
	Temp = Temp->NextGridThisLevel;
      }
    }

  } // end: if (CollapseTestRefineAtStart)

  /* If there is wind, initialize the exterior */
 
  if (CollapseTestWind) {
    Exterior.Prepare(TopGrid.GridData);

    const int MAX_BNDRY_VARS = 6;
    double mh = 1.67e-24;
    float mu=0.59;
    double kboltz = 1.38e-16;
    float InflowValue[MAX_BNDRY_VARS], Dummy[MAX_BNDRY_VARS];
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
    double MassUnits;
    if (GetUnits(&DensityUnits, &LengthUnits,&TemperatureUnits, &TimeUnits,
               &VelocityUnits, &MassUnits, MetaData.Time) == FAIL){
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
    }
    float EnergyUnits;
    float TempToEnergyConversion;
    EnergyUnits = POW(LengthUnits, 2.0) / POW(TimeUnits, 2.0);
    TempToEnergyConversion =  kboltz/((Gamma - 1.0)*mu*mh);
    TempToEnergyConversion /= EnergyUnits;  // this times temperature gives you energy units in ENZO UNITS (K -> Enzo)

    InflowValue[0] = CollapseTestInitialDensity;
    InflowValue[1] = CollapseTestInitialTemperature*TempToEnergyConversion;
    if (HydroMethod != 2) {
      InflowValue[1] = InflowValue[1] + 0.5*(POW(CollapseTestWindVelocity[0]/VelocityUnits,2)
                                                    + POW(CollapseTestWindVelocity[1]/VelocityUnits,2)
                                                    + POW(CollapseTestWindVelocity[2]/VelocityUnits,2));
    }
    InflowValue[2] = CollapseTestWindVelocity[0]/VelocityUnits;
    InflowValue[3] = CollapseTestWindVelocity[1]/VelocityUnits;
    InflowValue[4] = CollapseTestWindVelocity[2]/VelocityUnits;
    if (CollapseTestUseMetals)
      InflowValue[5] = 1.0e-10; ///need to be changed

    if (Exterior.InitializeExternalBoundaryFace(0, inflow, outflow, InflowValue,
                                                Dummy) == FAIL) {
      fprintf(stderr, "Error in InitializeExternalBoundaryFace.\n");
      return FAIL;
    }
    if (MetaData.TopGridRank > 1)
      Exterior.InitializeExternalBoundaryFace(1, periodic, periodic,
                                              Dummy, Dummy);
    if (MetaData.TopGridRank > 2)
      Exterior.InitializeExternalBoundaryFace(2, periodic, periodic,
                                              Dummy, Dummy);
  }




  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = (char*) GEName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;
  if (MultiSpecies) {
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
  }  // if Multispecies
  if (CollapseTestUseColour)
    DataLabel[count++] = (char*) ColourName;
  if (CollapseTestUseMetals)
    DataLabel[count++] = (char*) MetalName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
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
    fprintf(Outfptr, "CollapseTestWind               = $"ISYM"\n",
            CollapseTestWind);
    fprintf(Outfptr, "CollapseTestInitialTemperature = %"FSYM"\n",
	    CollapseTestInitialTemperature);
    fprintf(Outfptr, "CollapseTestInitialDensity     = %"FSYM"\n",
	    CollapseTestInitialDensity);
    fprintf(Outfptr, "CollapseTestUniformVelocity    = %"FSYM" %"FSYM" %"FSYM"\n",
	    CollapseTestUniformVelocity[0], CollapseTestUniformVelocity[1],
	    CollapseTestUniformVelocity[2]);
    fprintf(Outfptr, "CollapseTestWindVelocity    = %"FSYM" %"FSYM" %"FSYM"\n",
            CollapseTestWindVelocity[0], CollapseTestWindVelocity[1],
            CollapseTestWindVelocity[2]);
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
  }

  return SUCCESS;

}
