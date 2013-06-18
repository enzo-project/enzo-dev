/***********************************************************************
/
/  INITIALIZE PHOTON TEST
/
/  written by: Tom Abel
/  date:       Oct 2003-
/  modified1:
/
/  PURPOSE:
/    Set up a number of points sources. Modeled after GB's CollapseTest
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
#include "CosmologyParameters.h"


void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS);
int ReadPhotonSources(FILE *fptr, FLOAT CurrentTime);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

static float PhotonTestInitialFractionHII   = 1.2e-5;
static float PhotonTestInitialFractionHeII  = 1.0e-14;
static float PhotonTestInitialFractionHeIII = 1.0e-17;
static float PhotonTestInitialFractionHM    = 2.0e-9;
static float PhotonTestInitialFractionH2I   = 2.0e-20;
static float PhotonTestInitialFractionH2II  = 3.0e-14;

int PhotonTestInitialize(FILE *fptr, FILE *Outfptr, 
			 HierarchyEntry &TopGrid, TopGridData &MetaData,
			 bool Reinitialize)
{
  const char *DensName = "Density";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *ColourName = "colour";
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
  const char *kdissH2IName = "H2I_kdiss"; 
  const char *RadAccel1Name = "x-RadPressure";
  const char *RadAccel2Name = "y-RadPressure";
  const char *RadAccel3Name = "z-RadPressure";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  char *dummy = new char[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i, source;
  int   TotalRefinement;

  /* set default parameters */

  char *PhotonTestDensityFilename = NULL;
  int PhotonTestNumberOfSpheres = 1;
  int PhotonTestUseParticles    = FALSE;
  int PhotonTestUseColour       = FALSE;
  float PhotonTestInitialTemperature = 1000;
  int   PhotonTestSphereType[MAX_SPHERES],
    PhotonTestSphereConstantPressure[MAX_SPHERES],
    PhotonTestSphereSmoothSurface[MAX_SPHERES];
  float PhotonTestSphereDensity[MAX_SPHERES],
    PhotonTestSphereTemperature[MAX_SPHERES],
    PhotonTestSphereVelocity[MAX_SPHERES][MAX_DIMENSION],
    PhotonTestUniformVelocity[MAX_DIMENSION],
    PhotonTestFracKeplerianRot[MAX_SPHERES],
    PhotonTestSphereTurbulence[MAX_SPHERES],
    PhotonTestSphereCutOff[MAX_SPHERES],
    PhotonTestSphereAng1[MAX_SPHERES],
    PhotonTestSphereAng2[MAX_SPHERES],
    PhotonTestSphereSmoothRadius[MAX_SPHERES],
    PhotonTestSphereRadius[MAX_SPHERES],
    PhotonTestSphereCoreRadius[MAX_SPHERES],
    PhotonTestSphereHIIFraction[MAX_SPHERES],
    PhotonTestSphereHeIIFraction[MAX_SPHERES],
    PhotonTestSphereHeIIIFraction[MAX_SPHERES],
    PhotonTestSphereH2IFraction[MAX_SPHERES];
  int PhotonTestSphereNumShells[MAX_SPHERES];
  FLOAT PhotonTestSpherePosition[MAX_SPHERES][MAX_DIMENSION];

  float PhotonTestOmegaBaryonNow=0.05;
  int   PhotonTestRefineAtStart = 0;

  rewind(fptr);

  // Set default values
  if (debug)
    if (Reinitialize)
      fprintf(stderr, "PhotonTestInitialize: Reinitializing after root "
	      "grid split.\n");
    else
      fprintf(stderr, "PhotonTestInitialize: Set up test problem.\n");

  for (sphere = 0; sphere < MAX_SPHERES; sphere++) {
    PhotonTestSphereRadius[sphere]     = 0.5;
    PhotonTestSphereCoreRadius[sphere] = 0.1;
    PhotonTestSphereDensity[sphere]    = 1.0;
    PhotonTestSphereTemperature[sphere] = 1.0;
    PhotonTestFracKeplerianRot[sphere] = 0.0;
    PhotonTestSphereTurbulence[sphere] = 0.0;
    PhotonTestSphereCutOff[sphere] = 6.5;
    PhotonTestSphereAng1[sphere] = 0;
    PhotonTestSphereAng2[sphere] = 0;
    PhotonTestSphereNumShells[sphere] = 1;
    PhotonTestSphereSmoothRadius[sphere] = 1.2;
    PhotonTestSphereHIIFraction[sphere] = PhotonTestInitialFractionHII;
    PhotonTestSphereHeIIFraction[sphere] = PhotonTestInitialFractionHeII;
    PhotonTestSphereHeIIIFraction[sphere] = PhotonTestInitialFractionHeIII;
    PhotonTestSphereH2IFraction[sphere] = PhotonTestInitialFractionH2I;

    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      PhotonTestSpherePosition[sphere][dim] = 0.5*(DomainLeftEdge[dim] +
						     DomainRightEdge[dim]);
      PhotonTestSphereVelocity[sphere][dim] = 0;
    }
    PhotonTestSphereType[sphere]       = 0;
    PhotonTestSphereConstantPressure[sphere] = FALSE;
    PhotonTestSphereSmoothSurface[sphere] = FALSE;
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    PhotonTestUniformVelocity[dim] = 0;


  /* read input from file */
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;
    /* read parameters */

    ret += sscanf(line, "PhotonTestNumberOfSpheres = %"ISYM,
		  &PhotonTestNumberOfSpheres);
    ret += sscanf(line, "PhotonTestRefineAtStart = %"ISYM, 
		  &PhotonTestRefineAtStart);
    ret += sscanf(line, "PhotonTestUseParticles = %"ISYM, 
		  &PhotonTestUseParticles);
    ret += sscanf(line, "PhotonTestUseColour = %"ISYM, 
		  &PhotonTestUseColour);
    ret += sscanf(line, "PhotonTestInitialTemperature = %"FSYM, 
		  &PhotonTestInitialTemperature);
    if (sscanf(line, "PhotonTestDensityFilename = %s", dummy) == 1) {
      ret++;
      PhotonTestDensityFilename = dummy;
    }
    ret += sscanf(line, "PhotonTestUniformVelocity = %"FSYM" %"FSYM" %"FSYM, 
		  PhotonTestUniformVelocity, PhotonTestUniformVelocity+1,
		  PhotonTestUniformVelocity+2);
    if (sscanf(line, "PhotonTestSphereType[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereType[%"ISYM"] = %"ISYM, &sphere,
		    &PhotonTestSphereType[sphere]);
    if (sscanf(line, "PhotonTestSphereConstantPressure[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereConstantPressure[%"ISYM"] = %"ISYM, &sphere,
		    &PhotonTestSphereConstantPressure[sphere]);
    if (sscanf(line, "PhotonTestSphereSmoothSurface[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereSmoothSurface[%"ISYM"] = %"ISYM, &sphere,
		    &PhotonTestSphereSmoothSurface[sphere]);
    if (sscanf(line, "PhotonTestSphereSmoothRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereSmoothRadius[%"FSYM"] = %"FSYM, &sphere,
		    &PhotonTestSphereSmoothRadius[sphere]);
    if (sscanf(line, "PhotonTestSphereRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereRadius[%"ISYM"] = %"FSYM, &sphere,
		    &PhotonTestSphereRadius[sphere]);
    if (sscanf(line, "PhotonTestSphereCoreRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereCoreRadius[%"ISYM"] = %"FSYM, &sphere,
		    &PhotonTestSphereCoreRadius[sphere]);
    if (sscanf(line, "PhotonTestSphereDensity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereDensity[%"ISYM"] = %"FSYM, &sphere,
		    &PhotonTestSphereDensity[sphere]);
    if (sscanf(line, "PhotonTestSphereTemperature[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereTemperature[%"ISYM"] = %"FSYM, &sphere,
		    &PhotonTestSphereTemperature[sphere]);
    if (sscanf(line, "PhotonTestSpherePosition[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSpherePosition[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM, 
		    &sphere, &PhotonTestSpherePosition[sphere][0],
		    &PhotonTestSpherePosition[sphere][1],
		    &PhotonTestSpherePosition[sphere][2]);
    if (sscanf(line, "PhotonTestSphereVelocity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereVelocity[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
		    &sphere, &PhotonTestSphereVelocity[sphere][0],
		    &PhotonTestSphereVelocity[sphere][1],
		    &PhotonTestSphereVelocity[sphere][2]);
    if (sscanf(line, "PhotonTestFracKeplerianRot[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestFracKeplerianRot[%"ISYM"] = %"FSYM, &sphere,
                    &PhotonTestFracKeplerianRot[sphere]);
    if (sscanf(line, "PhotonTestSphereTurbulence[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereTurbulence[%"ISYM"] = %"FSYM, &sphere,
                    &PhotonTestSphereTurbulence[sphere]);
    if (sscanf(line, "PhotonTestSphereCutOff[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereCutOff[%"ISYM"] = %"FSYM, &sphere,
                    &PhotonTestSphereCutOff[sphere]);
    if (sscanf(line, "PhotonTestSphereAng1[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereAng1[%"ISYM"] = %"FSYM, &sphere,
                    &PhotonTestSphereAng1[sphere]);
    if (sscanf(line, "PhotonTestSphereAng2[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereAng2[%"ISYM"] = %"FSYM, &sphere,
                    &PhotonTestSphereAng2[sphere]);
    if (sscanf(line, "PhotonTestSphereNumShells[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereNumShells[%"ISYM"] = %"ISYM, &sphere,
                    &PhotonTestSphereNumShells[sphere]);
    if (sscanf(line, "PhotonTestSphereHIIFraction[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereHIIFraction[%"ISYM"] = %"FSYM, &sphere,
                    &PhotonTestSphereHIIFraction[sphere]);
    if (sscanf(line, "PhotonTestSphereHeIIFraction[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereHeIIFraction[%"ISYM"] = %"FSYM, &sphere,
                    &PhotonTestSphereHeIIFraction[sphere]);
    if (sscanf(line, "PhotonTestSphereHeIIIFraction[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereHeIIIFraction[%"ISYM"] = %"FSYM, &sphere,
                    &PhotonTestSphereHeIIIFraction[sphere]);
    if (sscanf(line, "PhotonTestSphereH2IFraction[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "PhotonTestSphereH2IFraction[%"ISYM"] = %"FSYM, &sphere,
                    &PhotonTestSphereH2IFraction[sphere]);

    ret += sscanf(line, "PhotonTestRefineAtStart = %"ISYM,
		  &PhotonTestRefineAtStart);
    ret += sscanf(line, "PhotonTimeStep = %"FSYM,
		  &dtPhoton);
    ret += sscanf(line, "PhotonTestOmegaBaryonNow = %"FSYM,
		  &PhotonTestOmegaBaryonNow);
    ret += sscanf(line, "PhotonTestInitialFractionHII = %"FSYM,
		  &PhotonTestInitialFractionHII);
    ret += sscanf(line, "PhotonTestInitialFractionHeII = %"FSYM,
		  &PhotonTestInitialFractionHeII);
    ret += sscanf(line, "PhotonTestInitialFractionHeIII = %"FSYM,
		  &PhotonTestInitialFractionHeIII);
    ret += sscanf(line, "PhotonTestInitialFractionHM = %"FSYM,
		  &PhotonTestInitialFractionHM);
    ret += sscanf(line, "PhotonTestInitialFractionH2I = %"FSYM,
		  &PhotonTestInitialFractionH2I);
    ret += sscanf(line, "PhotonTestInitialFractionH2II = %"FSYM,
		  &PhotonTestInitialFractionH2II);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && 
	(strstr(line, "PhotonTest") && (strstr(line, "Source") == NULL) 
	 && line[0] != '#'))
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "warning0: %"ISYM", the following parameter line was not interpreted:\n%s\n", ret, line);
    
  } // end input from parameter file

  /* define sources temporarily in case we want to refine by optical
     depth.  We read them in permanently in RadiativeTransferInitialize() */ 

  rewind(fptr);
  if (ProblemType == 50)
    if (ReadPhotonSources(fptr, MetaData.Time) == FAIL) {
      ENZO_FAIL("Error in ReadPhotonSources.\n");
    }
  
  PhotonTime = InitialTimeInCodeUnits;

  if (PhotonTestDensityFilename != NULL &&
      ParallelRootGridIO == TRUE)
    ENZO_FAIL("PRGIO and external density files aren't supported yet.");
  
  /* set up grid */

  int RefineByOpticalDepth = FALSE;
  TotalRefinement = (ParallelRootGridIO == TRUE && Reinitialize) ? -1 : 1;

  HierarchyEntry *CurrentGrid = &TopGrid;

  while (CurrentGrid) {
    CurrentGrid->GridData->PhotonTestInitializeGrid(
	     PhotonTestNumberOfSpheres, PhotonTestSphereRadius,
	     PhotonTestSphereCoreRadius, PhotonTestSphereDensity,
	     PhotonTestSphereTemperature,
	     PhotonTestSpherePosition, PhotonTestSphereVelocity,
             PhotonTestFracKeplerianRot, PhotonTestSphereTurbulence,
             PhotonTestSphereCutOff, PhotonTestSphereAng1,
             PhotonTestSphereAng2, PhotonTestSphereNumShells,
	     PhotonTestSphereType, PhotonTestSphereConstantPressure,
	     PhotonTestSphereSmoothSurface, PhotonTestSphereSmoothRadius,
	     PhotonTestSphereHIIFraction, PhotonTestSphereHeIIFraction,
	     PhotonTestSphereHeIIIFraction, PhotonTestSphereH2IFraction,
	     PhotonTestUseParticles,
             PhotonTestUniformVelocity, PhotonTestUseColour,
             PhotonTestInitialTemperature, 0, 
	     PhotonTestInitialFractionHII, PhotonTestInitialFractionHeII,
	     PhotonTestInitialFractionHeIII, PhotonTestInitialFractionHM,
	     PhotonTestInitialFractionH2I, PhotonTestInitialFractionH2II, 
	     RefineByOpticalDepth, TotalRefinement, PhotonTestDensityFilename);

    CurrentGrid = CurrentGrid->NextGridThisLevel;

  } // ENDWHILE

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }


  /* If requested, refine the grid to the desired level. */

  if (PhotonTestRefineAtStart &&
      ((ParallelRootGridIO == FALSE) ||
       (ParallelRootGridIO == TRUE && Reinitialize))) {

    if (PhotonTestDensityFilename != NULL)
      ENZO_FAIL("External density field not supported with RefineAtStart yet.");

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	ENZO_FAIL("Error in RebuildHierarchy.\n");
      }
      if (LevelArray[level+1] == NULL)
	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];

      if (Reinitialize) 
	TotalRefinement = -1;
      else
	TotalRefinement = nint(POW(RefineBy, level+1));

      while (Temp != NULL) {
	if (Temp->GridData->PhotonTestInitializeGrid(
	     PhotonTestNumberOfSpheres, PhotonTestSphereRadius,
	     PhotonTestSphereCoreRadius, PhotonTestSphereDensity,
	     PhotonTestSphereTemperature,
	     PhotonTestSpherePosition, PhotonTestSphereVelocity,
             PhotonTestFracKeplerianRot, PhotonTestSphereTurbulence,
             PhotonTestSphereCutOff, PhotonTestSphereAng1,
             PhotonTestSphereAng2, PhotonTestSphereNumShells,
	     PhotonTestSphereType, PhotonTestSphereConstantPressure,
	     PhotonTestSphereSmoothSurface, PhotonTestSphereSmoothRadius,
	     PhotonTestSphereHIIFraction, PhotonTestSphereHeIIFraction,
	     PhotonTestSphereHeIIIFraction, PhotonTestSphereH2IFraction,
	     PhotonTestUseParticles,
	     PhotonTestUniformVelocity, PhotonTestUseColour,
	     PhotonTestInitialTemperature, level+1, 
	     PhotonTestInitialFractionHII, PhotonTestInitialFractionHeII,
	     PhotonTestInitialFractionHeIII, PhotonTestInitialFractionHM,
	     PhotonTestInitialFractionH2I, PhotonTestInitialFractionH2II,
	     RefineByOpticalDepth, TotalRefinement, NULL) == FAIL) {
	  ENZO_FAIL("Error in PhotonTestInitializeGrid.\n");
	}
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* After we refine on criteria other than optical depth, refine
       again with a radiation field starting at the maximum level
       reached before */

    for (i = 0; i < MAX_FLAGGING_METHODS; i++)
      if (CellFlaggingMethod[i] == 9) {
	RefineByOpticalDepth = TRUE;
	break;
      }

    if (RefineByOpticalDepth) {

      int level2;
      for (level2 = level-1; level2 < MaximumRefinementLevel; level2++) {
	if (RebuildHierarchy(&MetaData, LevelArray, level2) == FAIL) {
	  ENZO_FAIL("Error in RebuildHierarchy.\n");
	}	

	if (LevelArray[level2+1] == NULL)
	  break;
	LevelHierarchyEntry *Temp = LevelArray[level2+1];

	if (Reinitialize) 
	  TotalRefinement = -1;
	else
	  TotalRefinement = nint(POW(RefineBy, level2+1));

	while (Temp != NULL) {
	  
	  if (Temp->GridData->PhotonTestInitializeGrid(
		    PhotonTestNumberOfSpheres, PhotonTestSphereRadius,
		    PhotonTestSphereCoreRadius, PhotonTestSphereDensity,
		    PhotonTestSphereTemperature,
		    PhotonTestSpherePosition, PhotonTestSphereVelocity,
		    PhotonTestFracKeplerianRot, PhotonTestSphereTurbulence,
		    PhotonTestSphereCutOff, PhotonTestSphereAng1,
		    PhotonTestSphereAng2, PhotonTestSphereNumShells,
		    PhotonTestSphereType, PhotonTestSphereConstantPressure,
		    PhotonTestSphereSmoothSurface, PhotonTestSphereSmoothRadius,
		    PhotonTestSphereHIIFraction, PhotonTestSphereHeIIFraction,
		    PhotonTestSphereHeIIIFraction, PhotonTestSphereH2IFraction,
		    PhotonTestUseParticles,
		    PhotonTestUniformVelocity, PhotonTestUseColour,
		    PhotonTestInitialTemperature, level2+1, 
		    PhotonTestInitialFractionHII, PhotonTestInitialFractionHeII,
		    PhotonTestInitialFractionHeIII, PhotonTestInitialFractionHM,
		    PhotonTestInitialFractionH2I, PhotonTestInitialFractionH2II,
		    RefineByOpticalDepth, TotalRefinement, NULL) == FAIL) {
	    ENZO_FAIL("Error in PhotonTestInitializeGrid.\n");
	  }
	  Temp = Temp->NextGridThisLevel;
	}
      } // end: loop over levels

      /* Loop back from the bottom, restoring the consistency among levels. */

      for (level = MaximumRefinementLevel; level > 0; level--) {
	LevelHierarchyEntry *Temp = LevelArray[level];
	while (Temp != NULL) {
	  if (Temp->GridData->ProjectSolutionToParentGrid(
				  *LevelArray[level-1]->GridData) == FAIL) {
	    ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.\n");
	  }
	  Temp = Temp->NextGridThisLevel;
	}
      }

    } /* ENDIF RefineByOpticalDepth */

  } // end: if (PhotonTestRefineAtStart)

  /* Delete the temporary sources */

  RadiationSourceEntry *RS = GlobalRadiationSources->NextSource;
  while (RS != NULL)
    RS = DeleteRadiationSource(RS);

  /* set up field names and units */

  if (!Reinitialize) {
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
  if (PhotonTestUseColour)
    DataLabel[count++] = (char*) ColourName;
  
  if (RadiativeTransfer)
    if (MultiSpecies) {
      DataLabel[count++]  = (char*) kphHIName;
      DataLabel[count++]  = (char*) gammaName;
      DataLabel[count++]  = (char*) kphHeIName;
      DataLabel[count++]  = (char*) kphHeIIName;
      if (MultiSpecies > 1) 
	DataLabel[count++]= (char*) kdissH2IName; 
    } // if RadiativeTransfer

  if (RadiationPressure) {
    DataLabel[count++]  = (char*) RadAccel1Name;
    DataLabel[count++]  = (char*) RadAccel2Name;
    DataLabel[count++]  = (char*) RadAccel3Name;
  }

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "PhotonTestNumberOfSpheres    = %"ISYM"\n",
	    PhotonTestNumberOfSpheres);
    fprintf(Outfptr, "PhotonTestRefineAtStart      = %"ISYM"\n",
	    PhotonTestRefineAtStart);
    fprintf(Outfptr, "PhotonTestUseParticles       = %"ISYM"\n",
	    PhotonTestUseParticles);
    fprintf(Outfptr, "PhotonTestUseColour          = %"ISYM"\n",
	    PhotonTestUseColour);
    fprintf(Outfptr, "PhotonTestInitialTemperature = %"FSYM"\n",
	    PhotonTestInitialTemperature);
    fprintf(Outfptr, "PhotonTestUniformVelocity    = %"FSYM" %"FSYM" %"FSYM"\n",
	    PhotonTestUniformVelocity[0], PhotonTestUniformVelocity[1],
	    PhotonTestUniformVelocity[2]);
    for (sphere = 0; sphere < PhotonTestNumberOfSpheres; sphere++) {
      fprintf(Outfptr, "PhotonTestSphereType[%"ISYM"] = %"ISYM"\n", sphere,
	      PhotonTestSphereType[sphere]);
      fprintf(Outfptr, "PhotonTestSphereConstantPressure[%"ISYM"] = %"ISYM"\n", sphere,
	      PhotonTestSphereConstantPressure[sphere]);
      fprintf(Outfptr, "PhotonTestSphereSmoothSurface[%"ISYM"] = %"ISYM"\n", sphere,
	      PhotonTestSphereSmoothSurface[sphere]);
      fprintf(Outfptr, "PhotonTestSphereSmoothRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      PhotonTestSphereSmoothRadius[sphere]);
      fprintf(Outfptr, "PhotonTestSphereRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      PhotonTestSphereRadius[sphere]);
      fprintf(Outfptr, "PhotonTestSphereCoreRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      PhotonTestSphereCoreRadius[sphere]);
      fprintf(Outfptr, "PhotonTestSphereDensity[%"ISYM"] = %"FSYM"\n", sphere,
	      PhotonTestSphereDensity[sphere]);
      fprintf(Outfptr, "PhotonTestSphereTemperature[%"ISYM"] = %"FSYM"\n", sphere,
	      PhotonTestSphereTemperature[sphere]);
      fprintf(Outfptr, "PhotonTestSpherePosition[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			PhotonTestSpherePosition[sphere]);
      fprintf(Outfptr, "PhotonTestSphereVelocity[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			PhotonTestSphereVelocity[sphere]);
      fprintf(Outfptr, "PhotonTestSphereHIIFraction[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      PhotonTestSphereHIIFraction[sphere]);
      fprintf(Outfptr, "PhotonTestSphereHeIIFraction[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      PhotonTestSphereHeIIFraction[sphere]);
      fprintf(Outfptr, "PhotonTestSphereHeIIIFraction[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      PhotonTestSphereHeIIIFraction[sphere]);
      fprintf(Outfptr, "PhotonTestSphereH2IFraction[%"ISYM"] = %"GOUTSYM"\n", sphere,
	      PhotonTestSphereH2IFraction[sphere]);
      fprintf(Outfptr, "PhotonTestFracKeplerianRot[%"ISYM"] = %"GOUTSYM"\n", sphere,
              PhotonTestFracKeplerianRot[sphere]);
      fprintf(Outfptr, "PhotonTestSphereTurbulence[%"ISYM"] = %"GOUTSYM"\n", sphere,
              PhotonTestSphereTurbulence[sphere]);
      fprintf(Outfptr, "PhotonTestSphereCutOff[%"ISYM"] = %"GOUTSYM"\n", sphere,
              PhotonTestSphereCutOff[sphere]);
      fprintf(Outfptr, "PhotonTestSphereAng1[%"ISYM"] = %"GOUTSYM"\n", sphere,
              PhotonTestSphereAng1[sphere]);
      fprintf(Outfptr, "PhotonTestSphereAng2[%"ISYM"] = %"GOUTSYM"\n", sphere,
              PhotonTestSphereAng2[sphere]);
      fprintf(Outfptr, "PhotonTestSphereNumShells[%"ISYM"] = %"ISYM"\n\n", sphere,
              PhotonTestSphereNumShells[sphere]);
    }
  }
  } // ENDIF !Reinitialize

  delete [] dummy;

  return SUCCESS;

}
