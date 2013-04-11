/***********************************************************************
/
/  INITIALIZE A Cool Core Cluster 
/
/  written by: Yuan Li and Greg Bryan
/  date:       Dec 2011 
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include "preincludes.h"
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
int CommunicationPartitionGrid(HierarchyEntry *Grid, int gridnum);


int ClusterInitialize(FILE *fptr, FILE *Outfptr, 
                           HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ColourName = "colour";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int ClusterNumberOfSpheres = 1;
  int ClusterRefineAtStart   = TRUE;
  int ClusterUseParticles    = FALSE;
  int ClusterUseColour       = FALSE;
  float ClusterInitialTemperature = 1000;
  float ClusterInitialSpinParameter = 0.05;
  int   ClusterSphereType[MAX_SPHERES];
  float ClusterSphereDensity[MAX_SPHERES],
        ClusterSphereTemperature[MAX_SPHERES],
        ClusterSphereVelocity[MAX_SPHERES][MAX_DIMENSION],
        ClusterUniformVelocity[MAX_DIMENSION];
  FLOAT ClusterSphereRadius[MAX_SPHERES],
        ClusterSphereCoreRadius[MAX_SPHERES],
        ClusterSpherePosition[MAX_SPHERES][MAX_DIMENSION];

  for (sphere = 0; sphere < MAX_SPHERES; sphere++) {
    ClusterSphereRadius[sphere]     = 1.0;
    ClusterSphereCoreRadius[sphere] = 0.1;
    ClusterSphereDensity[sphere]    = 1.0;
    ClusterSphereTemperature[sphere] = 1.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      ClusterSpherePosition[sphere][dim] = 0.5*(DomainLeftEdge[dim] +
                                                     DomainRightEdge[dim]);
      ClusterSphereVelocity[sphere][dim] = 0;
    }
    ClusterSphereType[sphere]       = 0;
  }
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    ClusterUniformVelocity[dim] = 0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "ClusterNumberOfSpheres = %"ISYM,
                  &ClusterNumberOfSpheres);
    ret += sscanf(line, "ClusterRefineAtStart = %"ISYM, 
                  &ClusterRefineAtStart);
    ret += sscanf(line, "ClusterUseParticles = %"ISYM, 
                  &ClusterUseParticles);
    ret += sscanf(line, "ClusterUseColour = %"ISYM, 
                  &ClusterUseColour);
    ret += sscanf(line, "ClusterInitialTemperature = %"FSYM, 
                  &ClusterInitialTemperature);
    ret += sscanf(line, "ClusterInitialSpinParameter = %"FSYM,
		  &ClusterInitialSpinParameter);
    ret += sscanf(line, "ClusterUniformVelocity = %"FSYM" %"FSYM" %"FSYM, 
                  ClusterUniformVelocity, ClusterUniformVelocity+1,
                  ClusterUniformVelocity+2);
    if (sscanf(line, "ClusterSphereType[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "ClusterSphereType[%"ISYM"] = %"ISYM, &sphere,
                    &ClusterSphereType[sphere]);
    if (sscanf(line, "ClusterSphereRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "ClusterSphereRadius[%"ISYM"] = %"PSYM, &sphere,
                    &ClusterSphereRadius[sphere]);
    if (sscanf(line, "ClusterSphereCoreRadius[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "ClusterSphereCoreRadius[%"ISYM"] = %"PSYM, &sphere,
                    &ClusterSphereCoreRadius[sphere]);
    if (sscanf(line, "ClusterSphereDensity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "ClusterSphereDensity[%"ISYM"] = %"FSYM, &sphere,
                    &ClusterSphereDensity[sphere]);
    if (sscanf(line, "ClusterSphereTemperature[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "ClusterSphereTemperature[%"ISYM"] = %"FSYM, &sphere,
                    &ClusterSphereTemperature[sphere]);
    if (sscanf(line, "ClusterSpherePosition[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "ClusterSpherePosition[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM, 
                    &sphere, &ClusterSpherePosition[sphere][0],
                    &ClusterSpherePosition[sphere][1],
                    &ClusterSpherePosition[sphere][2]);
    if (sscanf(line, "ClusterSphereVelocity[%"ISYM"]", &sphere) > 0)
      ret += sscanf(line, "ClusterSphereVelocity[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
                    &sphere, &ClusterSphereVelocity[sphere][0],
                    &ClusterSphereVelocity[sphere][1],
                    &ClusterSphereVelocity[sphere][2]);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "Cluster") 
        && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up grid */

  if (TopGrid.GridData->ClusterInitializeGrid(
             ClusterNumberOfSpheres, ClusterSphereRadius,
             ClusterSphereCoreRadius, ClusterSphereDensity,
             ClusterSphereTemperature, 
             ClusterSpherePosition, ClusterSphereVelocity,
             ClusterSphereType, ClusterUseParticles,
             ClusterUniformVelocity, ClusterUseColour,
             ClusterInitialTemperature, ClusterInitialSpinParameter, 0) == FAIL) {
    fprintf(stderr, "Error in ClusterInitializeGrid.\n");
    return FAIL;
  }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
        float(MetaData.TopGridDims[dim]);
  }

  Exterior.Prepare(TopGrid.GridData);

  float Dummy[MAX_DIMENSION];
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    if (Exterior.InitializeExternalBoundaryFace(dim, 
                                                MetaData.LeftFaceBoundaryCondition[dim],
                                                MetaData.RightFaceBoundaryCondition[dim],
                                                Dummy, Dummy)
        == FAIL) {
      fprintf(stderr, "Error in InitializeExternalBoundaryFace.\n");
      return FAIL;
    }

  /* Initialize particle boundary conditions. */

  Exterior.InitializeExternalBoundaryParticles(MetaData.ParticleBoundaryType);
  
  CommunicationPartitionGrid(&TopGrid, 0);

  /* If requested, refine the grid to the desired level. */

  if (ClusterRefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
        fprintf(stderr, "Error in RebuildHierarchy.\n");
        return FAIL;
      }
      if (LevelArray[level+1] == NULL)
        break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
        if (Temp->GridData->ClusterInitializeGrid(
             ClusterNumberOfSpheres, ClusterSphereRadius,
             ClusterSphereCoreRadius, ClusterSphereDensity,
             ClusterSphereTemperature, 
             ClusterSpherePosition, ClusterSphereVelocity,
             ClusterSphereType, ClusterUseParticles,
             ClusterUniformVelocity, ClusterUseColour,
             ClusterInitialTemperature, ClusterInitialSpinParameter, level+1) == FAIL) {
          fprintf(stderr, "Error in ClusterInitializeGrid.\n");
          return FAIL;
        }
        Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
        if (Temp->GridData->ProjectSolutionToParentGrid(*Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL) {
//                                 *LevelArray[level-1]->GridData) == FAIL) {
          fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
          return FAIL;
        }
        Temp = Temp->NextGridThisLevel;
      }
    }

  } // end: if (ClusterRefineAtStart)

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  if (ClusterUseColour)
    DataLabel[count++] = ColourName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "ClusterNumberOfSpheres    = %"ISYM"\n",
            ClusterNumberOfSpheres);
    fprintf(Outfptr, "ClusterRefineAtStart      = %"ISYM"\n",
            ClusterRefineAtStart);
    fprintf(Outfptr, "ClusterUseParticles       = %"ISYM"\n",
            ClusterUseParticles);
    fprintf(Outfptr, "ClusterUseColour          = %"ISYM"\n",
            ClusterUseColour);
    fprintf(Outfptr, "ClusterInitialTemperature = %"GOUTSYM"\n",
            ClusterInitialTemperature);
    fprintf(fptr, "ClusterInitialSpinParameter   = %"FSYM"\n",
	    ClusterInitialSpinParameter);
    fprintf(Outfptr, "ClusterUniformVelocity    = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
            ClusterUniformVelocity[0], ClusterUniformVelocity[1],
            ClusterUniformVelocity[2]);
    for (sphere = 0; sphere < ClusterNumberOfSpheres; sphere++) {
      fprintf(Outfptr, "ClusterSphereType[%"ISYM"] = %"ISYM"\n", sphere,
              ClusterSphereType[sphere]);
      fprintf(Outfptr, "ClusterSphereRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
              ClusterSphereRadius[sphere]);
      fprintf(Outfptr, "ClusterSphereCoreRadius[%"ISYM"] = %"GOUTSYM"\n", sphere,
              ClusterSphereCoreRadius[sphere]);
      fprintf(Outfptr, "ClusterSphereDensity[%"ISYM"] = %"GOUTSYM"\n", sphere,
              ClusterSphereDensity[sphere]);
      fprintf(Outfptr, "ClusterSphereTemperature[%"ISYM"] = %"GOUTSYM"\n", sphere,
              ClusterSphereTemperature[sphere]);
      fprintf(Outfptr, "ClusterSpherePosition[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
                        ClusterSpherePosition[sphere]);
      fprintf(Outfptr, "ClusterSphereVelocity[%"ISYM"] = ", sphere);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
                        ClusterSphereVelocity[sphere]);
    }
  }

  return SUCCESS;

}

