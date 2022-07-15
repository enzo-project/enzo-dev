/***********************************************************************
/
/  INITIALIZE A SINE WAVE DENSITY PROFILE TO TEST THE GRAVITY SOLVER
/
/  written by: JC Passy
/  date:       June, 2013
/
/  PURPOSE:
/    Set up a sine mass distribution to test the gravity solver.

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

/* function prototypes */
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
                     LevelHierarchyEntry *LevelArray[], int level);

#define MAX_INITIAL_GRIDS 10

int TestGravitySineWaveInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                                  TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *GPotName  = "Grav_Potential";

  /* parameter declarations */

  int   NumberOfSubgridZones[MAX_DIMENSION], SubgridDims[MAX_DIMENSION];
  float LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  
  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  i, dim, ret, SubgridsAreStatic;
 
  /* set default parameters */
  
  float TestGravitySineWaveAmplitude     = 1.0;
  float TestGravitySineWavePeriod        = 1.0;
  float TestGravitySineWaveAngle         = 0.0; // Inclination angle in degrees
  int TestGravitySineWaveRefineAtStart   = FALSE;

  /* Set default parameters: parameters, names and subgrid info */

  int   TestGravitySineWaveGridDimension[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  int   TestGravitySineWaveGridLevel[MAX_INITIAL_GRIDS];
  float TestGravitySineWaveGridLeftEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  float TestGravitySineWaveGridRightEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];

  for (i = 0; i < MAX_INITIAL_GRIDS; i++)
    TestGravitySineWaveGridLevel[i] = 1;

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    TestGravitySineWaveGridLeftEdge[0][dim] = DomainLeftEdge[dim];
    TestGravitySineWaveGridRightEdge[0][dim] = DomainRightEdge[dim];
    TestGravitySineWaveGridDimension[0][dim] = MetaData.TopGridDims[dim];
  }

  TestGravitySineWaveGridLevel[0] = 0;

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */
    ret += sscanf(line, "TestGravitySineWaveAmplitude = %"FSYM"", &TestGravitySineWaveAmplitude);
    ret += sscanf(line, "TestGravitySineWavePeriod = %"FSYM"", &TestGravitySineWavePeriod);
    ret += sscanf(line, "TestGravitySineWaveAngle = %"FSYM"", &TestGravitySineWaveAngle);
    ret += sscanf(line, "TestGravitySineWaveRefineAtStart = %"ISYM"", &TestGravitySineWaveRefineAtStart);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "TestGravitySineWave") &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, 
         "warning in TestGravitySineWaveInitialize.C: the following parameter line was not interpreted:\n%s\n",
	      line);

  } // end input from parameter file

  /* set the periodic boundaries */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }

  /* Set up grids from CollapseTestInitialize.C */

  if ((MetaData.StaticHierarchy) && (TestGravitySineWaveRefineAtStart))
    ENZO_FAIL("Error in TestGravitySineWaveInitialize.C: StaticHierarchy = %"ISYM", TestGravitySiveWaveRefineAtStart = %"ISYM" \n");
  
  SubgridsAreStatic = MetaData.StaticHierarchy;

  /* Initial grid */
  int level;
  level = 0;
  int TotalRefinement = nint(pow(FLOAT(RefineBy), TestGravitySineWaveGridLevel[level]));
  if (TopGrid.GridData->TestGravitySineWaveInitializeGrid(TestGravitySineWaveAmplitude,
                                                          TestGravitySineWavePeriod,
                                                          TestGravitySineWaveAngle,
                                                          SubgridsAreStatic,
                                                          TotalRefinement,
                                                          level
                                                          ) == FAIL){
    ENZO_FAIL("Error in grid->TestGravitySineWaveInitializeGrid.\n");
  }  

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
        float(MetaData.TopGridDims[dim]);
  }

  /* If requested, refine the grid to the desired level. */

  if (TestGravitySineWaveRefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL)
        ENZO_FAIL("Error in RebuildHierarchy.\n");

      if (LevelArray[level+1] == NULL)
        break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
        TotalRefinement = nint(pow(FLOAT(RefineBy),TestGravitySineWaveGridLevel[level+1]));
        if (Temp->GridData->TestGravitySineWaveInitializeGrid(TestGravitySineWaveAmplitude,
                                                              TestGravitySineWavePeriod,
                                                              TestGravitySineWaveAngle,
                                                              SubgridsAreStatic,
                                                              TotalRefinement,
                                                              level+1
                                                              ) == FAIL){
          ENZO_FAIL("Error in grid->TurbulenceICMInitializeGrid.\n");
        }
        Temp = Temp->NextGridThisLevel;
      }
    } /* end loop over levels */

    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
        if (Temp->GridData->ProjectSolutionToParentGrid(*LevelArray[level-1]->GridData) == FAIL)
          ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.\n");
        Temp = Temp->NextGridThisLevel;
      }
    }

  } /* End if (TestGravitySineWaveRefineAtStart) */

  /* set up fields name and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  if (WritePotential)
    DataLabel[count++] = GPotName;

  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  DataUnits[5] = NULL;
  DataUnits[6] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "TestGravitySineWaveAmplitude = %"FSYM"\n", TestGravitySineWaveAmplitude);
    fprintf(Outfptr, "TestGravitySineWavePeriod = %"FSYM"\n", TestGravitySineWavePeriod);
    fprintf(Outfptr, "TestGravitySineWaveAngle = %"FSYM"\n", TestGravitySineWaveAngle);
    fprintf(Outfptr, "\n");
  }

  return SUCCESS;
}
