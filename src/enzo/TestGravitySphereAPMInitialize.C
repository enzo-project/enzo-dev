/***********************************************************************
/
/  INITIALIZE A GRAVITY TEST
/
/  written by: Greg Bryan
/  date:       September, 1995
/  modified1: Jean-Claude Passy, June 2013
/
/  PURPOSE:
/    Set up a spherical mass distribution to test the gravity solver.
/
/  UPDATE: Following tests performed in Passy & Bryan (2014).
/          Possible distributions are:
/           - uniform density (type 0)
/           - isothermal spere (type 1)
/           - Plummer sphere (type 2)
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
#include "LevelHierarchy.h"
#include "TopGridData.h"

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
                     LevelHierarchyEntry *LevelArray[], int level);

int TestGravitySphereAPMInitialize(FILE *fptr, FILE *Outfptr,
                                   HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *GEName = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *GPotName = "Grav_Potential";

  /* declarations */

  char line[MAX_LINE_LENGTH];
  int dim, ret, level;

  /* Error check. */

  if (!SelfGravity)
    fprintf(stderr, "TestGravitySphere: gravity is off!?!");

  /* set default parameters */

  float TestGravitySphereInteriorDensity = 1.0;         // density inside sphere
  float TestGravitySphereExteriorDensity = tiny_number; // density outside sphere
  float TestGravitySphereRadius = 0.1;
  int TestGravitySphereType = 0; // uniform density
  int TestGravitySphereUseBaryons = TRUE;
  int TestGravitySphereRefineAtStart = FALSE;
  FLOAT TestGravitySphereCenter[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    TestGravitySphereCenter[dim] = 0.5 * (DomainLeftEdge[dim] + DomainRightEdge[dim]);

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
  {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "TestGravitySphereInteriorDensity = %" FSYM,
                  &TestGravitySphereInteriorDensity);
    ret += sscanf(line, "TestGravitySphereExteriorDensity = %" FSYM,
                  &TestGravitySphereExteriorDensity);
    ret += sscanf(line, "TestGravitySphereRadius = %" FSYM,
                  &TestGravitySphereRadius);
    ret += sscanf(line, "TestGravitySphereType = %" ISYM,
                  &TestGravitySphereType);
    ret += sscanf(line, "TestGravitySphereUseBaryons = %" ISYM,
                  &TestGravitySphereUseBaryons);
    ret += sscanf(line, "TestGravitySphereRefineAtStart = %" ISYM,
                  &TestGravitySphereRefineAtStart);
    ret += sscanf(line, "TestGravitySphereCenter = %" PSYM " %" PSYM " %" PSYM,
                  TestGravitySphereCenter, TestGravitySphereCenter + 1, TestGravitySphereCenter + 2);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "TestGravitySphere") && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up grid */

  if (TopGrid.GridData->TestGravitySphereAPMInitializeGrid(TestGravitySphereInteriorDensity,
                                                           TestGravitySphereExteriorDensity,
                                                           TestGravitySphereRadius,
                                                           TestGravitySphereType,
                                                           TestGravitySphereUseBaryons,
                                                           TestGravitySphereCenter) == FAIL)
  {
    ENZO_FAIL("Error in TestGravitySphereAPMInitializeGrid.");
  }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED)
  {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *= (DomainRightEdge[dim] - DomainLeftEdge[dim]) / float(MetaData.TopGridDims[dim]);
  }

  /* If requested, refine the grid to the desired level. */

  if (TestGravitySphereRefineAtStart)
  {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++)
    {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL)
      {
        ENZO_FAIL("Error in RebuildHierarchy.");
      }
      if (LevelArray[level + 1] == NULL)
        break;
      LevelHierarchyEntry *Temp = LevelArray[level + 1];
      while (Temp != NULL)
      {
        if (Temp->GridData->TestGravitySphereAPMInitializeGrid(TestGravitySphereInteriorDensity,
                                                               TestGravitySphereExteriorDensity,
                                                               TestGravitySphereRadius,
                                                               TestGravitySphereType,
                                                               TestGravitySphereUseBaryons,
                                                               TestGravitySphereCenter) == FAIL)
        {
          ENZO_FAIL("Error in TestGravitySphereAPMInitializeGrid.")
        }
        Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* Loop back from the bottom, restoring the consistency amoung levels. */

    for (level = MaximumRefinementLevel; level > 0; level--)
    {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL)
      {
        if (Temp->GridData->ProjectSolutionToParentGrid(*Temp->GridHierarchyEntry->ParentGrid->GridData) == FAIL)
          ENZO_FAIL("Error in grid->ProjectSolutionToParentGrid.");
        Temp = Temp->NextGridThisLevel;
      }
    }
  } // end: if (TestGravitySphereRefineAtStart)

  /* set up field names and units */

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

  for (int j = 0; j < count; j++)
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR)
  {
    fprintf(Outfptr, "TestGravitySphereInteriorDensity   = %" ESYM "\n",
            TestGravitySphereInteriorDensity);
    fprintf(Outfptr, "TestGravitySphereExteriorDensity   = %" ESYM "\n",
            TestGravitySphereExteriorDensity);
    fprintf(Outfptr, "TestGravitySphereRadius            = %" FSYM "\n",
            TestGravitySphereRadius);
    fprintf(Outfptr, "TestGravitySphereType              = %" ISYM "\n",
            TestGravitySphereType);
    fprintf(Outfptr, "TestGravitySphereUseBaryons        = %" ISYM "\n",
            TestGravitySphereUseBaryons);
    fprintf(Outfptr, "TestGravitySphereRefineAtStart     = %" ISYM "\n",
            TestGravitySphereRefineAtStart);
    fprintf(Outfptr, "TestGravitySphereCenter =  % " GOUTSYM " %" GOUTSYM " %" GOUTSYM "\n\n",
            TestGravitySphereCenter[0],
            TestGravitySphereCenter[1],
            TestGravitySphereCenter[2]);
  }

  return SUCCESS;
}
