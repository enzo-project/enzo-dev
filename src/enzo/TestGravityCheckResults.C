/***********************************************************************
/
/  CHECKS THE ACCURACY OF THE GRAVITY SOLVER
/
/  written by: Greg Bryan
/  date:       July, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
 
/* function prototypes */
 
 
char TGOutputFileName[] = "TestGravityCheckResults.out";
 
 
int TestGravityCheckResults(LevelHierarchyEntry *LevelArray[])
{
 
  /* declarations */
 
  int level;
  FILE *fptr;
  char name[MAX_LINE_LENGTH], proc[MAX_TASK_TAG_SIZE];
 
  /* Open output file. */
 
  strcpy(name, TGOutputFileName);
  if (NumberOfProcessors > 1) {
    sprintf(proc, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
    strcat(name, proc);
  }
  if ((fptr = fopen(name, "w")) == NULL) {
    fprintf(stderr, "TestGravityCheckResults:\n  Error opening file %s.\n",
	    name);
    exit(FAIL);
  }
 
  /* Output header. */
 
  fprintf(fptr, "#     r         f_tang         f_radial        f_analytic\n");
 
  /* For each grid on each level, check the results. */
 
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    LevelHierarchyEntry* Level = LevelArray[level];
 
    while (Level != NULL) {
 
      if (Level->GridData->TestGravityCheckResults(fptr,
				          LevelArray[0]->GridData) == FAIL) {
	fprintf(stderr, "Error in grid->TestGravityCheckResults\n");
	exit(FAIL);
      }
 
      /* Next grid on the this level. */
 
      Level = Level->NextGridThisLevel;
 
    }
 
  } // end: loop over levels
 
  /* Close output file. */
 
  fclose(fptr);
 
  return SUCCESS;
}
