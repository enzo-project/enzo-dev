/***********************************************************************
/
/  CHECKS THE ACCURACY OF THE GRAVITY SOLVER
/
/  written by: JC Passy
/  date:       March 2013
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
 
 
char TGOutputFileName2[] = "AccelerationField.out";
 
 
int OutputAccelerationField(HierarchyEntry *Grid, int level, int cycle)
{
 
  /* declarations */
 
  FILE *fptr;
  char name[MAX_LINE_LENGTH], proc[MAX_TASK_TAG_SIZE], cycle_name[MAX_TASK_TAG_SIZE];
 
  /* Open output file. */
 
  strcpy(name, TGOutputFileName2);
  sprintf(cycle_name, "%"TASK_TAG_FORMAT""ISYM, cycle);
  strcat(name, cycle_name);
  sprintf(proc, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
  strcat(name, "_p");
  strcat(name, proc);

  if ((fptr = fopen(name, "a")) == NULL) {
    fprintf(stderr, "CheckAccelerationField:\n  Error opening file %s.\n",
	    name);
    exit(FAIL);
  }
 
  /* For each grid on each level, check the results. */
 
  if (Grid->GridData->OutputAccelerationField(fptr, level) == FAIL)
    ENZO_FAIL("Error in grid->OutputAccelerationField\n");

  /* Close output file. */
  
  fclose(fptr);
 
  return SUCCESS;
}
