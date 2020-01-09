/***********************************************************************
/
/  CHECKS GRAVITATING MASS FIELDS
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
 
 
char TGOutputFileName3[] = "GMF.out";
char TGOutputFileName4[] = "GMFP.out";
 
 
int OutputGravitatingMassField(HierarchyEntry *Grid, int level,int cycle)
{
 
  /* declarations */
 
  FILE *fptr;
  char name[MAX_LINE_LENGTH], proc[MAX_TASK_TAG_SIZE], cycle_name[MAX_TASK_TAG_SIZE];;
  FILE *fptr2;
  char name2[MAX_LINE_LENGTH], proc2[MAX_TASK_TAG_SIZE], cycle_name2[MAX_TASK_TAG_SIZE];;
 
  /* Open output file. */
 
  strcpy(name, TGOutputFileName3);
  sprintf(cycle_name, "%"TASK_TAG_FORMAT""ISYM, cycle);
  strcat(name, cycle_name);
  sprintf(proc, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
  strcat(name, "_p");
  strcat(name, proc);

  if ((fptr = fopen(name, "a")) == NULL) {
    fprintf(stderr, "OutputGravitatingMassField:\n  Error opening file %s.\n",
	    name);
    exit(FAIL);
  }

  strcpy(name2, TGOutputFileName4);
  sprintf(cycle_name2, "%"TASK_TAG_FORMAT""ISYM, cycle);
  strcat(name2, cycle_name2);
  sprintf(proc2, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
  strcat(name2, "_p");
  strcat(name2, proc2);

  if ((fptr2 = fopen(name2, "a")) == NULL) {
    fprintf(stderr, "OutputGravitatingMassFieldParticles:\n  Error opening file %s.\n",
	    name2);
    exit(FAIL);
  }
 
  /* For each grid on each level, check the results. */
 
  if (Grid->GridData->OutputGravitatingMassField(fptr, fptr2, level) == FAIL)
    ENZO_FAIL("Error in grid->OutputGravitatingMassField\n");

  /* Close output file. */
  
  fclose(fptr);
  fclose(fptr2);
 
  return SUCCESS;
}
