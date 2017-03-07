#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
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

int InitializeTimeVaryingExternalAcceleration(float Time){

  // this doesn't concern us
  if (ExternalGravity < 2 || ExternalGravity > 3){
    return SUCCESS;
  }

  // if beyond the time off, do nothing
  if (Time > ExternalGravityTimeOff)
    return SUCCESS;

  // check if already initialized
  if (ExternalGravityTime != NULL){
    return SUCCESS;
  }

  // read from file
  char line[MAX_LINE_LENGTH];
  FILE *fptr = fopen("external_gravity.in", "r");
  if (fptr == NULL)
    ENZO_FAIL("Error opening time varying external gravity file, 'external_gravity.in' \n");

  ExternalGravityTime = new float[EXTERNAL_GRAVITY_ENTRIES];
  for (int dim = 0; dim < MAX_DIMENSION; dim++){
    ExternalGravityTimePositions[dim] = new float[EXTERNAL_GRAVITY_ENTRIES];
  }

  ExternalGravityNumberofTimePoints = 0;
  int i = 0, err;
  while( fgets(line, MAX_LINE_LENGTH, fptr) != NULL){
    if (line[0] != '#'){
      err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM,
                     &ExternalGravityTime[i],
                     &ExternalGravityTimePositions[0][i],
                     &ExternalGravityTimePositions[1][i],
                     &ExternalGravityTimePositions[2][i]);
    }
    i++;
  }
  ExternalGravityNumberofTimePoints = i - 1;
  fclose(fptr);

  // do a check
  if (   (ExternalGravityTimeOff - ExternalGravityTimeOn)
                 > ExternalGravityTime[ExternalGravityNumberofTimePoints - 1]){
    printf("WARNING: ExternalGravityTimeOff set to a value beyond the available interp bins. Correcting");

    ExternalGravityTimeOff = ExternalGravityTime[ExternalGravityNumberofTimePoints - 1] +
                            ExternalGravityTimeOn;
  }

  return SUCCESS;
}
