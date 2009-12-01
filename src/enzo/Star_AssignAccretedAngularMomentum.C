/***********************************************************************
/
/  ASSIGN ACCRETED ANGULAR MOMENTUM FROM THE FILE
/
/  written by: Ji-hoon Kim 
/  date:       November, 2009
/  modified1:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
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
#include "StarParticleData.h"

void Star::AssignAccretedAngularMomentum(void)
{

  if (CurrentGrid == NULL)
    return;

  FILE *fptr;
  int dim, dummy_int[2];
  float dummy[2], AccretedAngularMomentum[] = {0.0, 0.0, 0.0};
  char line[MAX_LINE_LENGTH];
  
  if ((fptr = fopen(MBHParticleIOFilename, "r")) == NULL) {
    fprintf(stderr, "Error opening file %s\n", MBHParticleIOFilename);
    fprintf(stderr, "Assume zero angular momentum accreted onto MBH so far.\n");
  } else {
    //naturally, the last line in the file matching the ID is used for angular momentum
    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) 
      if (line[0] != '#') {
	//order: time, regular star count, MBH id, MBH mass, MBH angular momentum 
	if (sscanf(line, " %"FSYM"  %"ISYM"  %"ISYM"  %"FSYM"  %"FSYM"  %"FSYM"  %"FSYM, 
		   &dummy[0], &dummy_int[0], &dummy_int[1], &dummy[1], 
		   &AccretedAngularMomentum[0], &AccretedAngularMomentum[1], 
		   &AccretedAngularMomentum[2]) != 7) {
	  fprintf(stderr, "File structure wrong: %s\n", MBHParticleIOFilename);
	  ENZO_FAIL("");
	}
	//go through the file and pick the right one that matches the ID
	if (Identifier == dummy_int[1])
	  for (dim = 0; dim < MAX_DIMENSION; dim++) 
	    accreted_angmom[dim] = AccretedAngularMomentum[dim];
      }
  }
  return;

}
