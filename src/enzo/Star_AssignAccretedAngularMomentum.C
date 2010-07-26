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
#include "LevelHierarchy.h"

void Star::AssignAccretedAngularMomentum(void)
{

  if (CurrentGrid == NULL)
    return;

  FILE *fptr;
  int dim, dummy_int[2];
  float dummy[2], AccretedAngularMomentum[] = {0.0, 0.0, 0.0};
  double dummy_double[2];
  char line[MAX_LINE_LENGTH];
  
  if ((fptr = fopen(MBHParticleIOFilename, "r")) == NULL) {
    fprintf(stderr, "Error opening file %s\n", MBHParticleIOFilename);
    fprintf(stderr, "Assume zero angular momentum accreted onto MBH so far.\n");
  } else {

    // naturally, the last line in the file matching the ID is used for angular momentum
    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) 
      if (line[0] != '#') {
	// order: time, regular star count, MBH id, MBH mass, MBH angular momentum 
	if (sscanf(line, " %"FSYM"  %"ISYM"  %"ISYM"  %lf  %"FSYM"  %"FSYM"  %"FSYM"  %lf", 
		   &dummy[0], &dummy_int[0], &dummy_int[1], &dummy_double[0], 
		   &AccretedAngularMomentum[0], &AccretedAngularMomentum[1], 
		   &AccretedAngularMomentum[2], &dummy_double[1]) != 8) {
	  ENZO_VFAIL("star::AAAM: File structure wrong: %s\n", MBHParticleIOFilename)
	}

	// go through the file and pick the right one that matches the ID
	if (Identifier == dummy_int[1]) {

	  // accreted angular momentum
	  for (dim = 0; dim < MAX_DIMENSION; dim++) 
	    accreted_angmom[dim] = AccretedAngularMomentum[dim];
	  // not ejected mass yet
	  NotEjectedMass = dummy_double[1];
//	  NotEjectedMass = fmod((double)(Mass * MBHFeedbackMassEjectionFraction), 
//				(double)(MBHFeedbackJetsThresholdMass));	  
	}
      }

    fclose(fptr);
  }

  return;

}
