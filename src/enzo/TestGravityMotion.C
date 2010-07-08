/***********************************************************************
/
/  INITIALIZE A TEST OF PARTICLE MOTION
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:
/
/  PURPOSE:
/    We set up a system in there is one or more particles in motion.
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
#include "TopGridData.h"
 
int TestGravityMotion(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		      TopGridData &MetaData)
{
 
  /* declarations */
 
  char  line[MAX_LINE_LENGTH];
  int   ret;
 
  /* Error check. */
 
  if (!SelfGravity)
    fprintf(stderr, "TestGravityMotion: gravity is off!?!");
 
  /* set default parameters */
 
  float TestGravityParticleVelocity  = 1.0;
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters */
 
    ret += sscanf(line, "TestGravityMotionParticleVelocity = %"PSYM,
		  &TestGravityParticleVelocity);
 
    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "TestGravityMotion")
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  } // end input from parameter file
 
  /* set up grid */
 
  if (TopGrid.GridData->TestGravityMotionInitializeGrid(
                                  TestGravityParticleVelocity) == FAIL) {
    ENZO_FAIL("Error in TestGravityMotionInitializeGrid.\n");
  }
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "TestGravityMotionParticleVelocity = %"GSYM"\n",
	    TestGravityParticleVelocity);
  }
 
  return SUCCESS;
 
}
