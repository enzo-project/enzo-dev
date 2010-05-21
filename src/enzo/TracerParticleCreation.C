/***********************************************************************
/
/  CHECKS PARAMETER FILE FOR TRACER PARTICLE CREATION PARAMETERS AND
/    IF PRESENT, CREATES TRACER PARTICLES.
/
/  written by: Greg Bryan
/  date:       April, 2004
/  modified1:  Robert Harkness
/  date:       September, 2004
/
/  PURPOSE:  This routine scans through the parameter file pointed to
/     by the file pointer passed in and looks for parameters specifying
/     that tracer particles should be created.  If present, it creates
/     the tracer particles.  Note that unlike most parameters, these
/     parameters are not recorded because tracer particle creation is
/     a one-time event.
/
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
 
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
 
 
 
 
int TracerParticleCreation(FILE *fptr, HierarchyEntry &TopGrid,
			   TopGridData &MetaData)
{
 
  char line[MAX_LINE_LENGTH];
  int dim;
 
  // Set default values for parameters
 
  //  Declared in global_data.h (RH)
  //  FLOAT TracerParticleCreationLeftEdge[MAX_DIMENSION];
  //  FLOAT TracerParticleCreationRightEdge[MAX_DIMENSION];
  //  FLOAT TracerParticleCreationSpacing = FLOAT_UNDEFINED;
 
  TracerParticleCreationSpacing = FLOAT_UNDEFINED;
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    TracerParticleCreationLeftEdge[dim] = FLOAT_UNDEFINED;
    TracerParticleCreationRightEdge[dim] = FLOAT_UNDEFINED;
  }
 
  // Read until out of lines
 
  rewind(fptr);
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    // Read tracer particle parameters
 
    sscanf(line, "TracerParticleCreation = %"ISYM, &MetaData.CycleNumber);
    sscanf(line, "TracerParticleCreationSpacing = %"PSYM,
	   &TracerParticleCreationSpacing);
    sscanf(line, "TracerParticleCreationLeftEdge = %"PSYM" %"PSYM" %"PSYM,
		  TracerParticleCreationLeftEdge,
		  TracerParticleCreationLeftEdge+1,
		  TracerParticleCreationLeftEdge+2);
    sscanf(line, "TracerParticleCreationRightEdge = %"PSYM" %"PSYM" %"PSYM,
		  TracerParticleCreationRightEdge,
		  TracerParticleCreationRightEdge+1,
		  TracerParticleCreationRightEdge+2);
 
  }
 
/*
  fprintf(stderr, "TracerParticleCreation = %"ISYM"\n", MetaData.CycleNumber);
  fprintf(stderr, "TracerParticleCreationSpacing = %"PSYM"\n", TracerParticleCreationSpacing);
  fprintf(stderr, "TracerParticleCreationLeftEdge = %"PSYM" %"PSYM" %"PSYM"\n",
                  TracerParticleCreationLeftEdge[0],
                  TracerParticleCreationLeftEdge[1],
                  TracerParticleCreationLeftEdge[2]);
  fprintf(stderr, "TracerParticleCreationRightEdge = %"PSYM" %"PSYM" %"PSYM"\n",
                  TracerParticleCreationRightEdge[0],
                  TracerParticleCreationRightEdge[1],
                  TracerParticleCreationRightEdge[2]);
*/
 
  // If spacing is non-zero, then create particles
 
  // This is possible here if ||rgio is OFF
  // Otherwise, the tracer particles should be added after
  // the regular particles have been read in.  Each of the
  // top level grids then has to check for overlap of this
  // region and create local particles and contribute to the
  // global sum for the total particle number
 
  if (ParallelRootGridIO == FALSE) {
    if (TracerParticleCreationSpacing > 0) {
      if (TopGrid.GridData->TracerParticleCreateParticles(
                                       TracerParticleCreationLeftEdge,
                                       TracerParticleCreationRightEdge,
                                       TracerParticleCreationSpacing,
                                       MetaData.NumberOfParticles) == FAIL) {
        ENZO_FAIL("Error in grid->TracerParticleCreateParticles.\n");

      }
    }
  }
 
  rewind(fptr);
 
  return SUCCESS;
}
