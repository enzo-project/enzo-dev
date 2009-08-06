/***********************************************************************
/
/  RETURNS A LIST OF MASSIVE PARTICLES
/
/  written by: John Wise
/  date:       April, 2009
/  modified1:
/
/  NOTES:  We look for massive particles in the current refine region. 
/
************************************************************************/
#include <stdlib.h>
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

int grid::FindMassiveParticles(float min_mass, int level, FLOAT *pos[], int &npart,
			       int CountOnly)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0)
    return SUCCESS;

  int i, dim;
  float AdjustedMass;
  FLOAT r;
  bool inside;

  // Factor of 1.5 to avoid any round-off errors
  AdjustedMass = 1.5 * min_mass * pow(8.0, level);

  if (CountOnly == TRUE) {

    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleType[i] == PARTICLE_TYPE_DARK_MATTER &&
	  ParticleMass[i] > AdjustedMass) {
	inside = true;
	for (dim = 0; dim < MAX_DIMENSION; dim++) {
	  inside &= (ParticlePosition[dim][i] >= RefineRegionLeftEdge[dim] &&
		     ParticlePosition[dim][i] <= RefineRegionRightEdge[dim]);
	  if (!inside) break;
	}
	if (inside) npart++;
      } // ENDIF low-resolution DM particle
  } // ENDIF CountOnly
  else {

    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleType[i] == PARTICLE_TYPE_DARK_MATTER &&
	  ParticleMass[i] > AdjustedMass) {
	inside = true;
	for (dim = 0; dim < MAX_DIMENSION; dim++) {
	  inside &= (ParticlePosition[dim][i] >= RefineRegionLeftEdge[dim] &&
		     ParticlePosition[dim][i] <= RefineRegionRightEdge[dim]);
	  if (!inside) break;
	}
	if (inside) {
	  for (dim = 0; dim < MAX_DIMENSION; dim++)
	    pos[dim][npart] = ParticlePosition[dim][i];
	  npart++;
	} // ENDIF inside
      } // ENDIF low-resolution DM particle
  } // ENDELSE CountOnly

  return SUCCESS;

}
