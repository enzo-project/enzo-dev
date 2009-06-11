/***********************************************************************
/
/  GRID CLASS (CHECK THE NUMERICAL GRAVITY FORCE AGAINST ANALYTIC)
/
/  written by: Greg Bryan
/  date:       July, 1995
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::TestGravityCheckResults(FILE *fptr, grid *TopGrid)
{
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  int dim, i;
  float dist[MAX_DIMENSION], Middle[MAX_DIMENSION], Width[MAX_DIMENSION];
  float r, fanalytic, fradial, ftang, pi = 3.14159;
 
  /* Set constants. */
 
  for (dim = 0; dim < GridRank; dim++) {
    Middle[dim] = 0.5*(DomainLeftEdge[dim] + DomainRightEdge[dim]);
    Width[dim]  =      DomainRightEdge[dim] - DomainLeftEdge[dim] ;
  }
 
  /* Set top grid cell size. */
 
  float TopGridCellWidth = TopGrid->CellWidth[0][0];
  dtFixed = max(dtFixed, TopGrid->dtFixed);
 
  /* Loop over particles, computing radial distance from the center and
     comparing the analyic force to the compute acceleration (determined
     by assuming the velocity is due soley to the analytic acceleration). */
 
  for (i = 1; i < NumberOfParticles; i++) {
 
    /* Compute distance. */
 
    r       = 0.0;
    fradial = 0.0;
    for (dim = 0; dim < GridRank; dim++) {
 
      dist[dim] = ParticlePosition[dim][i] - Middle[dim];
      if (fabs(dist[dim]) > Middle[dim] &&
	  GravityBoundaryType != TopGridIsolated)
	dist[dim] = -(Width[dim] - fabs(dist[dim]))*sign(dist[dim]);
 
      /* Compute distance. */
 
      r += dist[dim]*dist[dim];
 
      /* Compute radial component. */
 
      fradial += ParticleVelocity[dim][i]/dtFixed * dist[dim];
 
    }
 
    /* Normalize radial force components. */
 
    r = sqrt(r);
    fradial /= r;
 
    /* Compute tangential component. */
 
    ftang = 0.0;
    for (dim = 0; dim < GridRank; dim++)
      ftang += POW(ParticleVelocity[dim][i]/dtFixed, float(2.0));
    ftang = sqrt(max(ftang - fradial*fradial, 0.0));
 
    /* Compute analytic acceleration. */
 
    if (GridRank == 1)
      fanalytic = 2.0*pi*TopGridCellWidth;
    if (GridRank == 2)
      fanalytic = 2.0/r*TopGridCellWidth*TopGridCellWidth;
    if (GridRank == 3)
      fanalytic = 1.0/(r*r)*POW(TopGridCellWidth, 3);
 
    /* Output results. */
 
    fprintf(fptr, "%"FSYM"  %e   %e   %e\n", r/TopGridCellWidth, ftang,
	    -fradial, fanalytic);
 
  } // end loop over particles.
 
 
  /* If requested, set up the baryon field. */
 
  if (NumberOfBaryonFields > 0) {
 
    /* INCOMPLETE */
 
  }
 
  return SUCCESS;
}
