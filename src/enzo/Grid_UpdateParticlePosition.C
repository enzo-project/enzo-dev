/***********************************************************************
/
/  GRID CLASS (UPDATE PARTICLE POSITION FOR VELOCITY)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
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
 
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
 
int grid::UpdateParticlePosition(float TimeStep)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if (NumberOfParticles == 0) return SUCCESS;
 
  FLOAT a = 1.0, dadt;
  int i, dim;
 
//  if (debug)
//    printf("UpdateParticlePosition: moving %"ISYM" particles forward by %"FSYM".\n",
//	   NumberOfParticles, TimeStep);
 
  /* If using comoving coordinates, divide the acceleration by a(t) first.
     (We use abs(TimeStep) because this routine is occasionally used to
     move particles forward dt and then immediately reverse this (-dt);
     using abs(dt) keeps things consistent). */
 
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time + 0.5*fabs(TimeStep), &a, &dadt)
	== FAIL) {
      fprintf(stderr, "Error in CsomologyComputeExpansionFactors.\n");
      ENZO_FAIL("Error in: "__FILE__);
    }
 
  /* Loop over dimensions. */
 
  for (dim = 0; dim < GridRank; dim++) {
 
    /* Error check. */
 
    if (ParticleVelocity[dim] == NULL) {
      fprintf(stderr, "No ParticleVelocity present.\n");
      ENZO_FAIL("Error in: "__FILE__);
    }
 
    /* update velocities. */
 
    float Coefficient = TimeStep/a;
    for (i = 0; i < NumberOfParticles; i++)
      ParticlePosition[dim][i] += Coefficient*ParticleVelocity[dim][i];
 
    /* wrap particle positions for periodic case.
       (now done in CommunicationTransferParticles) */
 
#ifdef UNUSED
    FLOAT Width = DomainRightEdge[dim] - DomainLeftEdge[dim];
    for (i = 0; i < NumberOfParticles; i++) {
      if (ParticlePosition[dim][i] > DomainRightEdge[dim])
	ParticlePosition[dim][i] -= Width;
      if (ParticlePosition[dim][i] < DomainLeftEdge[dim])
	ParticlePosition[dim][i] += Width;
    }
#endif /* UNUSED */
 
  }
 
  return SUCCESS;
}
