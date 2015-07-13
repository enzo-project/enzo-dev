/***********************************************************************
/
/  GRID CLASS (UPDATE PARTICLE VELOCITY FROM ACCELERATIONS)
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
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
#define VELOCITY_METHOD3
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
 
int grid::UpdateParticleVelocity(float TimeStep)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if (NumberOfParticles == 0 || ParticleAcceleration[0] == NULL) return SUCCESS;
 
  FLOAT a = 1.0, dadt;
#if defined(VELOCITY_METHOD1) || defined(VELOCITY_METHOD2)
  float VelocityMidStep;
#endif
  int i;
 
  /* If using comoving coordinates, divide by a(t) first. */
 
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time + TimeStep, &a, &dadt)
	== FAIL) {
            ENZO_FAIL("Error in CsomologyComputeExpansionFactors.");
    }
 
  /* Loop over dimensions. */
 
  for (int dim = 0; dim < GridRank; dim++) {
 
    /* Error check. */
 
    if (ParticleAcceleration[dim] == NULL) {
            ENZO_FAIL("No ParticleAcceleration present.");
    }
 
    /* Update velocities.  */
 
    if (ComovingCoordinates) {
 
      FLOAT coef = 0.5*dadt/a*TimeStep;
      FLOAT coef1 = 1.0 - coef;
      FLOAT coef2 = 1.0 / (1.0 + coef);
 
      /* If using comoving coordinates, subtract the (time-centered)
	 drag-like term and add the acceleration. The acceleration has
	 already been divided by a(t). */
 
      for (i = 0; i < NumberOfParticles; i++) {
 
#ifdef VELOCITY_METHOD1
 
        /* i) partially time-centered. */
 
	VelocityMidStep = ParticleVelocity[dim][i] +
	                  ParticleAcceleration[dim][i]*0.5*TimeStep;
 
	ParticleVelocity[dim][i] +=
	  (-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;
 
#endif /* VELOCITY_METHOD1 */
 
#ifdef VELOCITY_METHOD2
 
        /* ii) partially backward. */
 
	VelocityMidStep = ParticleVelocity[dim][i] ;
 
	ParticleVelocity[dim][i] +=
	  (-VelocityMidStep*dadt/a + ParticleAcceleration[dim][i]) * TimeStep;
 
#endif /* VELOCITY_METHOD2 */
 
#ifdef VELOCITY_METHOD3
 
        /* iii) Semi-implicit way */
 
        ParticleVelocity[dim][i] = (coef1*ParticleVelocity[dim][i] +
                                    ParticleAcceleration[dim][i]*TimeStep)*coef2;

 
#endif /* VELOCITY_METHOD3 */
 
      }
    }
    else
 
      /* Otherwise, just add the acceleration. */
 
      for (i = 0; i < NumberOfParticles; i++)
	ParticleVelocity[dim][i] += ParticleAcceleration[dim][i] * TimeStep;
 
  }


  if (ProblemType == 29)
    for (i = 0; i < NumberOfParticles; i++)
      printf("id=%"PISYM"  %"PSYM" %"PSYM" %"PSYM"  %"ESYM" %"ESYM" %"ESYM" \n", ParticleNumber[i],
	     ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i],
             ParticleVelocity[0][i], ParticleVelocity[1][i], ParticleVelocity[2][i]);

 
  return SUCCESS;
}
