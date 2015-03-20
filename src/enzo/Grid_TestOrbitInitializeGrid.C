/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GRAVITY ORBIT TEST)
/
/  written by: Greg Bryan
/  date:       August, 2006
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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::TestOrbitInitializeGrid(int NumberOfTestParticles,
				  FLOAT TestRadius,
				  float CentralMass,
				  float TestMass,
				  int UseBaryons)
{
  /* declarations */

  int dim, i;

  NumberOfParticles = NumberOfTestParticles + 1;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  /* Get Units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1;
  double MassUnits = 1;
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Set number of particles for this grid and allocate space. */

  this->AllocateNewParticles(NumberOfParticles);

  /* Set particle IDs and types */

  for (i = 0; i < NumberOfParticles; i++) {
    ParticleNumber[i] = i;

    /* if MRPs are activated, set the particle type for test particles to MUST_REFINE */
    if (MustRefineParticlesRefineToLevel > 0 && i > 0) {
      ParticleType[i] = PARTICLE_TYPE_MUST_REFINE;
    }else{
      ParticleType[i] = PARTICLE_TYPE_DARK_MATTER;
    }
  }

  /* Set central particle. */

  for (dim = 0; dim < GridRank; dim++) {
    ParticlePosition[dim][0] = 0.5*(DomainLeftEdge[dim]+DomainRightEdge[dim]);
    ParticleVelocity[dim][0] = 0;
  }
  ParticleMass[0] = CentralMass/pow(CellWidth[0][0], 3);  // what is actually stored in ParticleMass is density so divide by volume

  /* Compute the circular velocity, first in cgs units, and then divide
     by code units to get the correct velocity. */

  double BigGee = GravitationalConstant/(4.0*M_PI);  // big G is the constant/4pi
  double MassCGS = CentralMass*MassUnits;

/* JRT 09/13/06  replace following line by */

/* float circular_velocity = sqrt(BigGee*MassCGS/(TestRadius*LengthUnits)) / VelocityUnits; */
  double MassTotal = (CentralMass+TestMass) ;
  float circular_velocity = sqrt(float(BigGee*MassTotal)*MassUnits/(float(TestRadius)*LengthUnits)) / VelocityUnits;
  float CentralVelocity = circular_velocity  * TestMass / float(MassTotal) ;
  float TestVelocity = circular_velocity  * CentralMass / float(MassTotal) ;
/* JRT 09/13/06  replace following line by */

  /* Create test particle(s). */

  if (NumberOfTestParticles != 1) {
    ENZO_FAIL("Only 1 test particle may be created (for now).\n");

  }

  /* This is an orbit in the x-y plane. */

  ParticlePosition[0][1] = ParticlePosition[0][0] - TestRadius;
  ParticlePosition[1][1] = ParticlePosition[1][0];
  ParticlePosition[2][1] = ParticlePosition[2][0];

  ParticleVelocity[0][1] = 0;

/* JRT 09/13/06 */
/*ParticleVelocity[1][1] = circular_velocity; */
  ParticleVelocity[1][1] = TestVelocity;
/* JRT 09/13/06 */
  ParticleVelocity[2][1] = 0;

  ParticleMass[1] = TestMass/pow(CellWidth[0][0], 3);  // what is actually stored in ParticleMass is density so divide by volume

/* JRT 09/13/06  */
  
  ParticleVelocity[1][0] = - CentralVelocity;
/* JRT 09/13/06  */


  printf("The particle masses are:\n");
  printf("   (central)   %e\n",ParticleMass[0]);
  printf("   (test)      %e\n",ParticleMass[1]);

  printf("The particle positions are:\n");
  printf("  (central)   %e %e %e\n",ParticlePosition[0][0],ParticlePosition[1][0],ParticlePosition[2][0] );
  printf("  (test)      %e %e %e\n",ParticlePosition[0][1],ParticlePosition[1][1],ParticlePosition[2][1] );

  printf("The particle velocities are:\n");
  printf("  (central)   %e %e %e\n",ParticleVelocity[0][0],ParticleVelocity[1][0],ParticleVelocity[2][0] );
  printf("  (test)      %e %e %e\n",ParticleVelocity[0][1],ParticleVelocity[1][1],ParticleVelocity[2][1] );

  if(UseBaryons)
    printf("\n    ** Baryon fields have been turned on! **\n");

  fflush(stdout);


  return SUCCESS;
}

