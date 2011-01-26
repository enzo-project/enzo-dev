/***********************************************************************
/
/  GRID CLASS (DEPOSIT LOCAL PARTICLE POSITIONS ONTO THE SPECIFIED FIELD)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  May, 2009 by John Wise
/                Only consider particles local on the processor.  This
/                happens in FindSubgrids (RebuildHierarchy) when we
/                keep the particles on the subgrid's local processor.
/
/  PURPOSE:
/     This routine deposits the particle living in this grid into either
/       the GravitatingMassField or the GravitatingMassFieldParticles of
/       itself, depending on the value of DepositField.
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
 
/* This controls the maximum particle mass which will be deposited in
   the MASS_FLAGGING_FIELD.  Only set in Grid_SetFlaggingField. */
 
extern float DepositParticleMaximumParticleMass;
 
int grid::DepositParticlePositionsLocal(FLOAT DepositTime, int DepositField)
{
 
  /* Declarations. */
 
  int dim, i;
  float MassFactor = 1.0, *ParticleMassTemp, *ParticleMassPointer;
 
  /* If there are no particles, don't deposit anything. */
 
  if (NumberOfParticles == 0)
    return SUCCESS;
 
  /* If the Target is this grid and the DepositField is MassFlaggingField,
     then multiply the Particle density by the volume to get the mass. */
 
  if (DepositField == MASS_FLAGGING_FIELD || 
      DepositField == PARTICLE_MASS_FLAGGING_FIELD)
    for (dim = 0; dim < GridRank; dim++)
      MassFactor *= CellWidth[dim][0];
 
  /* If required, Change the mass of particles in this grid. */
 
  if (MassFactor != 1.0) {
    ParticleMassTemp = new float[NumberOfParticles];
    for (i = 0; i < NumberOfParticles; i++)
      ParticleMassTemp[i] = ParticleMass[i]*MassFactor;
    ParticleMassPointer = ParticleMassTemp;
  } else
    ParticleMassPointer = ParticleMass;
 
  /* If the target field is MASS_FLAGGING_FIELD, then set masses of
     particles which are too large to zero (to prevent run-away refinement). */
 
  if ((DepositField == MASS_FLAGGING_FIELD ||
       DepositField == PARTICLE_MASS_FLAGGING_FIELD) &&
      DepositParticleMaximumParticleMass > 0 && MassFactor != 1.0)
    for (i = 0; i < NumberOfParticles; i++)
      ParticleMassPointer[i] = min(DepositParticleMaximumParticleMass,
				   ParticleMassPointer[i]);
 
  /* Compute difference between current time and DepositTime. */
 
  float TimeDifference = DepositTime - Time;
 
  /* Move particles to positions at Time + TimeDifference. */
 
  this->UpdateParticlePosition(TimeDifference);
 
  /* Deposit particles. */
 
//  fprintf(stderr, "----DPP Call this->DepositPositions with NP = %"ISYM"\n", NumberOfParticles);
 
  if (this->DepositPositions(ParticlePosition, ParticleMassPointer,
			     NumberOfParticles, DepositField) == FAIL) {
    ENZO_FAIL("Error in grid->DepositPositions\n");
  }
 
  /* If necessary, delete the particle mass temporary. */
 
  if (MassFactor != 1.0)

    delete [] ParticleMassTemp;
 
  /* Return particles to positions at Time. */
 
  this->UpdateParticlePosition(-TimeDifference);
 
  return SUCCESS;
}
