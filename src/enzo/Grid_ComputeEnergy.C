/***********************************************************************
/
/  GRID CLASS (SUM THE ENERGY TERMS OVER THE GRID)
/
/  written by: Greg Bryan
/  date:       March, 1995
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
 
 
int grid::ComputeEnergy(float EnergySum[])
{
 
  /* declarations */
 
  int i, j, k, n, dim, index;
 
  float CellVolume = 1, KineticEnergy, mass;
  for (dim = 0; dim < GridRank; dim++)
    CellVolume *= CellWidth[dim][0];
 
  /* 0,1) baryon energy terms. */
 
  if (NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }
 
    /* Sum over mesh. */
 
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	index = (k*GridDimension[1] + j)*GridDimension[0] +
	  GridStartIndex[0];
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
 
	  mass = BaryonField[DensNum][index]*CellVolume;
 
	  /* Kinetic energy. */
 
	  KineticEnergy = BaryonField[Vel1Num][index]*
	                  BaryonField[Vel1Num][index];
	  if (GridRank > 0)
	    KineticEnergy += BaryonField[Vel2Num][index]*
	                     BaryonField[Vel2Num][index];
	  if (GridRank > 1)
	    KineticEnergy += BaryonField[Vel3Num][index]*
	                     BaryonField[Vel3Num][index];
 
	  KineticEnergy *= 0.5*mass;
 
	  EnergySum[0] += KineticEnergy;
 
	  /* Thermal energy. */
 
	  if (DualEnergyFormalism)
	    EnergySum[1] += BaryonField[GENum][index]*mass;
	  else
	    EnergySum[1] += max(BaryonField[TENum][index]*mass -
				KineticEnergy, 0);
 
	  /* Gravitational potential energy (the if statement is to
	     elliminate regions covered by other grids). */
 
#ifdef UNUSED
	  if (SelfGravity && ComputePotential)
//	    if (BaryonField[NumberOfBaryonFields][index] == 0)
	      EnergySum[3] += mass*AccelerationField[GridRank][index];
#endif /* UNUSED */
	
	  }
      }
 
  } // end: if (NumberOfBaryonFields > 0)
 
  /* 2) particle kinetic energy. */
 
  for (dim = 0; dim < GridRank; dim++)
    for (n = 0; n < NumberOfParticles; n++)
      EnergySum[2] += 0.5*ParticleVelocity[dim][n]*ParticleVelocity[dim][n]*
	ParticleMass[n]*CellVolume;
 
  /* 4) Particle gravitational potential energy. */
 
#ifdef UNUSED
  int FieldPosition[MAX_DIMENSION];
  if (SelfGravity && ComputePotential) {
  }
#endif /* UNUSED */
 
  // BUG?? PotentialSum never defined??
  if (ComputePotential)

    EnergySum[3] += PotentialSum;
 
  return SUCCESS;
}
 
