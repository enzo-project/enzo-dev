/***********************************************************************
/
/  GRID: RETURN PARTICLE INFORMATIONS
/
/  written by: Peng Wang
/  date:       January, 2009
/  modified1:
/
/  PURPOSE: return particle informations in structure array
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


int grid::ReturnParticleEntry(ParticleEntry *List)
{

  for (int i = 0; i < NumberOfParticles; i++) {
    List->Number = ParticleNumber[i];
    List->Type = ParticleType[i];
    List->Mass = ParticleMass[i]*pow(CellWidth[0][0],3);
    for (int dim = 0; dim < GridRank; dim++) {
      List->Position[dim] = ParticlePosition[dim][i];
      List->Velocity[dim] = ParticleVelocity[dim][i];
    }
    for (int n = 0; n < NumberOfParticleAttributes; n++) 
      List->Attribute[n] = ParticleAttribute[n][i];
    List++;
  }

  return NumberOfParticles;

}

void grid::RemoveMergedParticles(ParticleEntry *List, const int &Size, int *Flag)
{

  for (int i = 0; i < NumberOfParticles; i++) {

    /* for every particle, check the global list whether it will be merged */

    int id = ParticleNumber[i];
    for (int j = 0; j < Size; j++) {
      if (id == List[j].Number) {
	if (Flag[j] >= 0) {
	  if (ParticleType[i] ==  PARTICLE_TYPE_MUST_REFINE) 
	    ParticleMass[i] = FLOAT_UNDEFINED;
	  else 
	    ParticleMass[i] = tiny_number;	  
	  ParticleType[i] = PARTICLE_TYPE_DARK_MATTER;
	}
      }
      
    }
    
  }
}
