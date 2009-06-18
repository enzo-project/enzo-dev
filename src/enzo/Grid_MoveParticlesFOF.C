/***********************************************************************
/
/  GRID CLASS (MOVE ALL PARTICLES FROM SPECIFIED GRID TO FOF STRUCTURE)
/
/  written by: John Wise
/  date:       June, 2009
/  modified1:  
/
/  PURPOSE:
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

#include "FOF_allvars.h"
 
int grid::MoveParticlesFOF(int level, int GridNum, FOF_particle_data* &P, 
			   int &Index, FOFData &AllVars, float VelocityUnits, 
			   double MassUnits, int CopyDirection)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int i, j, dim, Index0;
  double VelocityConv, MassConv;

  VelocityConv = VelocityUnits / AllVars.UnitVelocity_in_cm_per_s;
  MassConv = MassUnits / pow(8.0, level) / AllVars.UnitMass_in_g;

  if (CopyDirection == COPY_OUT) {

    Index0 = Index;

    for (dim = 0; dim < GridRank; dim++) {
      Index = Index0;
      for (i = 0; i < NumberOfParticles; i++, Index++) {
	P[Index].Pos[dim] = AllVars.BoxSize * ParticlePosition[dim][i];
	P[Index].Vel[dim] = ParticleVelocity[dim][i] * VelocityConv;
      }
    }

    for (j = 0; j < NumberOfParticleAttributes; j++) {
      Index = Index0;
      for (i = 0; i < NumberOfParticles; i++, Index++)
	P[Index].Attr[j] = ParticleAttribute[j][i];
    }

    Index = Index0;
    for (i = 0; i < NumberOfParticles; i++) {

      P[Index].slab = (int) (ParticlePosition[0][i] * NumberOfProcessors);

      P[Index].Mass = ParticleMass[i] * MassConv;

      P[Index].PartID = ParticleNumber[i];
      P[Index].Type = ParticleType[i];
      P[Index].level = level;
      P[Index].GridID = GridNum;

      P[Index].Energy = 0.0;
      P[Index].Rho = 0.0;

    }

    this->NumberOfParticles = 0;
    this->DeleteParticles();

  } // ENDIF (COPY_OUT)

  else {

    i = NumberOfParticles;
    for (dim = 0; dim < GridRank; dim++) {
      ParticlePosition[dim][i] = P[Index].Pos[dim] / AllVars.BoxSize;
      ParticleVelocity[dim][i] = P[Index].Vel[dim] / VelocityConv;
    }

    for (j = 0; j < NumberOfParticleAttributes; j++)
      ParticleAttribute[j][i] = P[Index].Attr[j];

    ParticleMass[i] = P[Index].Mass / MassConv;
    ParticleType[i] = P[Index].Type;
    ParticleNumber[i] = P[Index].PartID;

    NumberOfParticles++;

  } // ENDELSE (COPY_IN)

  return SUCCESS;

}
