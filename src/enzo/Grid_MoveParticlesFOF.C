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
 
int grid::MoveParticlesFOF(int level, FOF_particle_data* &P, 
			   int &Index, FOFData AllVars, float VelocityUnits, 
			   double MassUnits, int CopyDirection)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int i, j, dim, Index0;
  double VelocityConv, MassConv;
  double VelConvInv, MassConvInv;
  double BoxSizeInv;

  VelocityConv = VelocityUnits / AllVars.UnitVelocity_in_cm_per_s;
  MassConv = MassUnits / pow(8.0, level) / AllVars.UnitMass_in_g;

  VelConvInv = 1.0 / VelocityConv;
  MassConvInv = 1.0 / MassConv;
  BoxSizeInv = 1.0 / AllVars.BoxSize;

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
    for (i = 0; i < NumberOfParticles; i++, Index++) {

      P[Index].slab = (int) (ParticlePosition[0][i] * NumberOfProcessors);

      P[Index].Mass = ParticleMass[i] * MassConv;

      P[Index].PartID = ParticleNumber[i];
      P[Index].Type = ParticleType[i];
      //P[Index].level = level;
      //P[Index].GridID = GridNum;

      P[Index].Energy = 0.0;
      P[Index].Rho = 0.0;

    }

    this->DeleteParticles();

  } // ENDIF (COPY_OUT)

  else {

    /* When we create the smoothed DM fields, we copy particles across
       periodic boundaries.  Don't copy back these duplicated
       particles. */

    int npart = AllVars.Nlocal;//slab[MyProcessorNumber];
    for (i = 0; i < AllVars.Nlocal; i++)
      if (P[i].PartID < 0 || P[i].slab != MyProcessorNumber) npart--;

    /* Only move particles in this slab -- exclude the shadows. */

    this->AllocateNewParticles(npart);

    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < AllVars.Nlocal; i++)
	if (P[i].slab == MyProcessorNumber && P[i].PartID >= 0) {
	  ParticlePosition[dim][i] = P[i].Pos[dim] * BoxSizeInv;
	  ParticleVelocity[dim][i] = P[i].Vel[dim] * VelConvInv;
	}

    for (j = 0; j < NumberOfParticleAttributes; j++)
      for (i = 0; i < AllVars.Nlocal; i++)
	if (P[i].slab == MyProcessorNumber && P[i].PartID >= 0)
	  ParticleAttribute[j][i] = P[i].Attr[j];

    for (i = 0; i < AllVars.Nlocal; i++) 
      if (P[i].slab == MyProcessorNumber && P[i].PartID >= 0) {
	ParticleMass[i] = P[i].Mass * MassConvInv;
	ParticleType[i] = P[i].Type;
	ParticleNumber[i] = P[i].PartID;
      }

    NumberOfParticles = npart;

  } // ENDELSE (COPY_IN)

  return SUCCESS;

}
