/***********************************************************************
/
/  REFLECTS CHANGES IN STAR PARTICLE IN NORMAL PARTICLE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:  November, 2011 (JHW) -- porting to active particles
/
************************************************************************/
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "ActiveParticle.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

void ActiveParticleType::MirrorToParticle(void)
{

  if (CurrentGrid == NULL)
    return;

  const double Msun = 1.989e33;
  int i, dim, place = -1;
  float MassConversion;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, CurrentGrid->Time);

  double dx = LengthUnits * CurrentGrid->CellWidth[0][0];
  MassConversion = (float) (dx*dx*dx * double(DensityUnits) / Msun);

  // Find where this star particle is stored in main arrays
  for (i = 0; i < CurrentGrid->NumberOfParticles; i++) 
    if (CurrentGrid->ParticleNumber[i] == this->Identifier) {
      place = i;
      break;
    }

  if (place < 0) {
    printf("star::MTP: CurrentGrid->NumberOfParticles = %d, level = %d, "
	   "place =%d, Mass = %d, GridID = %d\n", 
	   CurrentGrid->NumberOfParticles, level, place, Mass, GridID); 
    printf("star::MTP: LeftEdge // RightEdge = %"PSYM" %"PSYM" %"PSYM
	   " // %"PSYM" %"PSYM" %"PSYM"\n",
	   CurrentGrid->GridLeftEdge[0], CurrentGrid->GridLeftEdge[1], 
	   CurrentGrid->GridLeftEdge[2], 
	   CurrentGrid->GridRightEdge[0], CurrentGrid->GridRightEdge[1], 
	   CurrentGrid->GridRightEdge[2]);
    this->PrintInfo();
    ENZO_FAIL("Cannot find matching particle for this star.");
  }
  //assert(place >= 0);

  // Change all particle data in favor of updated Star particle
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    CurrentGrid->ParticlePosition[dim][place] = this->pos[dim];
    CurrentGrid->ParticleVelocity[dim][place] = this->vel[dim];
  }
  CurrentGrid->ParticleMass[place] = (float)(this->Mass) / MassConversion;
  CurrentGrid->ParticleType[place] = this->type;
  CurrentGrid->ParticleAttribute[0][place] = this->BirthTime;
  CurrentGrid->ParticleAttribute[1][place] = this->DynamicalTime;
  CurrentGrid->ParticleAttribute[2][place] = this->Metallicity;
  
  return;
}
