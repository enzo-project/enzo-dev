/***********************************************************************
/
/  SUBTRACT ACCRETED MASS FROM THE GRID (AND CHANGES TO MOMENTUM)
/
/  written by: Ji-hoon Kim
/  date:       January, 2010
/  modified1: 
/
/ PURPOSE: To subtract mass from the grid, we must consider multiple grids
/          since sometimes the accretion radius often exceeds the grid
/          boundaries.  
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int StarParticleSubtractAccretedMass(TopGridData *MetaData, 
				     LevelHierarchyEntry *LevelArray[], int level, 
				     Star *&AllStars)
{

  const double pc = 3.086e18, Msun = 1.989e33, pMass = 1.673e-24, 
    gravConst = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13;

  Star *cstar;
  int i, l, dim, temp_int, SkipMassRemoval, SphereContained,
      SphereContainedNextLevel, dummy;
  float influenceRadius, RootCellWidth, SNe_dt, mdot;
  float dtForThisStar;
  double EjectaThermalEnergy, EjectaDensity, EjectaMetalDensity, EjectaVolume;
  FLOAT Time;
  LevelHierarchyEntry *Temp;

  if (AllStars == NULL)
    return SUCCESS;

  /* Get time and SNe timestep */

  Temp = LevelArray[level];
  Time = Temp->GridData->ReturnTime();
  if (LastSupernovaTime < 0)
    SNe_dt = 0.0;
  else
    SNe_dt = Time - LastSupernovaTime;
  LastSupernovaTime = Time;
  RootCellWidth = 1.0 / MetaData->TopGridDims[0];

  /* Set the units. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

//    if (debug)
//      cstar->PrintInfo();

    if (cstar->ReturnType() != BlackHole)
      continue;

    /* Now let us do the job! */

    switch (cstar->ReturnType()) {
    case BlackHole:

      /* If BlackHole, subtract mass from a single cell (old way) so as 
	 not to bother John or others at the moment! */ 

      if (cstar->SubtractAccretedMassFromCell() == FAIL) {
	fprintf(stderr, "Error in star::SubtractAccretedMass.\n");
	ENZO_FAIL("");
      }

      break;

    } // ENDSWITCH ReturnType()

  } // ENDFOR stars

  return SUCCESS;

}
