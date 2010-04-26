/***********************************************************************
/
/  SET MINIMUM REFINEMENT LEVEL FOR METALLICITY
/
/  written by: John Wise
/  date:       April, 2010
/  modified1:
/
/  PURPOSE: 
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
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
#include "CommunicationUtilities.h"

double CalculateBlastWaveRadius(double Mass, double n0, double Time);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int StarParticleSetRefinementLevel(Star *AllStars)
{

  if (PopIIISupernovaMustRefine == FALSE)
    return SUCCESS;

  // End metallicity refinement after 1% of the star's lifetime
  const float EndRefineAtTime = 1.01;
  const float AmbientDensity = 1.0;

  FLOAT GridTime;
  int RefinementLevel = -1;
  float factor, Diameter, DesiredResolution;
  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits;
  Star *cstar;

  for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

    // CurrentGrid is NULL if the grid is not on this processor.
    if (cstar->ReturnCurrentGrid() == NULL)
      continue;

    // Find Pop III stars that have just gone supernova (lifetime,
    // mass range, near-zero mass)

    if (cstar->ReturnType() == PopIII) {
      GridTime = cstar->ReturnCurrentGrid()->ReturnTime();
      factor = (GridTime - cstar->ReturnBirthtime()) / cstar->ReturnLifetime();
      if (factor > 1 && factor < EndRefineAtTime && 
	  cstar->ReturnMass() <= 1e-9) {

	GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, GridTime);

	Diameter = 2.0 *
	  float(CalculateBlastWaveRadius(cstar->ReturnFinalMass(), 
					 AmbientDensity,
					 (factor-1)*TimeUnits) 
		/ double(LengthUnits));

	// Given how many cells we want to resolve the diameter,
	// calculate the refinement level
	DesiredResolution = Diameter / PopIIISupernovaMustRefineResolution;

	RefinementLevel = 
	  nint(ceil(logf(TopGridDx[0]/DesiredResolution) / logf(RefineBy)));

	printf("Diameter = %g, factor = %f, Need_dx = %g (lvl %d), "
	       "MustRefineTo = %d\n",
	       Diameter, factor, DesiredResolution, cstar->ReturnLevel(),
	       RefinementLevel);
	RefinementLevel = min(max(RefinementLevel, 0), RefinementLevel);

      } // ENDIF recent SN
    } // ENDIF PopIII star

  } // ENDFOR stars

  MetallicityRefinementMinLevel = CommunicationMaxValue(RefinementLevel);

  return SUCCESS;

}
