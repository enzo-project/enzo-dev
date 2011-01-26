/***********************************************************************
/
/  ADJUST MustRefineParticlesRefineToLevel PARAMETER IF REQUESTED
/
/  written by: Ji-hoon Kim
/  date:       April, 2010
/  modified1:
/
/  PURPOSE: If requested (MustRefineParticlesRefineToLevelAutoAdjust > 0),
/           adjust MustRefineParticlesRefineToLevel parameter 
/           When nonzero integer value (in pc) is given for this parameter, 
/           AdjustMustRefineParticlesRefineToLevel will be called at every 
/           rootgrid timestep to automatically change MustRefineParticlesRefineToLevel
/           so that the cell size at that level matches the corresponding 
/           physical length of MustRefineParticlesRefineToLevelAutoAdjust.
/
/           For example, if you have an AMR run with 16 comoving Mpc boxsize 
/           (topgrid = 128^3), and MustRefineParticlesRefineToLevelAutoAdjust 
/           = 120 (pc), then
/
/           - at z=7, MustRefineParticlesRefineToLevel 
/                     = floor [ ln(16 Mpc/(1+7)/128/120 pc) / ln(2) ] = 7
/           - at z=0, MustRefineParticlesRefineToLevel 
/                     = floor [ ln(16 Mpc/(1+0)/128/120 pc) / ln(2) ] = 10
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
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
#include "Star.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int AdjustMustRefineParticlesRefineToLevel(TopGridData *MetaData, int EL_level)
{

  if (MustRefineParticlesRefineToLevelAutoAdjust == FALSE || EL_level != 0)
    return SUCCESS;

  if (!ComovingCoordinates)
    return SUCCESS;

  if (!(MetaData->TopGridRank == 3 &&
	MetaData->TopGridDims[0] == MetaData->TopGridDims[1] &&
	MetaData->TopGridDims[1] == MetaData->TopGridDims[2])) {
    ENZO_FAIL("Your setup of assymmetrical topgrid is never tested.\n");
  }
      
  const double pc = 3.086e18, ln2 = 0.69314718;
  int MustRefineParticlesRefineToLevel_prev = MustRefineParticlesRefineToLevel;

  /* Compute Units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
    VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, MetaData->Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Compute new MustRefineParticlesRefineToLevel */

  //   Note:  LengthUnits / TopGridDims[0] / (2^MustRefineParticlesRefineToLevel) 
  //          = MustRefineParticlesRefineToLevelAutoAdjust

  MustRefineParticlesRefineToLevel = 
    (int) floor( log( LengthUnits / MetaData->TopGridDims[0] 
		     / MustRefineParticlesRefineToLevelAutoAdjust / pc) / ln2);

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    printf("AdjustMustRefineParticlesRefineToLevel: Changed MustRefineParticlesRefineToLevel from %"ISYM" to %"ISYM"\n", 
	   MustRefineParticlesRefineToLevel_prev, MustRefineParticlesRefineToLevel);
//  printf("AdjustMustRefineParticlesRefineToLevel: LengthUnits = %g\n", LengthUnits);
  }

  return SUCCESS;

}
