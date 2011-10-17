/***********************************************************************
/
/  COMPUTE DEDNER WAVE SPEEDS
/
/  written by: Tom Abel
/  date:       October 2009
/
/  ======================================================================= 
/  This routine computes the wave speeds used for the Dedner formalism.
/
/       Reference: e.g.  Matsumoto, PASJ, 2007, 59, 905 
/
/  Output: C_h and C_p global variables defined in global_data.h
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
#include "TopGridData.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "../hydro_rk/tools.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int ComputeDednerWaveSpeeds(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], 
			    int level, FLOAT dt0)
{
  /* Count the number of grids on this level. */
 
  if (HydroMethod != MHD_RK) 
    return SUCCESS;

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, MagneticUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, 1.0);
  MassUnits = DensityUnits*pow(LengthUnits,3);

  int lmax;
  LevelHierarchyEntry *Temp;
  // Determine current maximum level
  for (lmax = MAX_DEPTH_OF_HIERARCHY-1; lmax >= 0; lmax--) {
    Temp = LevelArray[lmax];
    if (Temp != NULL) 
      break;
  }

  //      lmax = 0; // <- Pengs version had lmax = 6

  // using this for a cosmology run ... 
  lmax = MaximumRefinementLevel;
  FLOAT dx0, dy0, dz0, h_min;
  
  dx0 = (DomainRightEdge[0] - DomainLeftEdge[0]) / MetaData->TopGridDims[0];
  dy0 = (MetaData->TopGridRank > 1) ? 
    (DomainRightEdge[1] - DomainLeftEdge[1]) / MetaData->TopGridDims[1] : 1e8;
  dz0 = (MetaData->TopGridRank > 2) ? 
    (DomainRightEdge[2] - DomainLeftEdge[2]) / MetaData->TopGridDims[2] : 1e8;
  h_min = my_MIN(dx0, dy0, dz0);
  h_min /= pow(RefineBy, lmax);
  C_h = 0.3*MetaData->CourantSafetyNumber*(h_min/dt0);
  //  C_h = min( C_h, 1e6/VelocityUnits); // never faster than __ cm/s (for very small dt0 a problems)
  if (EOSType == 3)  // for isothermal runs just use the constant sound speed
    C_h = EOSSoundSpeed;
  if (EOSType == 4 || EOSType == 5)  // for isothermal runs just use the constant sound speed
    C_h = EOSSoundSpeed;

  C_p = sqrt(0.18*DivBDampingLength*C_h);

  return SUCCESS;
}
 
  
