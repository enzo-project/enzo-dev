/***********************************************************************
/
/  CALCULATE FEEDBACK SPHERE PARAMETERS
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "StarParticleData.h"

void Star::CalculateFeedbackParameters(float &Radius, 
				       float RootCellWidth,
				       float SNe_dt, double &EjectaDensity,
				       double &EjectaThermalEnergy,
				       double &EjectaMetalDensity,
				       float DensityUnits, float LengthUnits, 
				       float TemperatureUnits, float TimeUnits,
				       float VelocityUnits)
{

  // Parameters for the Stroemgen sphere in Whalen et al. (2004)
  const float	BirthRadius	  = 50;		// pc
  const float	WhalenTemperature = 20000;	// K
  const float	WhalenDensity	  = 1;	// cm^-3
  const float	WhalenMaxVelocity = 35;		// km/s

  const double pc = 3.086e18, Msun = 1.989e33, pMass = 1.673e-24, 
    gravConst = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13;

  float StarLevelCellWidth;
  double EjectaVolume, SNEnergy, HeliumCoreMass, Delta_SF;

  Radius = 0.0;
  EjectaDensity = 0.0;
  EjectaThermalEnergy = 0.0;
  EjectaMetalDensity = 0.0;
  StarLevelCellWidth = RootCellWidth / powf(float(RefineBy), float(this->level));

  switch (this->FeedbackFlag) {
  case SUPERNOVA:  // pair-instability SNe
    Radius = PopIIISupernovaRadius * pc / LengthUnits;
    Radius = max(Radius, 3.5*StarLevelCellWidth);
    EjectaVolume = 4.0/3.0 * 3.14159 * pow(PopIIISupernovaRadius*pc, 3);
    EjectaDensity = Mass * Msun / EjectaVolume / DensityUnits;
    HeliumCoreMass = (13./24.) * (Mass - 20);
    SNEnergy = (5.0 + 1.304 * (HeliumCoreMass - 64)) * 1e51;
    EjectaMetalDensity = HeliumCoreMass * Msun / EjectaVolume / 
      DensityUnits;
    EjectaThermalEnergy = SNEnergy / (Mass * Msun) / VelocityUnits /
      VelocityUnits;

    // Exaggerate influence radius because the blastwave will enter
    // into some of the surrounding parent grids within the next
    // timestep if we inject the energy into a small radius.
    Radius *= 8.0;
    break;

  case STROEMGREN:
    Radius = BirthRadius * pc / LengthUnits;
    Radius = max(Radius, 1.5*StarLevelCellWidth);
    EjectaDensity = WhalenDensity * pMass / DensityUnits;
    EjectaThermalEnergy =
      WhalenTemperature / (TemperatureUnits * (Gamma-1.0) * 0.6);
    break;

  case FORMATION:
    Radius = 0;
    break;

  case CONT_SUPERNOVA:
    // inject energy into a sphere
    Radius = StarClusterSNRadius * pc / LengthUnits;
    Radius = max(Radius, 2*StarLevelCellWidth);

    // Release SNe energy constantly over 16 Myr (t = 4-20 Myr).
    Delta_SF = Mass * SNe_dt * TimeUnits / (16.0*Myr);
    EjectaVolume = 4.0/3.0 * 3.14159 * pow(Radius*LengthUnits, 3);
    EjectaDensity = Delta_SF * Msun / EjectaVolume / DensityUnits;
    EjectaMetalDensity = EjectaDensity * StarMetalYield;
    EjectaThermalEnergy = StarClusterSNEnergy / Msun /
      (VelocityUnits * VelocityUnits);
    break;

  } // ENDSWITCH FeedbackFlag
  
  return;
}
