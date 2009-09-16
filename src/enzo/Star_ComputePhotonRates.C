/***********************************************************************
/
/  COMPUTE STELLAR PHOTON EMISSION RATES
/
/  written by: John Wise
/  date:       November, 2005
/  modified1:
/
/  ---------- SPECIES --------
/  0 : HI
/  1 : HeI
/  2 : HeII
/  3 : Lyman-Werner (H2)
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "StarParticleData.h"

int Star::ComputePhotonRates(float E[], double Q[])
{

  float x, x2, _mass, EnergyFractionLW, MeanEnergy, XrayLuminosityFraction;
  x = log10(this->Mass);
  x2 = x*x;

  switch(this->type) {

    /* Luminosities from Schaerer (2002) */

  case PopIII:
    E[0] = 28.0;
    E[1] = 30.0;
    E[2] = 58.0;
    E[3] = 12.8;
    _mass = max(min(this->Mass, 500), 5);
    if (_mass > 9 && _mass <= 500) {
      Q[0] = pow(10.0, 43.61 + 4.9*x   - 0.83*x2);
      Q[1] = pow(10.0, 42.51 + 5.69*x  - 1.01*x2);
      Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
      Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
    } else if (_mass > 5 && _mass <= 9) {
      Q[0] = pow(10.0, 39.29 + 8.55*x);
      Q[1] = pow(10.0, 29.24 + 18.49*x);
      Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
      Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
    } // ENDELSE
    break;

    /* Average energy from Schaerer (2003) */

  case PopII:
    EnergyFractionLW = 0.01;
    E[0] = 21.5; // eV (good for a standard, low-Z IMF)
    E[1] = 0.0;
    E[2] = 0.0;
    E[3] = 12.8;
    Q[0] = StarClusterIonizingLuminosity * this->Mass;
    Q[1] = 0.0;
    Q[2] = 0.0;
    Q[3] = EnergyFractionLW * Q[0];
    break;

    /* Approximation to the multi-color disk and power law of an
       accreting BH (Kuhlen & Madau 2004; Alvarez et al. 2009) */

  case BlackHole:
    XrayLuminosityFraction = 0.43;
    EnergyFractionLW = 1.51e-3;
    MeanEnergy = 93.0;  // eV
    E[0] = 460.0;
    E[1] = 0.0;
    E[2] = 0.0;
    E[3] = 12.8;
    Q[0] = 3.54e58 * PopIIIBHLuminosityEfficiency * XrayLuminosityFraction *
      this->DeltaMass / E[0];
    Q[1] = 0.0;
    Q[2] = 0.0;
    Q[3] = EnergyFractionLW * (E[0]/MeanEnergy) * Q[0];
    break;

    /* Approximation to the multi-color disk and power law of an
       accreting massive BH */

  case MBH:
    XrayLuminosityFraction = 0.43;
    EnergyFractionLW = 1.51e-3;
    MeanEnergy = 93.0;  // eV
    E[0] = 460.0;
    E[1] = 0.0;
    E[2] = 0.0;
    E[3] = 12.8;
    Q[0] = 3.54e58 * MBHFeedbackRadiativeEfficiency * XrayLuminosityFraction *
      this->DeltaMass / E[0];
    Q[1] = 0.0;
    Q[2] = 0.0;
    Q[3] = EnergyFractionLW * (E[0]/MeanEnergy) * Q[0];

    //fprintf(stdout, "this->DeltaMass = %g, Q[0]=%g\n", this->DeltaMass, Q[0]); 
    break;

  default:
    fprintf(stderr, "Star type = %"ISYM" not understood.\n", this->type);
    ENZO_FAIL("");
  } // ENDSWITCH

  return SUCCESS;
}
