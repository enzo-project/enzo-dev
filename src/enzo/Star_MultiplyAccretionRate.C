/***********************************************************************
/
/  CORRECT THE ACCRETION RATE BY MULTIPLYING BY A CONSTANT
/
/  written by: Ji-hoon Kim
/  date:       June, 2010
/
/  modified1: 
/         
/  PURPOSE:    In the case of MBHAccretingMassRatio = BONDI_ACCRETION_
/              CORRECT_NUMERICAL, recalibrate (correct) the accreting 
/              rate by multiplying RecalibrateAccretingMassRatio. Do 
/              remember that the new accretion rate should be bounded
/              by the Eddington limit.
/
/  CAUTION:    Here this->naccretions = 1 is assumed as in 
/              Star_CalculateMassAccretion.C  
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

void Star::MultiplyAccretionRate(float &RecalibrateAccretingMassRatio)
{

  const double Grav = 6.673e-8, k_b = 1.38e-16, m_h = 1.673e-24;
  const double Msun = 1.989e33, yr = 3.1557e7, sigma_T = 6.65e-25, c = 3.0e10;

  float mdot, mdot_Edd;

  // Multiplication
  mdot = accretion_rate[0] * RecalibrateAccretingMassRatio; 

  // Below is exactly the same as in Star_CalculateMassAccretion
  mdot_Edd = 4.0 * PI * Grav * this->Mass * m_h /
    max(MBHFeedbackRadiativeEfficiency, 0.1) / sigma_T / c; 

  accretion_rate[0] = min(mdot, mdot_Edd);

//  fprintf(stderr, "RecalibrateAccretingMassRatio = %g\n", 
//	  RecalibrateAccretingMassRatio);  

  return;
}




