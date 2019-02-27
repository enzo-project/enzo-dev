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
#include "phys_constants.h"

void Star::MultiplyAccretionRate(float &RecalibrateAccretingMassRatio)
{

  float mdot, mdot_Edd;

  // Multiplication
  mdot = accretion_rate[0] * RecalibrateAccretingMassRatio; 

  // Below is exactly the same as in Star_CalculateMassAccretion
  mdot_Edd = 4.0 * PI * GravConst * this->Mass * mh /
    max(MBHFeedbackRadiativeEfficiency, 0.1) / sigma_thompson / clight; 

  accretion_rate[0] = min(mdot, mdot_Edd);

//  fprintf(stderr, "RecalibrateAccretingMassRatio = %g\n", 
//	  RecalibrateAccretingMassRatio);  

  return;
}




