/***********************************************************************
/
/  CALCULATE MASS LOSS FROM SUPERNOVAE
/
/  written by: John Wise
/  date:       January, 2014
/  modified1: 
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
#include "CosmologyParameters.h"

#define LIFETIME_IN_TDYN 12.0

float Star::CalculateMassLoss(const float dt)
{

  float MassLoss = 0.0;
  float M0, Mform, xv1, xv2, tdyn;

  if (this->type != NormalStar)
    return 0.0;

#ifdef TRANSFER
  tdyn = this->LifeTime / LIFETIME_IN_TDYN;
  xv1 = (PhotonTime - this->BirthTime) / tdyn;
  xv2 = (PhotonTime + dt - this->BirthTime) / tdyn;
  M0 = this->Mass / (1.0 - StarMassEjectionFraction * 
		     (1.0 - (1.0 + xv1) * exp(-xv1)));
  Mform = M0 * ((1.0 + xv1) * exp(-xv1) -
		       (1.0 + xv2) * exp(-xv2));
  Mform = max(min(Mform, this->Mass), 0.0);
  MassLoss = StarMassEjectionFraction * Mform;
#endif

  return MassLoss;
}
