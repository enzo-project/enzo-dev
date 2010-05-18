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

bool Star::ApplyFeedbackTrue(float dt)
{

  bool result = true;
  bool *rules;
  int i;

  const int NumberOfRules = 3;

  rules = new bool[NumberOfRules];

  rules[0] = !(FeedbackFlag == NO_FEEDBACK);
  rules[1] = !(FeedbackFlag == STROEMGREN && RadiativeTransfer);
  //rules[2] = !(FeedbackFlag == CONT_SUPERNOVA && dt <= 0);
  rules[2] = true;

  for (i = 0; i < NumberOfRules; i++)
    result &= rules[i];

  delete [] rules;
  
  return result;
}
