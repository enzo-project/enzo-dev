/***********************************************************************
/
/  GRID CLASS (HANDLE CALLING AND SOLVING SHOCK ANALYSIS)
/
/  written by: Samuel Skillman
/  date:       July, 2009
/  modified1:
/
/  PURPOSE: Move logic for shock module selection here
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::ShocksHandler()
{
  if (!CRModel) return SUCCESS; 
  int shock_status;

  switch(ShockMethod){
  case 0:
    shock_status = this->FindShocks();
    break;
  case 1:
    shock_status = this->FindTempSplitShocks();
    break;
  case 2:
    shock_status = this->FindVelShocks();
    break;
  case 3:
    shock_status = this->FindVelSplitShocks();
    break;
  default:
    shock_status = this->FindShocks();
  }

  if(shock_status == FAIL){
    ENZO_FAIL("Error in grid->ShocksHandler.");
  }
  return SUCCESS;
}
