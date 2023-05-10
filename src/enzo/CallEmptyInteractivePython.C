/***********************************************************************
/
/  CALL LIBYT AT FIXED TIMESTEPS
/
/  written by: Shin-Rong Tsai, Matthew Turk
/  date:       April, 2023
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#ifdef USE_LIBYT
#include "libyt.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "TopGridData.h"

void CommunicationBarrier();

int CallEmptyInteractivePython(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level, int from_topgrid)
{
#ifndef USE_LIBYT
    return SUCCESS;
#else

  if (yt_run_InteractiveMode("LIBYT_STOP") != YT_SUCCESS) {
      printf("yt_interactive_mode failed!\n");
      return 1;
  }

  CommunicationBarrier();
  return SUCCESS;
#endif
}
