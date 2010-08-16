/***********************************************************************
/
/  CHECK FOR QUEUE RE-SUBMIT
/
/  written by: John Wise
/  date:       July, 2008
/  modified1:
/
/  PURPOSE: Stop job and resubmit to queue if the next topgrid timestep 
/           will exceed the wallclock limit (StopCPUTime).
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
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

double ReturnWallTime();
int CheckForResubmit(TopGridData &MetaData, int &Stop)
{

  if (MetaData.ResubmitOn == FALSE)
    return SUCCESS;

  char *cmd = new char[512];
  //double CurrentCPUTime = ReturnWallTime() - MetaData.StartCPUTime;

  if (MetaData.CPUTime + MetaData.LastCycleCPUTime > MetaData.StopCPUTime) {
    if (debug)
      printf("Next topgrid timestep will exceed StopCPUTime.  Stopping.\n"
	     "Executing resubmission script, %s\n", MetaData.ResubmitCommand);
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      if (MetaData.DataDumpDir != NULL)
	sprintf(cmd, "%s/%s %"ISYM" %s%"CYCLE_TAG_FORMAT""ISYM"/%s%"CYCLE_TAG_FORMAT""ISYM,
		MetaData.GlobalDir, MetaData.ResubmitCommand, NumberOfProcessors, 
		MetaData.DataDumpDir, MetaData.DataDumpNumber-1, 
		MetaData.DataDumpName, MetaData.DataDumpNumber-1);
      else
	sprintf(cmd, "%s/%s %"ISYM" %s%"CYCLE_TAG_FORMAT""ISYM, 
		MetaData.GlobalDir, MetaData.ResubmitCommand, NumberOfProcessors, 
		MetaData.DataDumpName, MetaData.DataDumpNumber-1);
      printf("command: %s\n", cmd);
      system(cmd);
    }
    Stop = TRUE;
  } // ENDIF

  delete [] cmd;

  return SUCCESS;
}
