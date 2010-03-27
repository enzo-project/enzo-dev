/***********************************************************************
/
/  CHECK FOR OUTPUT
/
/  written by: Greg Bryan
/  date:       January, 1996
/  modified:   Robert Harkness
/  date:       January, 2007
/              Mods for group in-core i/o
/  date:       May, 2008
/              Remove Dan Reynold's iso_grav code
/
/  PURPOSE:
/    This routine checks a number of criteria for output and then calls
/      the appropriate routine.
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
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
 
/* function prototypes */
//#ifdef USE_HDF5_GROUPS
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1, int RestartDump = FALSE);
//#else
/* 
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1);
*/
//#endif

double ReturnWallTime(void);
void my_exit(int status);


int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData,
		   ExternalBoundary *Exterior, int &WroteData)
{
 
  /* Declarations. */
 
  char *Name;
  int i, Number;
  WroteData = FALSE;
 
  /* Check for output: restart-based. */

  char *param;
  FILE *pfptr;
  double CurrentCPUTime = ReturnWallTime() - MetaData.StartCPUTime;
  float FractionalCPUTime = 1.0 - MetaData.LastCycleCPUTime / MetaData.StopCPUTime;

  if ((CurrentCPUTime >= MetaData.dtRestartDump && 
       MetaData.dtRestartDump > 0 ) ||
      (MetaData.CycleNumber - MetaData.CycleLastRestartDump >= 
       MetaData.CycleSkipRestartDump &&
       MetaData.CycleSkipRestartDump > 0)) {

    MetaData.CycleLastRestartDump = MetaData.CycleNumber;

    if (debug) printf("Writing restart dump.\n");
    Group_WriteAllData(MetaData.RestartDumpName, MetaData.RestartDumpNumber++,
		       TopGrid, MetaData, Exterior);

    /* On the root processor, write the restart parameter filename to
       a file that will be read by a (batch) script to restart enzo.
       We cannot call another MPI application from here. */

    if (MyProcessorNumber == ROOT_PROCESSOR) {
      param = new char[512];
      if (MetaData.RestartDumpDir != NULL)
	sprintf(param, "%s%"CYCLE_TAG_FORMAT""ISYM"/%s%"CYCLE_TAG_FORMAT""ISYM,
		MetaData.RestartDumpDir, MetaData.RestartDumpNumber-1,
		MetaData.RestartDumpName, MetaData.RestartDumpNumber-1);
      else
	sprintf(param, "%s%"CYCLE_TAG_FORMAT""ISYM,
		MetaData.RestartDumpName, MetaData.RestartDumpNumber-1);

      if ((pfptr = fopen("RestartParamFile", "w")) == NULL)
	ENZO_FAIL("Error opening RestartParamFile");
      fprintf(pfptr, "%s", param);
      fclose(pfptr);
      
      delete [] param;
    } // ENDIF ROOT_PROCESSOR

    WroteData = TRUE;
    return SUCCESS;
  }
    
  /* Check for output: time-based. */
 
  if (MetaData.Time >= MetaData.TimeLastDataDump + MetaData.dtDataDump
      && MetaData.dtDataDump > 0.0) {
    MetaData.TimeLastDataDump += MetaData.dtDataDump;

    //#ifdef USE_HDF5_GROUPS
    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		       TopGrid, MetaData, Exterior);
// #else
//     if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
// 		     TopGrid, MetaData, Exterior) == FAIL) {
// 	fprintf(stderr, "Error in WriteAllData.\n");
// 	ENZO_FAIL("");
//     }
// #endif

    WroteData = TRUE;
  }
 
  /* Check for output: cycle-based. */
 
  if (MetaData.CycleNumber >= MetaData.CycleLastDataDump +
                              MetaData.CycleSkipDataDump   &&
      MetaData.CycleSkipDataDump > 0) {
    MetaData.CycleLastDataDump += MetaData.CycleSkipDataDump;

    //#ifdef USE_HDF5_GROUPS
    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		       TopGrid, MetaData, Exterior);
// #else
//     if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
// 		     TopGrid, MetaData, Exterior) == FAIL) {
// 	fprintf(stderr, "Error in WriteAllData.\n");
// 	ENZO_FAIL("");
//     }
// #endif

    WroteData = TRUE;
  }
 
  /* Check for output: CPU time-based.  

     If there is less time until 0.9*StopCPUTime than the last cycle's
     CPU time, output the data to ensure we get the last data dump for
     restarting! */

  //float FractionalCPUTime = 0.9;
  if (debug)
    printf("CPUTime-output: Frac = %"FSYM", Current = %lg (%lg), Stop = %"FSYM", "
	   "Last = %lg\n",
	   FractionalCPUTime, ReturnWallTime(), CurrentCPUTime, 
	   MetaData.StopCPUTime, MetaData.LastCycleCPUTime);
  if (CurrentCPUTime + MetaData.LastCycleCPUTime > 
      FractionalCPUTime*MetaData.StopCPUTime && MetaData.StartCPUTime > 0 &&
      WroteData == FALSE) {
    MetaData.CycleLastDataDump = MetaData.CycleNumber;
    if (debug) printf("CPUtime-based output!\n");
//#ifdef USE_HDF5_GROUPS
    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		       TopGrid, MetaData, Exterior);
// #else
//     if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
// 		     TopGrid, MetaData, Exterior) == FAIL) {
// 	fprintf(stderr, "Error in WriteAllData.\n");
// 	ENZO_FAIL("");
//     }
// #endif
    WroteData = TRUE;
  } // ENDIF

  /* Check for output: redshift-based. */
 
  if (ComovingCoordinates)
    for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
      if (CosmologyOutputRedshift[i] != -1)
	if (MetaData.Time >= CosmologyOutputRedshiftTime[i]) {
	  CosmologyOutputRedshift[i] = -1; // done, turn it off
	  if (CosmologyOutputRedshiftName[i] == NULL) {
	    Name   = MetaData.RedshiftDumpName;
	    Number = i;   // append number to end of name
	  }
	  else {
	    Name   = CosmologyOutputRedshiftName[i];
	    Number = -1;  // Don't append number (####) to end of name
	  }

	  //#ifdef USE_HDF5_GROUPS
	  Group_WriteAllData(Name, Number, TopGrid, MetaData, Exterior);
// #else
// 	  if (WriteAllData(Name, Number, TopGrid, MetaData, Exterior) == FAIL) {
// 	    fprintf(stderr, "Error in WriteAllData.\n");
// 	    ENZO_FAIL("");
// 	  }
// #endif

	  WroteData = TRUE;
	}

  /* Check for output: when the MBH jets haven't been ejected for too long 
                       this is currently a test - Ji-hoon Kim, Mar.2010 */  //#####
 
  if ((MBHFeedback == 2 || MBHFeedback == 3) && 
      OutputWhenJetsHaveNotEjected == TRUE) {

    fprintf(stdout, "CheckForOutput: MBH_JETS - file output complete; restart with the dump!\n");
    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		       TopGrid, MetaData, Exterior);

    OutputWhenJetsHaveNotEjected = FALSE;
    WroteData = TRUE;
    my_exit(EXIT_SUCCESS);

  }
  
  return SUCCESS;
}
