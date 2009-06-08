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
 
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1);

int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		       TopGridData &MetaData, ExternalBoundary *Exterior,
		       FLOAT WriteTime = -1);
double ReturnWallTime(void);


int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData,
		   ExternalBoundary *Exterior, int &WroteData)
{
 
  /* Declarations. */
 
  char *Name;
  int i, Number;
  WroteData = FALSE;
 
  /* Check for output: time-based. */
 
  if (MetaData.Time >= MetaData.TimeLastDataDump + MetaData.dtDataDump
      && MetaData.dtDataDump > 0.0) {
    MetaData.TimeLastDataDump += MetaData.dtDataDump;

#ifdef USE_HDF5_GROUPS
    if (Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	return FAIL;
    }
#else
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
    }
#endif

    WroteData = TRUE;
  }
 
  /* Check for output: cycle-based. */
 
  if (MetaData.CycleNumber >= MetaData.CycleLastDataDump +
                              MetaData.CycleSkipDataDump   &&
      MetaData.CycleSkipDataDump > 0) {
    MetaData.CycleLastDataDump += MetaData.CycleSkipDataDump;

#ifdef USE_HDF5_GROUPS
    if (Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	return FAIL;
    }
#else
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
    }
#endif

    WroteData = TRUE;
  }
 
  /* Check for output: CPU time-based.  

     If there is less time until 0.9*StopCPUTime than the last cycle's
     CPU time, output the data to ensure we get the last data dump for
     restarting! */

  double CurrentCPUTime = ReturnWallTime() - MetaData.StartCPUTime;
  float FractionalCPUTime = 1.0 - MetaData.LastCycleCPUTime / MetaData.StopCPUTime;
  //float FractionalCPUTime = 0.9;
  if (debug)
    printf("CPUTime-output: Frac = %"FSYM", Current = %lg (%lg), Stop = %"FSYM", "
	   "Last = %lg\n",
	   FractionalCPUTime, ReturnWallTime(), CurrentCPUTime, 
	   MetaData.StopCPUTime, MetaData.LastCycleCPUTime);
  if (CurrentCPUTime + MetaData.LastCycleCPUTime > 
      FractionalCPUTime*MetaData.StopCPUTime && MetaData.StartCPUTime > 0) {
    MetaData.CycleLastDataDump = MetaData.CycleNumber;
    if (debug) printf("CPUtime-based output!\n");
#ifdef USE_HDF5_GROUPS
    if (Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	return FAIL;
    }
#else
    if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		     TopGrid, MetaData, Exterior) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
    }
#endif
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

#ifdef USE_HDF5_GROUPS
	  if (Group_WriteAllData(Name, Number, TopGrid, MetaData, Exterior) == FAIL) {
	    fprintf(stderr, "Error in Group_WriteAllData.\n");
	    return FAIL;
	  }
#else
	  if (WriteAllData(Name, Number, TopGrid, MetaData, Exterior) == FAIL) {
	    fprintf(stderr, "Error in WriteAllData.\n");
	    return FAIL;
	  }
#endif

	  WroteData = TRUE;
	}
 
  return SUCCESS;
}
