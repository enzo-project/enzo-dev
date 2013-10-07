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
#include "preincludes.h"
 
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
#include "CommunicationUtilities.h"
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif
 
/* function prototypes */
//#ifdef USE_HDF5_GROUPS
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
	         ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1, int RestartDump = FALSE);
//#else
/* 
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
	         ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1);
*/
//#endif

double ReturnWallTime(void);
void my_exit(int status);


int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData,
		   ExternalBoundary *Exterior, 
#ifdef TRANSFER
	           ImplicitProblemABC *ImplicitSolver,
#endif
		   int Restart)
{
 
  /* Declarations. */
 
  char *Name;
  int i, Number;
  MetaData.WroteData = FALSE;
  double SavedCPUTime;

  /* Check for output: restart-based. */

  char *param;
  FILE *pfptr;

  if (Restart == TRUE && MetaData.WroteData == FALSE) {

    MetaData.CycleLastRestartDump = MetaData.CycleNumber;

    if (debug) printf("Writing restart dump.\n");
    Group_WriteAllData(MetaData.RestartDumpName, MetaData.RestartDumpNumber++,
		       TopGrid, MetaData, Exterior
#ifdef TRANSFER
		       , ImplicitSolver
#endif
		       );


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

    MetaData.WroteData = TRUE;
    return SUCCESS;
  }
    
  /* Check for output: CPU time-based.  

     If there is less time until (StopCPUTime - LastCycleCPUTime) than
     the last cycle's CPU time, output the data to ensure we get the
     last data dump for restarting! */

  float FractionalCPUTime = 1.0 - MetaData.LastCycleCPUTime / MetaData.StopCPUTime;

  if (debug)
    printf("CPUTime-output: Frac = %"FSYM", Current = %lg (%lg), Stop = %"FSYM", "
	   "Last = %lg\n",
	   FractionalCPUTime, ReturnWallTime()-MetaData.StartCPUTime, 
	   MetaData.CPUTime, 
	   MetaData.StopCPUTime, MetaData.LastCycleCPUTime);
  if (MetaData.CPUTime + MetaData.LastCycleCPUTime > 
      FractionalCPUTime*MetaData.StopCPUTime && MetaData.StartCPUTime > 0 &&
      MetaData.WroteData == FALSE) {
    MetaData.CycleLastDataDump = MetaData.CycleNumber;
    SavedCPUTime = MetaData.CPUTime;
    MetaData.CPUTime = 0.0;
    if (debug) printf("CPUtime-based output!\n");
    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		       TopGrid, MetaData, Exterior
#ifdef TRANSFER
		       , ImplicitSolver
#endif
		       );
    MetaData.CPUTime = SavedCPUTime;
    MetaData.WroteData = TRUE;
  } // ENDIF

  /* Check for output: time-based. */
 
  if (MetaData.Time >= MetaData.TimeLastDataDump + MetaData.dtDataDump
      && MetaData.dtDataDump > 0.0) {
    SavedCPUTime = MetaData.CPUTime;
    MetaData.CPUTime = 0.0;
    MetaData.TimeLastDataDump += MetaData.dtDataDump;

    //#ifdef USE_HDF5_GROUPS
    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		       TopGrid, MetaData, Exterior
#ifdef TRANSFER
		       , ImplicitSolver
#endif
		       );
// #else
//     if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
// 		        TopGrid, MetaData, Exterior
//#ifdef TRANSFER
//			, ImplicitSolver
//#endif
//		        ) == FAIL {
// 	ENZO_FAIL("Error in WriteAllData.\n");
//     }
// #endif

    MetaData.CPUTime = SavedCPUTime;
    MetaData.WroteData = TRUE;
  }
 
  /* Check for output: cycle-based. */
 
  if (MetaData.CycleNumber >= MetaData.CycleLastDataDump +
                              MetaData.CycleSkipDataDump   &&
      MetaData.CycleSkipDataDump > 0) {
    SavedCPUTime = MetaData.CPUTime;
    MetaData.CPUTime = 0.0;
    MetaData.CycleLastDataDump += MetaData.CycleSkipDataDump;

    //#ifdef USE_HDF5_GROUPS
    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		       TopGrid, MetaData, Exterior
#ifdef TRANSFER
		       , ImplicitSolver
#endif
		       );
// #else
//     if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
// 		        TopGrid, MetaData, Exterior
//#ifdef TRANSFER
//			, ImplicitSolver
//#endif
//                      ) == FAIL) {
// 	ENZO_FAIL("Error in WriteAllData.\n");
//     }
// #endif

    MetaData.CPUTime = SavedCPUTime;
    MetaData.WroteData = TRUE;
  }
 
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

	  SavedCPUTime = MetaData.CPUTime;
	  MetaData.CPUTime = 0.0;
	  //#ifdef USE_HDF5_GROUPS
	  Group_WriteAllData(Name, Number, TopGrid, MetaData, Exterior
#ifdef TRANSFER
			     , ImplicitSolver
#endif
			     );
// #else
// 	  if (WriteAllData(Name, Number, TopGrid, MetaData, Exterior
//#ifdef TRANSFER
//			   , ImplicitSolver
//#endif
//                         ) == FAIL) {
// 	    ENZO_FAIL("Error in WriteAllData.\n");
// 	  }
// #endif

	  MetaData.CPUTime = SavedCPUTime;
	  MetaData.WroteData = TRUE;
	}

#ifdef UNUSED
  /* Check for output: when the MBH jets haven't been ejected for too long 
                       this is currently a test - Ji-hoon Kim, Mar.2010 */  
 
  if ((MBHFeedback >= 2 && MBHFeedback <= 5) && 
      OutputWhenJetsHaveNotEjected == TRUE) {

    fprintf(stdout, "CheckForOutput: MBH_JETS - file output complete; restart with the dump!\n");
    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber++,
		       TopGrid, MetaData, Exterior);

    OutputWhenJetsHaveNotEjected = FALSE;
    MetaData.WroteData = TRUE;
    my_exit(EXIT_SUCCESS);

  }
#endif   

  return SUCCESS;
}
