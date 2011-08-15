/***********************************************************************
/
/  COMPUTE THE SMOOTHED DM FIELD AND OUTPUT
/
/  written by: John Wise
/  date:       May, 2010
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/

#include "preincludes.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
 
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
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif
#include "CosmologyParameters.h"

int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
		       ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1, int RestartDump = FALSE);

int OutputSmoothedDarkMatterOnly(char *ParameterFile,
				 LevelHierarchyEntry *LevelArray[], 
				 HierarchyEntry *TopGrid,
				 TopGridData &MetaData,
				 ExternalBoundary &Exterior
#ifdef TRANSFER
		       , ImplicitProblemABC *ImplicitSolver
#endif
                )
{

  int level;

  /* Determine the parameter name prefix */

  char *cptr;
  int DumpNumber;
  char DumpName[MAX_LINE_LENGTH];
  char NumberString[MAX_LINE_LENGTH];
  if ( (cptr = strstr(ParameterFile, MetaData.DataDumpName)) ) {
    strcpy(DumpName, MetaData.DataDumpName);
  }
  else if ( (cptr = strstr(ParameterFile, MetaData.RestartDumpName)) ) {
    strcpy(DumpName, MetaData.RestartDumpName);
  }
  else if ( (cptr = strstr(ParameterFile, MetaData.RedshiftDumpName)) ) {
    strcpy(DumpName, MetaData.RedshiftDumpName);
  }
  else
    ENZO_FAIL("Cannot determine output type.");

  /* Extract output number */
  
  strncpy(NumberString, ParameterFile+2, 4);
  NumberString[4] = '\0';
  DumpNumber = atoi(NumberString);

  /* If we're not using parallel root grid I/O and we're parallel, we
     need to rebuild the hierarchy with the multiple root grids. */

  if (!ParallelRootGridIO && NumberOfProcessors > 1)
    RebuildHierarchy(&MetaData, LevelArray, 0);

  // Negative number means that it'll be reset to zero after
  // calculating the DM field so it doesn't propagate to later runs
  if (OutputSmoothedDarkMatter == FALSE)
    OutputSmoothedDarkMatter = -2;

  Group_WriteAllData(DumpName, DumpNumber, TopGrid, MetaData, &Exterior
#ifdef TRANSFER
		       , ImplicitSolver
#endif
               );
  
  return SUCCESS;

}
