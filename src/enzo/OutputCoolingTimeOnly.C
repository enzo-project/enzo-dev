/***********************************************************************
/
/  COMPUTE THE COOLING TIME AND OUTPUT
/
/  written by: John Wise
/  date:       April, 2011
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

int Group_WriteAllData(char *basename, int filenumber,
		       HierarchyEntry *TopGrid, TopGridData &MetaData,
		       ExternalBoundary *Exterior, 
#ifdef TRANSFER
		       ImplicitProblemABC *ImplicitSolver,
#endif
		       FLOAT WriteTime = -1, 
		       int CheckpointDump = FALSE);

int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
#ifdef TRANSFER
int RadiativeTransferReadParameters(FILE *fptr);
#endif

int OutputCoolingTimeOnly(char *ParameterFile,
			  LevelHierarchyEntry *LevelArray[], 
			  HierarchyEntry *TopGrid,
			  TopGridData &MetaData,
			  ExternalBoundary &Exterior
#ifdef TRANSFER
			  , ImplicitProblemABC *ImplicitSolver
#endif
			  )
{

  int i, level;
  const float When = 0.5;
  LevelHierarchyEntry *Temp;
  HierarchyEntry **Grids;

  /* Exit if OutputCoolingTime is already on. */

  if (OutputCoolingTime == TRUE) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      printf("OutputCoolingTime is already on.  "
	     "Change to 0 in %s to recompute.\n", ParameterFile);
    return SUCCESS;
  }

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

  /* Initialize radiative transfer parameters, if needed */

#ifdef TRANSFER
  FILE *fptr;

  if ((fptr = fopen(ParameterFile, "r")) == NULL) {
    ENZO_VFAIL("Error opening ParameterFile %s\n", ParameterFile)
  }

  RadiativeTransferReadParameters(fptr);

  fclose(fptr);
#endif /* TRANSFER */

  // Negative number to indicate that this won't propagate to the
  // parameter, and only compute the cooling time.
  OutputCoolingTime = -1;

  Group_WriteAllData(DumpName, DumpNumber, TopGrid, MetaData, &Exterior
#ifdef TRANSFER
		       , ImplicitSolver
#endif
               );
  

  return SUCCESS;

}
