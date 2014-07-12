/***********************************************************************
/
/  COMPUTE THE POTENTIAL FIELD AND OUTPUT
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

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

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
#ifdef FAST_SIB
int CreateSiblingList(HierarchyEntry **Grids, int NumberOfGrids, 
		      SiblingGridList *SiblingList, int StaticLevelZero,
		      TopGridData *MetaData, int level);
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			int level, TopGridData *MetaData, FLOAT When);
#else  // !FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[], 
                        int level, TopGridData *MetaData, FLOAT When);
#endif  // end FAST_SIB

int OutputPotentialFieldOnly(char *ParameterFile,
			     LevelHierarchyEntry *LevelArray[], 
			     HierarchyEntry *TopGrid,
			     TopGridData &MetaData,
			     ExternalBoundary &Exterior,
#ifdef TRANSFER
		         ImplicitProblemABC *ImplicitSolver,
#endif
			     int OutputDM)
{

  int level;
  const float When = 0.5;
  LevelHierarchyEntry *Temp;
  HierarchyEntry **Grids;

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

  strcat(MetaData.GlobalDir, "/GravPotential");

  WritePotential = TRUE;
  CopyGravPotential = TRUE;
  OutputParticleTypeGrouping = FALSE;

  /* Add the potential field pointer to BaryonField and field names */

  int TypesToAdd[] = {GravPotential};
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->AddFields(TypesToAdd, 1);

  Exterior.AddField(GravPotential);

  int n = LevelArray[0]->GridData->ReturnNumberOfBaryonFields();
  DataLabel[n-1] = "PotentialField";

  /* If we're not using parallel root grid I/O and we're parallel, we
     need to rebuild the hierarchy with the multiple root grids. */

  if (!ParallelRootGridIO && NumberOfProcessors > 1)
    RebuildHierarchy(&MetaData, LevelArray, 0);

  /* Compute the potential field on all levels */

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    if (LevelArray[level] != NULL) {

#ifdef FAST_SIB
      int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
      SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];
      CreateSiblingList(Grids, NumberOfGrids, SiblingList, FALSE, &MetaData,
			level);
      PrepareDensityField(LevelArray,  level, &MetaData, When);
      delete [] SiblingList;
      delete [] Grids;
#else
      PrepareDensityField(LevelArray, level, &MetaData, When);
#endif

      if (level > 0)
	for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel) {
	  Temp->GridData->SolveForPotential(level);
	} // ENDFOR grids

    } // ENDIF grids on level
  } // ENDFOR levels

  WritePotential = FALSE;
  CopyGravPotential = FALSE;

  if (OutputDM == TRUE)
    OutputSmoothedDarkMatter = -2;

  Group_WriteAllData(DumpName, DumpNumber, TopGrid, MetaData, &Exterior
#ifdef TRANSFER
		       , ImplicitSolver
#endif
               );
  

  return SUCCESS;

}
