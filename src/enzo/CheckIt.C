/***********************************************************************
/
/  Check Config
/
/  written by: Robert Harkness
/  date:       August 19th 2006
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "units.h"
#include "flowdefs.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#undef DEFINE_STORAGE
 
// Function prototypes
 
int InitializeNew(  char *filename, HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary &Exterior, float *Initialdt);
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
int EvolveHierarchy(                HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, LevelHierarchyEntry *Array[],
		    float Initialdt);
void ExtractSection(HierarchyEntry &TopGrid, TopGridData &tgd,
		    LevelHierarchyEntry *Array[], ExternalBoundary *Exterior,
		    int ExtractStart[], int ExtractEnd[],
		    FLOAT ProjectStartCoordinates[],
		    FLOAT ProjectEndCoordinates[], int ExtractLevel);
int OutputLevelInformation(FILE *fptr, TopGridData &tgd,
			   LevelHierarchyEntry *Array[]);
int ProjectToPlane(TopGridData &MetaData, LevelHierarchyEntry *LevelArray[],
		   int ProjectStart[], int ProjectEnd[],
		   FLOAT ProjectStartCoordinates[],
		   FLOAT ProjectEndCoordinates[], int ProjectLevel,
		   int ProjectionDimension, char *ProjectionFileName,
		   int ProjectionSmooth, ExternalBoundary *Exterior);
int OutputAsParticleData(TopGridData &MetaData,
			 LevelHierarchyEntry *LevelArray[],
			 int RegionStart[], int RegionEnd[],
			 FLOAT RegionStartCoordinates[],
			 FLOAT RegionEndCoordinates[], int RegionLevel,
			 char *OutputFileName);
int InterpretCommandLine(int argc, char *argv[], char *myname,
			 int &restart, int &debug, int &extract,
			 int &InformationOutput,
			 int &OutputAsParticleDataFlag,
			 int &project, int &ProjectionDimension,
			 int &ProjectionSmooth,
			 char *ParameterFile[],
			 int RegionStart[], int RegionEnd[],
			 FLOAT RegionStartCoordinates[],
			 FLOAT RegionEndCoordinates[],
			 int &Level, int MyProcessorNumber);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
int CommunicationInitialize(int *argc, char **argv[]);
int CommunicationFinalize();
int CommunicationPartitionGrid(HierarchyEntry *Grid, int gridnum);
int ENZO_OptionsinEffect(void);

#ifdef TASKMAP
int GetNodeFreeMemory(void);
#endif

void my_exit(int status);
 
#ifdef MEM_TRACE
Eint64 mused(void);
#endif 
 
 
Eint32 main(Eint32 argc, char *argv[])
{

  // Initialize Communications

  int int_argc;
  int_argc = argc;
 
  CommunicationInitialize(&int_argc, &argv);
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("ENZO V64 19th August 2006\n");

#ifdef OOC_BOUNDARY
  ExternalBoundaryIO = TRUE;
  ExternalBoundaryTypeIO = FALSE;
  ExternalBoundaryValueIO = FALSE;
#else
  ExternalBoundaryIO = FALSE;
  ExternalBoundaryTypeIO = FALSE;
  ExternalBoundaryValueIO = FALSE;
#endif

  ENZO_OptionsinEffect();

  my_exit(EXIT_SUCCESS);

}

 
void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}
