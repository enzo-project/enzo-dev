/***********************************************************************
/
/  AMR MAIN CODE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       August 12th 2006
/              May 13th 2008
/
/  PURPOSE:
/    This is main() for the amr code.  It interprets the arguments and
/    then calls the appropriate routines depending on whether we are
/    doing a new run, a restart, an extraction, or a projection.
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
 
#include "ErrorExceptions.h"
#include "svn_version.def"
#include "performance.h"
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
#include "communication.h"
#include "CommunicationUtilities.h"
#ifdef TRANSFER
#include "PhotonCommunication.h"
#endif
#undef DEFINE_STORAGE
#ifdef USE_PYTHON
int InitializePythonInterface(int argc, char **argv);
#endif
 
// Function prototypes
 
int InitializeNew(  char *filename, HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary &Exterior, float *Initialdt);
int InitializeMovieFile(TopGridData &MetaData, HierarchyEntry &TopGrid);

int InitializeLocal(int restart, HierarchyEntry &TopGrid, 
		    TopGridData &MetaData);

int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
int Group_ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);

int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior, 
		    LevelHierarchyEntry *Array[], float Initialdt);

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
int CommunicationInitialize(Eint32 *argc, char **argv[]);
int CommunicationFinalize();

int CommunicationPartitionGrid(HierarchyEntry *Grid, int gridnum);
void CommunicationAbort(int);
int ENZO_OptionsinEffect(void);

#ifdef TASKMAP
int GetNodeFreeMemory(void);
#endif

#ifdef TRANSFER
int RadiativeTransferInitialize(char *ParameterFile, TopGridData &MetaData,
				ExternalBoundary &Exterior, 
				LevelHierarchyEntry *LevelArray[]);
#endif

#ifdef USE_JBPERF
void jbPerfInitialize (int max_level);
#endif

void my_exit(int status);

 
#ifdef MEM_TRACE
Eint64 mused(void);
#endif 


 
//  ENZO Main Program

 
Eint32 main(Eint32 argc, char *argv[])
{

#ifdef MEM_TRACE
    Eint64 MemInUse;
#endif

  int i;

  // Initialize Communications

  CommunicationInitialize(&argc, &argv);

#ifdef DEBUG_MPI
  int impi = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready for debugger attach\n", getpid(), hostname);
  fflush(stdout);
  while (impi == 0)
    sleep(5);
#endif
  

  int int_argc;
  int_argc = argc;
 
  if (MyProcessorNumber == ROOT_PROCESSOR &&
      ENZO_SVN_REVISION != 0) {
    printf("=========================\n");
    printf("Enzo SVN Branch   %s\n",ENZO_SVN_BRANCH);
    printf("Enzo SVN Revision %d\n",ENZO_SVN_REVISION);
    printf("=========================\n");
    fflush(stdout);
  }
  // Performance Monitoring

#ifdef USE_MPI
  double t_init0, t_init1;

  t_init0 = MPI_Wtime();
#endif /* USE_MPI */

#ifdef USE_JBPERF

    // Initialize jbPerf performance collecting

  jbPerfInitialize(MaximumRefinementLevel);

#endif


  // Task Mapping

  for (i = 0; i < MAX_NUMBER_OF_TASKS; i++ ) {
    TaskMemory[i] = -1;
  }

#ifdef TASKMAP
  GetNodeFreeMemory();
#endif


  // Begin 

#ifdef MPI_INSTRUMENTATION
  Start_Wall_Time = MPI_Wtime();
#endif

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


  // Main declarations
 
  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

 
  // Initialize
 
  int restart                  = FALSE,
      OutputAsParticleDataFlag = FALSE,
      InformationOutput        = FALSE,
      project                  = FALSE,
      ProjectionDimension      = INT_UNDEFINED,
      ProjectionSmooth         = FALSE;
  debug                        = FALSE;
  extract                      = FALSE;
  flow_trace_on                = FALSE;
  char *ParameterFile          = NULL;
  char *myname                 = argv[0];

  int RegionStart[MAX_DIMENSION],
      RegionEnd[MAX_DIMENSION],
      RegionLevel;
  FLOAT RegionStartCoordinates[MAX_DIMENSION],
        RegionEndCoordinates[MAX_DIMENSION];

  float Initialdt              = 0;
 
#ifdef FLOW_TRACE
  char pid[MAX_TASK_TAG_SIZE];
  char flow_trace_log[MAX_NAME_LENGTH];
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
  strcpy(flow_trace_log, "FlowTrace.");
  strcat(flow_trace_log, pid);
  flow_trace_fptr = fopen( flow_trace_log, "w" );
  flow_trace_level = 0;
  flow_trace_on = TRUE;
#endif
 
  if (flow_trace_on) flow_trace1("ENZO");
 
#ifdef MPI_INSTRUMENTATION
  char perfname[MAX_NAME_LENGTH];
 
  flagging_count = 0;
  moving_count = 0;
  in_count = 0;
  out_count = 0;
  flagging_pct = 0.0;
  moving_pct = 0.0;
 
  GlobalCommunication = 0.0;
  RecvComm = 0.0;
  WaitComm = 0.0;
 
  for (i=0; i < MAX_COUNTERS; i++) {
    timer[i] = 0.0;
    counter[i] = 0;
  }
 
  sprintf(perfname, "perfdata_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  perfname[strlen(perfname)] = '\0';
 
  sprintf(tracename, "trace_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  tracename[strlen(tracename)] = '\0';

#ifdef MPI_TRACE
  tracePtr = fopen(tracename, "w");
  traceMPI = TRUE;
#else
  traceMPI = FALSE;
  tracePtr = NULL;
#endif /* MPI_TRACE */

#endif /* MPI_INSTRUMENTATION */

#ifdef MEM_TRACE
  sprintf(memtracename, "mem_%"ISYM".%"ISYM, NumberOfProcessors, MyProcessorNumber);
  memtracename[strlen(memtracename)] = '\0';
  memtracePtr = fopen(memtracename, "w");
  traceMEM = TRUE;
#endif /* MEM_TRACE */


  // START
 
  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    LevelArray[level] = NULL;
  }
 

  // Interpret command-line arguments
 
  if (InterpretCommandLine(int_argc, argv, myname, restart, debug, extract,
			   InformationOutput, OutputAsParticleDataFlag,
			   project, ProjectionDimension, ProjectionSmooth,
			   &ParameterFile,
			   RegionStart, RegionEnd,
			   RegionStartCoordinates, RegionEndCoordinates,
			   RegionLevel, MyProcessorNumber) == FAIL) {
    my_exit(EXIT_FAILURE);
  }

 
  // If we need to read the parameter file as a restart file, do it now
 
  if (restart || OutputAsParticleDataFlag || extract || InformationOutput || project) {
 
    SetDefaultGlobalValues(MetaData);
 
    if (debug) printf("Reading parameter file %s\n", ParameterFile);


  // First expect to read in packed-HDF5

#ifdef USE_HDF5_GROUPS
    if (Group_ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR) {
	fprintf(stderr, "Error in Group_ReadAllData %s\n", ParameterFile);
	fprintf(stderr, "Probably not in a packed-HDF5 format. Trying other read routines.\n");
      }
#endif
      // If not packed-HDF5, then try usual HDF5 or HDF4
      if (ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior) == FAIL) {
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  fprintf(stderr, "Error in ReadAllData %s.\n", ParameterFile);
	my_exit(EXIT_FAILURE);
      }
#ifdef USE_HDF5_GROUPS
    }
#endif

    /*
#if defined(USE_HDF5_GROUPS) && !defined(USE_HDF4)
    if (Group_ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Error in Group_ReadAllData %s.\n", ParameterFile);
      my_exit(EXIT_FAILURE);
    }
#else 
    if (ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Error in ReadAllData %s.\n", ParameterFile);
      my_exit(EXIT_FAILURE);
    }
#endif
    */

 
    if (!ParallelRootGridIO && restart && TopGrid.NextGridThisLevel == NULL) {
      CommunicationPartitionGrid(&TopGrid, 0);  // partition top grid if necessary
    }
 
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "Successfully read ParameterFile %s.\n", ParameterFile);
    }
 
    AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels

    // Initialize streaming files if necessary
    InitializeMovieFile(MetaData, TopGrid);
 
  }


  // Information dump
 
  if (InformationOutput) {
    OutputLevelInformation(stdout, MetaData, LevelArray);
    my_exit(EXIT_SUCCESS);
  }
 
 
  // Extract a grid subsection (doesn't return)
 
  if (extract) {
    ExtractSection(TopGrid, MetaData, LevelArray, &Exterior,
		   RegionStart, RegionEnd,
		   RegionStartCoordinates, RegionEndCoordinates, RegionLevel);
  }
 
 
  // Output fields as particle data
 
  if (OutputAsParticleDataFlag) {
    if (OutputAsParticleData(MetaData, LevelArray, RegionStart, RegionEnd,
			     RegionStartCoordinates, RegionEndCoordinates,
			     RegionLevel, "amr.particles") == FAIL)
      my_exit(EXIT_FAILURE);
    else
      my_exit(EXIT_SUCCESS);
  } 

 
  // Project 3D field to 2D plane
 
  if (project) {
    if (ProjectToPlane(MetaData, LevelArray, RegionStart, RegionEnd,
		     RegionStartCoordinates, RegionEndCoordinates,
		     RegionLevel, ProjectionDimension, "amr.project",
		     ProjectionSmooth, &Exterior) == FAIL)
      my_exit(EXIT_FAILURE);
    else
      my_exit(EXIT_SUCCESS);
  } 

 
  // Normal start: Open and read parameter file

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Call initialize %16"ISYM" \n", MemInUse);
#endif

 
  if (!restart) {

    if (InitializeNew(ParameterFile, TopGrid, MetaData, Exterior, &Initialdt) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Error in Parameter File %s.\n", ParameterFile);
      my_exit(EXIT_FAILURE);
    }
    else {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "Successfully read in parameter file %s.\n", ParameterFile);
      AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels
    }
 
#ifdef USE_MPI
    CommunicationBarrier();
    t_init1 = MPI_Wtime();
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "INITIALIZATION TIME = %16.8e\n", (t_init1-t_init0));
    CommunicationBarrier();
#endif /* USE_MPI */

  }


/*

  // Perform local initialization (even for restart, may just return, depends on problem)

  if (InitializeLocal(restart, TopGrid, MetaData) == FAIL) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "Error in Local Initialization.\n");
    my_exit(EXIT_FAILURE); 
  }

*/ 

  /* Initialize the radiative transfer */

#ifdef TRANSFER
  if (RadiativeTransferInitialize(ParameterFile, MetaData, Exterior, 
				  LevelArray) == FAIL) {
    fprintf(stderr, "Error in RadiativeTransferInitialize.\n");
    my_exit(EXIT_FAILURE);
  }
#endif


#ifdef MEM_TRACE
  MemInUse = mused();
  fprintf(memtracePtr, "Call evolve hierarchy %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif

#ifdef USE_PYTHON
  // We initialize our Python interface now
  if(debug)fprintf(stdout, "Initializing Python interface\n");
  InitializePythonInterface(argc, argv);
#endif 


 
  // Call the main evolution routine
 
  try {
  if (EvolveHierarchy(TopGrid, MetaData, &Exterior, LevelArray, Initialdt) == FAIL) {
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "Error in EvolveHierarchy.\n");
    }
    // Close the streaming files
    if (MovieSkipTimestep != INT_UNDEFINED) {
      fprintf(stderr, "Closing movie file.\n");
      MetaData.AmiraGrid.AMRHDF5Close();
    }
    my_exit(EXIT_FAILURE);
  }
  else
  {
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "Successful run, exiting.\n");
    }
  }
  } catch(EnzoFatalException&) {
    fprintf(stderr, "Failure reported on processor %"ISYM"\n",
                MyProcessorNumber);
    CommunicationAbort(EXIT_FAILURE);
  }

 
  if (flow_trace_on) flow_trace2("ENZO");




  // Terminate
  
#ifdef MPI_INSTRUMENTATION

  End_Wall_Time = MPI_Wtime();
  WallTime = End_Wall_Time - Start_Wall_Time;
 
  filePtr = fopen(perfname, "w");
  fprintf(filePtr, "Elapsed wall time:                   %12.6e\n", WallTime);
  fprintf(filePtr, "Communication time:                  %12.6e\n", CommunicationTime);
  fprintf(filePtr, "Global communication time:           %12.6e\n", GlobalCommunication);
  fprintf(filePtr, "Receive communication time:          %12.6e\n", RecvComm);
  fprintf(filePtr, "Waiting communication time:          %12.6e\n", WaitComm);
 
  fprintf(filePtr, "\n\n");
  fprintf(filePtr, "Transferring region       (%8"ISYM" times) %12.6e\n", counter[5],  timer[5]);
  fprintf(filePtr, "Sending particles         (%8"ISYM" times) %12.6e\n", counter[7],  timer[7]);
  fprintf(filePtr, "Transferring particles    (%8"ISYM" times) %12.6e\n", counter[9],  timer[9]);
  fprintf(filePtr, "Transferring Fluxes       (%8"ISYM" times) %12.6e\n", counter[12], timer[12]);
  fprintf(filePtr, "ShareGrids                (%8"ISYM" times) %12.6e\n", counter[13], timer[13]);
  fprintf(filePtr, "Transpose                 (%8"ISYM" times) %12.6e\n", counter[14], timer[14]);
  fprintf(filePtr, "BroadcastValue            (%8"ISYM" times) %12.6e\n", counter[15], timer[15]);
  fprintf(filePtr, "MinValue                  (%8"ISYM" times) %12.6e\n", counter[16], timer[16]);
  fprintf(filePtr, "UpdateStarParticleCount   (%8"ISYM" times) %12.6e\n", counter[11], timer[11]);
 
  fprintf(filePtr, "\n\n");
  fprintf(filePtr, "RebuildHierarchy          (%8"ISYM" times) %12.6e\n", counter[1],  timer[1]);
  fprintf(filePtr, "RebuildHierarchy interval (%8"ISYM" times) %12.6e\n", counter[0],  timer[0]);
  fprintf(filePtr, "Load balancing            (%8"ISYM" times) %12.6e\n", counter[2],  timer[2]);
  fprintf(filePtr, "Region transfer size      (%8"ISYM" times) %12.6e\n", counter[5],  timer[6]);
  fprintf(filePtr, "Particles sent            (%8"ISYM" times) %12.6e\n", counter[7],  timer[8]);
  fprintf(filePtr, "Particle transfer size    (%8"ISYM" times) %12.6e\n", counter[9],  timer[10]);
 
  fprintf(filePtr, "\n\n");
  fprintf(filePtr, "Number of load balancing calls %"ISYM"/%"ISYM" (LOAD_BALANCE_RATIO=%"FSYM")\n",counter[3], counter[2], timer[3]);
  fprintf(filePtr, "Number of flagging cells  (%8"ISYM" times) %12.6e\n", counter[4], timer[4]);
 
  fprintf(filePtr, "\n\n");

  if ( flagging_count != 0 )
    fprintf(filePtr, "Average percentage of flagging cells %12.6e(= %12.6e/%"ISYM")\n",
            flagging_pct/flagging_count, flagging_pct, flagging_count);
  else
    fprintf(filePtr, "Average percentage of flagging cells 0\n");
 
  if ( moving_count != 0 )
    fprintf(filePtr, "Average percentage of moving cells   %12.6e(= %12.6e/%"ISYM")\n",
            moving_pct/moving_count, moving_pct,moving_count);
  else
    fprintf(filePtr, "Average percentage of moving cells 0\n");
 
  fclose(filePtr);
 
#ifdef MPI_TRACE
  fclose(tracePtr);
#endif /* MPI_TRACE */

#endif /* MPI_INSTRUMENTATION */

#ifdef MEM_TRACE
  fclose(memtracePtr);
#endif /* MEM_TRACE */
 
#ifdef FLOW_TRACE
    fclose(flow_trace_fptr);
#endif
 
  my_exit(EXIT_SUCCESS);
 
}
 
void my_exit(int status)
{
  // Exit gracefully if successful; abort on error

  if (status == EXIT_SUCCESS) {

    if (MyProcessorNumber==0) {
      fprintf (stdout,"%s:%d Exiting.\n", __FILE__,__LINE__);
    }

    CommunicationFinalize();

    exit(status);

  } else if (status == EXIT_FAILURE) {

    fprintf (stderr,"%s:%d %"ISYM" ABORT ON EXIT_FAILURE!\n",
	     __FILE__,__LINE__,MyProcessorNumber);

    ENZO_FAIL("my_exit has been called.");
    CommunicationAbort(status);

  } else {

    fprintf (stderr,"%s:%d %"ISYM" ABORT ON UNKNOWN EXIT VALUE %"ISYM"!\n",
	     __FILE__,__LINE__,MyProcessorNumber,status);

    ENZO_FAIL("my_exit has been called without a known exit value.");
    CommunicationAbort(status);

  }
}

