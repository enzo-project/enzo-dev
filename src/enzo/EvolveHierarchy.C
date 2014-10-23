/***********************************************************************
/
/  EVOLVE HIERARCHY FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  February, 1995 by GB
/              Changed to reflect changes in EvolveGrid & EvolveTopGrid.
/  modified2:  July, 1995 by GB
/              Changed to reflect new routine EvolveLevel.
/  modified3:  February, 2006 by Daniel Reynolds
/              Updated call interface to ComputePotentialFieldLevelZero
/  modified4:  February, 2007 by Robert Harkness
/              Group and in-core i/o
/  modified5:  December, 2007 by Robert Harkness
/              Remove calls to ComputePotential to allow use of the
/              Fast Sibling Locator to remedy Ncpu^2 scaling
/  modified6:  February, 2008 by Robert Harkness
/              Conditional calls to MPI_Barrier to force MPICH progress
/  modified7:  October, 2009 by Ji-hoon Kim
/              Added particle splitter routine
/
/  PURPOSE:
/    This routine is responsible for the evolution of the grid hierarchy.
/    It assumes the hierarchy is already constructed and the grids
/    initialized.  The routine then loops over time until one of the
/    stopping criteria is reached.  This routine also handles data dumps,
/    history dumps and restarts dumps (although the later two have not
/    yet been implemented).
/
************************************************************************/
#include "preincludes.h"
 
#ifdef USE_MPI
#include <mpi.h>
#endif
 
#include <stdio.h>

#include "EnzoTiming.h" 
#include "performance.h"
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
#include "communication.h"
#include "CommunicationUtilities.h"
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif
 
// function prototypes
 
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior
#ifdef TRANSFER
		, ImplicitProblemABC *ImplicitSolver
#endif
    ,SiblingGridList *SiblingGridListStorage[]
		);

int EvolveLevel_RK2(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
                    int level, float dtLevelAbove, ExternalBoundary *Exterior, 
#ifdef TRANSFER
		    ImplicitProblemABC *ImplicitSolver, 
#endif
		    FLOAT dt0 ,SiblingGridList *SiblingGridListStorage[]);

int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, 
#ifdef TRANSFER
		 ImplicitProblemABC *ImplicitSolver,
#endif		 
		 FLOAT WriteTime = -1);

int Group_WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, 
#ifdef TRANSFER
		 ImplicitProblemABC *ImplicitSolver,
#endif		 
		 FLOAT WriteTime = -1,
		 int RestartDump = FALSE);

int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData,
			 LevelHierarchyEntry *LevelArray[], int level);
int TestGravityCheckResults(LevelHierarchyEntry *LevelArray[]);
int TestGravitySphereCheckResults(LevelHierarchyEntry *LevelArray[]);
int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData,
		   ExternalBoundary *Exterior, 
#ifdef TRANSFER
		   ImplicitProblemABC *ImplicitSolver,
#endif		 
		   int Restart = FALSE);
int CheckForTimeAction(LevelHierarchyEntry *LevelArray[],
		       TopGridData &MetaData);
int CheckForResubmit(TopGridData &MetaData, int &Stop);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int OutputLevelInformation(FILE *fptr, TopGridData &MetaData,
			   LevelHierarchyEntry *LevelArray[]);
int PrepareGravitatingMassField(HierarchyEntry *Grid, TopGridData *MetaData,
				LevelHierarchyEntry *LevelArray[], int level);
int ReduceFragmentation(HierarchyEntry &TopGrid, TopGridData &MetaData,
			ExternalBoundary *Exterior,
			LevelHierarchyEntry *LevelArray[]);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);
double ReturnWallTime(void);
int Enzo_Dims_create(int nnodes, int ndims, int *dims);
int FOF(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], 
	int WroteData, int FOFOnly=FALSE);
int StarParticleCountOnly(LevelHierarchyEntry *LevelArray[]);
int CommunicationLoadBalanceRootGrids(LevelHierarchyEntry *LevelArray[], 
				      int TopGridRank, int CycleNumber);
int CommunicationBroadcastValue(Eint32 *Value, int BroadcastProcessor);
int CommunicationBroadcastValue(Eint64 *Value, int BroadcastProcessor);
int ParticleSplitter(LevelHierarchyEntry *LevelArray[], int ThisLevel,
		     TopGridData *MetaData); 
int MagneticFieldResetter(LevelHierarchyEntry *LevelArray[], int ThisLevel,
			  TopGridData *MetaData); 
void PrintMemoryUsage(char *str);
int SetEvolveRefineRegion(FLOAT time);

#ifdef MEM_TRACE
Eint64 mused(void);
#endif
#ifdef USE_PYTHON
int CallPython(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level, int from_topgrid);
#endif

 
 
#define NO_REDUCE_FRAGMENTATION
 


 
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &MetaData,
                    ExternalBoundary *Exterior,
#ifdef TRANSFER
		    ImplicitProblemABC *ImplicitSolver,
#endif
		    LevelHierarchyEntry *LevelArray[], float Initialdt)
{
 
  float dt;
 
  int i, dim, Stop = FALSE;
  int StoppedByOutput = FALSE;
  int Restart = FALSE;
  double tlev0, tlev1, treb0, treb1, tloop0, tloop1, tentry, texit;
  LevelHierarchyEntry *Temp;
  double LastCPUTime;

  LCAPERF_BEGIN("EL");
  LCAPERF_START("EvolveHierarchy");

#ifdef USE_MPI
  tentry = MPI_Wtime();
#endif
 
  if (MetaData.Time        >= MetaData.StopTime ) Stop = TRUE;
  if (MetaData.CycleNumber >= MetaData.StopCycle) Stop = TRUE;
  MetaData.StartCPUTime = LastCPUTime = ReturnWallTime();
  MetaData.LastCycleCPUTime = 0.0;

  // Reset CPUTime, if it's very large (absolute UNIX time), which
  // was the default from before.
  if (MetaData.CPUTime > 1e2*MetaData.StopCPUTime)
    MetaData.CPUTime = 0.0;
 
  /* Attach RandomForcingFields to BaryonFields temporarily to apply BCs */
 
  if (RandomForcing) { //AK
    Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->AppendForcingToBaryonFields();
      Temp = Temp->NextGridThisLevel;
    }
    Exterior->AppendForcingToBaryonFields();
  }
 
  /* Set top grid boundary conditions. */

  Temp = LevelArray[0];

  PrintMemoryUsage("Enter EH");

#ifdef FORCE_MSG_PROGRESS
  CommunicationBarrier();
#endif

  CommunicationReceiveIndex = 0;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;

  while (Temp != NULL) {
    if (Temp->GridData->SetExternalBoundaryValues(Exterior) == FAIL) {
      //      ENZO_FAIL("Error in grid->SetExternalBoundaryValues.\n");
      Exterior->Prepare(Temp->GridData);

    }
    if (CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, 0)
	== FAIL)
      ENZO_FAIL("Error in CopyOverlappingZones.");
    Temp = Temp->NextGridThisLevel;
  }
  
  CommunicationDirection = COMMUNICATION_SEND;

  Temp = LevelArray[0];
  while (Temp != NULL) {
    if (CopyOverlappingZones(Temp->GridData, &MetaData, LevelArray, 0)
	== FAIL)
      ENZO_FAIL("Error in CopyOverlappingZones.");
    Temp = Temp->NextGridThisLevel;
  }

#ifdef FORCE_MSG_PROGRESS 
  CommunicationBarrier();
#endif

  CommunicationReceiveHandler();

#ifdef FORCE_MSG_PROGRESS
  CommunicationBarrier();
#endif

  PrintMemoryUsage("Bdry set");
 
  /* Remove RandomForcingFields from BaryonFields when BCs are set. */
 
  if (RandomForcing) { //AK
    LevelHierarchyEntry *Temp = LevelArray[0];
    while (Temp != NULL) {
      Temp->GridData->DetachForcingFromBaryonFields();
      Temp = Temp->NextGridThisLevel;
    }
    Exterior->DetachForcingFromBaryonFields();
  }
 
  /* Check for output. */
 
  CheckForOutput(&TopGrid, MetaData, Exterior, 
#ifdef TRANSFER
		 ImplicitSolver,
#endif		 
		 Restart);

  PrintMemoryUsage("Output");
 
  /* Compute the acceleration field so ComputeTimeStep can find dtAccel.
     (Actually, this is a huge pain-in-the-ass, so only do it if the
      problem really requires it). */
 
/*
  if (ProblemType == 21) {
    PrepareGravitatingMassField(&TopGrid, &MetaData, LevelArray, 0);
    ComputePotentialFieldLevelZero(&MetaData, Grids, NumberOfGrids);
    TopGrid.GridData->ComputeAccelerationField(GRIDS);
  }
*/


  /* Do the first grid regeneration. */
 
  if(CheckpointRestart == FALSE) {
    RebuildHierarchy(&MetaData, LevelArray, 0);
  }

  PrintMemoryUsage("1st rebuild");
 
  /* Particle Splitter. Split particles into 13 (=1+12) child particles */
  
  if (MetaData.FirstTimestepAfterRestart == TRUE &&
      ParticleSplitterIterations > 0)
    ParticleSplitter(LevelArray, 0, &MetaData);

  /* Reset magnetic fields if requested. */
  
  if (MetaData.FirstTimestepAfterRestart == TRUE &&
      ResetMagneticField == TRUE)
    MagneticFieldResetter(LevelArray, 0, &MetaData);

  /* Open the OutputLevelInformation file. */
 
  FILE *LevelInfofptr;
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    LevelInfofptr = fopen("OutputLevelInformation.out", "w");
    fclose(LevelInfofptr);
  }

  /* For top-level timestepping with radiative star particles, we want
     to restrict the timesteps.  Collect info here. */

  StarParticleCountOnly(LevelArray);
 
#ifdef USE_LCAPERF
  Eint32 lcaperf_iter;
#endif

  LCAPERF_STOP("EvolveHierarchy");
  LCAPERF_END("EH");

  SiblingGridList *SiblingGridListStorage[MAX_DEPTH_OF_HIERARCHY];
  for( int level=0; level < MAX_DEPTH_OF_HIERARCHY; level++ ){
    SiblingGridListStorage[level] = NULL;
  }


  /* ====== MAIN LOOP ===== */

  bool FirstLoop = true;
  while (!Stop) {

  TIMER_START("Total");

#ifdef USE_LCAPERF
    lcaperf_iter = MetaData.CycleNumber;
    static bool isFirstCall = true;
    if ((lcaperf_iter % LCAPERF_DUMP_FREQUENCY)==0 || isFirstCall) lcaperf.begin("EL");
    isFirstCall = false;
    lcaperf.attribute ("timestep",&lcaperf_iter, LCAPERF_INT);
    lcaperf.start("EL");
#endif

#ifdef USE_MPI
    tloop0 = MPI_Wtime();
#endif

#ifdef MEM_TRACE
    fprintf(memtracePtr, "==== CYCLE %"ISYM" ====\n", MetaData.CycleNumber);
#endif    
    PrintMemoryUsage("Top");

    /* Load balance the root grids if this isn't the initial call */

    if ((CheckpointRestart == FALSE) && (!FirstLoop))
      CommunicationLoadBalanceRootGrids(LevelArray, MetaData.TopGridRank, 
					MetaData.CycleNumber);

    /* Output level information to log file. */
 
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      LevelInfofptr = fopen("OutputLevelInformation.out", "a");
      if (LevelInfofptr == NULL)
        ENZO_FAIL("Can't open OutputLevelInformation.out!");
    }

    // OutputLevelInformation() only needs to be called by all processors
    // when lcaperf is enabled.

    OutputLevelInformation(LevelInfofptr, MetaData, LevelArray);

    if (MyProcessorNumber == ROOT_PROCESSOR)
      fclose(LevelInfofptr);
 
    /* Compute minimum timestep on the top level. */
 
    float dtProc   = huge_number;
    Temp = LevelArray[0];
 
    // Start skipping
    if(CheckpointRestart == FALSE) {
      while (Temp != NULL) {
        float dtProcTemp = Temp->GridData->ComputeTimeStep();
        dtProc = min(dtProc, dtProcTemp);
        Temp = Temp->NextGridThisLevel;
      }

      dt = RootGridCourantSafetyNumber*CommunicationMinValue(dtProc);
      dt = min(MetaData.MaximumTopGridTimeStep, dt);

      if (debug) fprintf(stderr, "dt, Initialdt: %g %g \n", dt, Initialdt);
      if (Initialdt != 0) {
      
	dt = min(dt, Initialdt);
	if (debug) fprintf(stderr, "dt, Initialdt: %g %g \n", dt, Initialdt);
        Initialdt = 0;
      }

      /* Make sure timestep doesn't go past an output. */

      if (ComovingCoordinates)
        for (i = 0; i < MAX_NUMBER_OF_OUTPUT_REDSHIFTS; i++)
          if (CosmologyOutputRedshift[i] != -1)
            dt = min(1.0001*(CosmologyOutputRedshiftTime[i]-MetaData.Time), dt);
      for (i = 0; i < MAX_TIME_ACTIONS; i++)
        if (TimeActionTime[i] > 0 && TimeActionType[i] > 0)
          dt = min(1.0001*(TimeActionTime[i] - MetaData.Time), dt);
      if (MetaData.dtDataDump > 0.0) {
        while (MetaData.TimeLastDataDump+MetaData.dtDataDump < MetaData.Time)
          MetaData.TimeLastDataDump += MetaData.dtDataDump;
        dt = min(1.0001*(MetaData.TimeLastDataDump + MetaData.dtDataDump -
              MetaData.Time), dt);
      }

      /* Set the time step.  If it will cause Time += dt > StopTime, then
         set dt = StopTime - Time */

      dt = min(MetaData.StopTime - MetaData.Time, dt);
    } else { 
      dt = dtThisLevel[0]; 
    }

    /* Set the time step.  If it will cause Time += dt > StopTime, then
       set dt = StopTime - Time */
 
    dt = min(MetaData.StopTime - MetaData.Time, dt);
    Temp = LevelArray[0];
    // Stop skipping

    // Set dt from stored in CheckpointRestart
    while (Temp != NULL) {
      Temp->GridData->SetTimeStep(dt);
      Temp = Temp->NextGridThisLevel;
    }
 
#ifdef CONFIG_TESTING
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      printf("enzo-test: MetaData.CycleNumber %"ISYM"\n", MetaData.CycleNumber);
      printf("enzo-test: dt %.14g\n",dt);
      printf("enzo-test: MetaData.Time %"GOUTSYM"\n", MetaData.Time);
      fflush(stdout);
    }
#endif

    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "TopGrid dt = %"ESYM"     time = %"GOUTSYM"    cycle = %"ISYM,
	     dt, MetaData.Time, MetaData.CycleNumber);

      if (ComovingCoordinates) {
	FLOAT a, dadt;
	CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
	fprintf(stderr, "    z = %"GOUTSYM, (1 + InitialRedshift)/a - 1);
      }
      fprintf(stderr, "\n");
    }
    //}
 
    /* Inline halo finder */

    FOF(&MetaData, LevelArray, MetaData.WroteData);

    /* If provided, set RefineRegion from evolving RefineRegion */
    if ((RefineRegionTimeType == 1) || (RefineRegionTimeType == 0)) {
        if (SetEvolveRefineRegion(MetaData.Time) == FAIL) 
	  ENZO_FAIL("Error in SetEvolveRefineRegion.");
    }

    /* Evolve the top grid (and hence the entire hierarchy). */

#ifdef USE_MPI 
    CommunicationBarrier();
    tlev0 = MPI_Wtime();
#endif

    /* Zeroing out the rootgrid Emissivity before EvolveLevel is called 
       so when rootgrid emissivity values are calculated they are put in 
       clean rootgrid array */
#ifdef EMISSIVITY
/*
    if(StarMakerEmissivityField > 0){
      LevelHierarchyEntry *RootTemp;
      RootTemp = LevelArray[0];
      while (RootTemp != NULL) {
	RootTemp->GridData->ClearEmissivity();
	RootTemp = RootTemp->NextGridThisLevel;
      }
    }
*/
#endif
 
    if (HydroMethod == PPM_DirectEuler || HydroMethod == Zeus_Hydro || 
	HydroMethod == PPM_LagrangeRemap || HydroMethod == HydroMethodUndefined ||
	HydroMethod == MHD_Li || HydroMethod == NoHydro ||
	HydroMethod < 0) {
      if (EvolveLevel(&MetaData, LevelArray, 0, dt, Exterior
#ifdef TRANSFER
		      , ImplicitSolver
#endif
          ,SiblingGridListStorage
		      ) == FAIL) {
        if (NumberOfProcessors == 1) {
          fprintf(stderr, "Error in EvolveLevel.\n");
          fprintf(stderr, "--> Dumping data (output number %d).\n",
                  MetaData.DataDumpNumber);
	Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
		     &TopGrid, MetaData, Exterior
#ifdef TRANSFER
		     , ImplicitSolver
#endif		 
		     );
        }
        return FAIL;
      }
    } else {
      if (HydroMethod == HD_RK || HydroMethod == MHD_RK)
	if (EvolveLevel_RK2(&MetaData, LevelArray, 0, dt, Exterior, 
#ifdef TRANSFER
			    ImplicitSolver, 
#endif
			    dt, SiblingGridListStorage) == FAIL) {
	  if (NumberOfProcessors == 1) {
	    fprintf(stderr, "Error in EvolveLevel_RK2.\n");
	    fprintf(stderr, "--> Dumping data (output number %d).\n",
		    MetaData.DataDumpNumber);
	    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
			       &TopGrid, MetaData, Exterior
#ifdef TRANSFER
			       , ImplicitSolver
#endif		 
			       );
	  }
        return FAIL;
      }
    }



#ifdef USE_MPI 
    CommunicationBarrier();
    tlev1 = MPI_Wtime();
#endif
  
    /* Add time and check stopping criteria (steps #21 & #22)
       (note the topgrid is also keeping its own time but this statement will
       keep the two in synch). */
 
    MetaData.Time += dt;
    MetaData.CycleNumber++;
    MetaData.LastCycleCPUTime = ReturnWallTime() - LastCPUTime;
    MetaData.CPUTime += MetaData.LastCycleCPUTime;
    LastCPUTime = ReturnWallTime();

    if (MyProcessorNumber == ROOT_PROCESSOR) {
	
    if (MetaData.Time >= MetaData.StopTime) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on top grid time limit.\n");
      Stop = TRUE;
    } else
    if (MetaData.CycleNumber >= MetaData.StopCycle) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on top grid cycle limit.\n");
      Stop = TRUE;
    } else
    if (MetaData.CPUTime >= MetaData.StopCPUTime) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on CPU time limit.\n");
      Stop = TRUE;
    } else
    if ((ReturnWallTime() - MetaData.StartCPUTime >= MetaData.dtRestartDump &&
	 MetaData.dtRestartDump > 0) ||
	(MetaData.CycleNumber - MetaData.CycleLastRestartDump >= 
	 MetaData.CycleSkipRestartDump &&
	 MetaData.CycleSkipRestartDump > 0)) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping to restart.\n");
      Stop = TRUE;
      Restart = TRUE;
    }
    } // ENDIF ROOT_PROCESSOR

    CommunicationBroadcastValue(&Stop, ROOT_PROCESSOR);
    CommunicationBroadcastValue(&Restart, ROOT_PROCESSOR);

    /* If not restarting, rebuild the grids from level 0. */

#ifdef USE_MPI
    treb0 = MPI_Wtime();
#endif

    PrintMemoryUsage("Pre loop rebuild");
 
    if (ProblemType != 25 && Restart == FALSE)
      RebuildHierarchy(&MetaData, LevelArray, 0);

    PrintMemoryUsage("Post loop rebuild");

#ifdef USE_MPI
    treb1 = MPI_Wtime();
#endif
 
    /* Check for time-actions. */
 
    CheckForTimeAction(LevelArray, MetaData);
 
    /* Check for output. */
 
    CheckForOutput(&TopGrid, MetaData, Exterior, 
#ifdef TRANSFER
		   ImplicitSolver,
#endif		 
		   Restart);

    /* Call inline analysis. */

#ifdef USE_PYTHON
    LCAPERF_START("CallPython");
    CallPython(LevelArray, &MetaData, 0, 1);
    LCAPERF_STOP("CallPython");
#endif

    /* Check for resubmission */
    
    if (!Restart)
      CheckForResubmit(MetaData, Stop);

    /* If stopping, inline halo finder one more time */

    if (Stop && !Restart)
      FOF(&MetaData, LevelArray, TRUE);

    /* Try to cut down on memory fragmentation. */
 
#ifdef REDUCE_FRAGMENTATION
 
    if (MetaData.WroteData && !Stop)
      ReduceFragmentation(TopGrid, MetaData, Exterior, LevelArray);
 
#endif /* REDUCE_FRAGMENTATION */

#ifdef USE_LCAPERF
    lcaperf.stop("EL");
    if (((lcaperf_iter+1) % LCAPERF_DUMP_FREQUENCY)==0) lcaperf.end("EL");
#endif

    PrintMemoryUsage("Bot");

  for ( i = 0; i < MAX_NUMBER_OF_TASKS; i++ ) {
    TaskMemory[i] = -1;
  }

#ifdef MEM_TRACE

  /*
  MPI_Datatype DataTypeInt = (sizeof(Eint64) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  MPI_Arg ThisTask;
  MPI_Arg TaskCount;
  MPI_Arg Count = 1;
  MPI_Arg stat;

  stat = MPI_Comm_size(MPI_COMM_WORLD, &TaskCount);
  stat = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  stat = MPI_Allgather(&MemInUse, Count, DataTypeInt, TaskMemory, Count, DataTypeInt, MPI_COMM_WORLD);

  if (ThisTask == 0 ) {
    for ( i = 0; i < TaskCount; i++) {
      fprintf(stderr, "TaskMemory : Task %"ISYM"  Memory %"ISYM"\n", i, TaskMemory[i]);
    }
  }
*/

#endif /* MEM_TRACE */
#ifdef USE_MPI
    tloop1 = MPI_Wtime();
#endif

    FILE *evlog;

    if (MyProcessorNumber == ROOT_PROCESSOR) {
      evlog = fopen("Evtime", "a");
      fprintf(evlog, "%8"ISYM"  %16.9e  %16.9e  %16.9e\n", MetaData.CycleNumber, tlev1-tlev0, treb1-treb0, tloop1-tloop0);
      fclose(evlog);
    }

#ifdef MEM_TRACE
    Eint64 MemInUse;
    if (MetaData.WroteData) {
      MemInUse = mused();
      MemInUse = CommunicationMaxValue(MemInUse);
      if (MemInUse > MemoryLimit) {
        if (MyProcessorNumber == ROOT_PROCESSOR)
          printf("Stopping due to memory limit.\n");
        Stop = TRUE;
      }
    }
#endif

    TIMER_STOP("Total");
    if ((MetaData.CycleNumber-1) % TimingCycleSkip == 0)
		  TIMER_WRITE(MetaData.CycleNumber);

    FirstLoop = false;
 
    /* If simulation is set to stop after writing a set number of outputs, check that here. */

    if (MetaData.NumberOfOutputsBeforeExit && MetaData.WroteData) {
      MetaData.OutputsLeftBeforeExit--;
      if (MetaData.OutputsLeftBeforeExit <= 0) {
        if (MyProcessorNumber == ROOT_PROCESSOR) {
          fprintf(stderr, "Exiting after writing %"ISYM" datadumps.\n",
                  MetaData.NumberOfOutputsBeforeExit);
        }      
        Stop = TRUE;
        StoppedByOutput = TRUE;
      }
    }

  } // ===== end of main loop ====

#ifdef USE_LCAPERF
  if (((lcaperf_iter+1) % LCAPERF_DUMP_FREQUENCY)!=0) lcaperf.end("EL");
  lcaperf.attribute ("timestep",0, LCAPERF_NULL);
#endif

  MetaData.CPUTime = ReturnWallTime() - MetaData.StartCPUTime;
 
  /* Done, so report on current time, etc. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("Time     = %9"FSYM"   CycleNumber = %6"ISYM"    Wallclock   = %9"FSYM"\n",
	   MetaData.Time, MetaData.CycleNumber, MetaData.CPUTime);
    printf("StopTime = %9"FSYM"   StopCycle   = %6"ISYM"\n",
	   MetaData.StopTime, MetaData.StopCycle);
  }
 
  /* If we are running problem 23, TestGravity, then check the results. */
 
  if (ProblemType == 23)
    TestGravityCheckResults(LevelArray);
  if (ProblemType == 25 && NumberOfProcessors == 0)
    TestGravitySphereCheckResults(LevelArray);
 
  /* if we are doing data dumps, then do one last one */
 
  if ((MetaData.dtDataDump != 0.0 || MetaData.CycleSkipDataDump != 0) &&
      !MetaData.WroteData)
    //#ifdef USE_HDF5_GROUPS
    if (Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
			   &TopGrid, MetaData, Exterior, 
#ifdef TRANSFER
			   ImplicitSolver, 
#endif		 
			   -666) == FAIL)
      ENZO_FAIL("Error in Group_WriteAllData.");
// #else
//     if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
// 		     &TopGrid, MetaData, Exterior, 
//#ifdef TRANSFER
//		     ImplicitSolver, 
//#endif		 
//                   -666) == FAIL) {
//       ENZO_FAIL("Error in WriteAllData.\n");
//     }
// #endif
 
  /* Write a file to indicate that we're finished. */

  FILE *Exit_fptr;
  if (!Restart && !StoppedByOutput && MyProcessorNumber == ROOT_PROCESSOR) {
    if ((Exit_fptr = fopen("RunFinished", "w")) == NULL)
      ENZO_FAIL("Error opening RunFinished.");
    fprintf(Exit_fptr, "Finished on cycle %"ISYM"\n", MetaData.CycleNumber);
    fclose(Exit_fptr);
  }

  if (NumberOfProcessors > 1)

    printf("Communication: processor %"ISYM" CommunicationTime = %"FSYM"\n",
	   MyProcessorNumber, CommunicationTime);
 
  /* done */

#ifdef USE_MPI
  texit = MPI_Wtime();
#endif
 
  return SUCCESS;
}
