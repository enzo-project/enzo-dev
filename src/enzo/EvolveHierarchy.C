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
 
#ifdef USE_MPI
#include <mpi.h>
#endif
 
#include <stdio.h>
 
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
 
// function prototypes
 
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior);

int EvolveLevel_RK2(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
                    int level, float dtLevelAbove, ExternalBoundary *Exterior, FLOAT dt0);

int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, FLOAT WriteTime = -1);

int Group_WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, FLOAT WriteTime = -1,
         int RestartDump = FALSE);

int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData,
			 LevelHierarchyEntry *LevelArray[], int level);
int TestGravityCheckResults(LevelHierarchyEntry *LevelArray[]);
int TestGravitySphereCheckResults(LevelHierarchyEntry *LevelArray[]);
int CheckForOutput(HierarchyEntry *TopGrid, TopGridData &MetaData,
		   ExternalBoundary *Exterior, int &WroteData);
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
	int WroteData);
int StarParticleCountOnly(LevelHierarchyEntry *LevelArray[]);
int CommunicationLoadBalanceRootGrids(LevelHierarchyEntry *LevelArray[], 
				      int TopGridRank, int CycleNumber);
int ParticleSplitter(LevelHierarchyEntry *LevelArray[], int ThisLevel,
		     TopGridData *MetaData); 

#ifdef USE_PYTHON
int CallPython();
#endif

#ifdef MEM_TRACE
Eint64 mused(void);
#endif
 
 
#define NO_REDUCE_FRAGMENTATION
 
 
 
 
int EvolveHierarchy(HierarchyEntry &TopGrid, TopGridData &MetaData,
                    ExternalBoundary *Exterior,
		    LevelHierarchyEntry *LevelArray[], float Initialdt)
{
 
  float dt;
 
  int i, dim, Stop = FALSE, WroteData;
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
  MetaData.StartCPUTime = MetaData.CPUTime = LastCPUTime = ReturnWallTime();
  MetaData.LastCycleCPUTime = 0.0;
 
#ifdef MEM_TRACE
  Eint64 MemInUse;
#endif
 
  /* Double-check if the topgrid is evenly divided if we're using the
     optimized version of CommunicationTransferParticles. */

#ifdef OPTIMIZED_CTP
  int NumberOfGrids = 0, Layout[MAX_DIMENSION] = {1,1,1};
  Temp = LevelArray[0];
  while (Temp != NULL) {
    NumberOfGrids++;
    Temp = Temp->NextGridThisLevel;
  }
  Enzo_Dims_create(NumberOfGrids, MetaData.TopGridRank, Layout);
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    if (MetaData.TopGridDims[dim] % Layout[MAX_DIMENSION-1-dim] != 0) {
      if (debug)
	fprintf(stderr, "ERROR: "
		"\tTo use the optimized CommunicationTransferParticles routine,\n"
		"\tthe top grid must be split evenly, "
		"i.e. mod(Dims[i], Layout[i]) != 0\n"
		"\t==> dimension %"ISYM": Dims = %"ISYM", Layout = %"ISYM"\n",
		dim, MetaData.TopGridDims[dim], Layout[MAX_DIMENSION-1-dim]);
      ENZO_FAIL("");
    }
#endif /* OPTIMIZED_CTP */

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

#ifdef MEM_TRACE
  MemInUse = mused();
  fprintf(memtracePtr, "Enter EH %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif

#ifdef FORCE_MSG_PROGRESS
  CommunicationBarrier();
#endif

  CommunicationReceiveIndex = 0;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;

  while (Temp != NULL) {
    if (Temp->GridData->SetExternalBoundaryValues(Exterior) == FAIL) {
      fprintf(stderr, "Error in grid->SetExternalBoundaryValues.\n");
      //      ENZO_FAIL("");
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
 
#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Bdry set %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif
 
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
 
  CheckForOutput(&TopGrid, MetaData, Exterior, WroteData);

#ifdef MEM_TRACE
  MemInUse = mused();
  fprintf(memtracePtr, "Output %8"ISYM"  %16"ISYM" \n", 
	  MetaData.CycleNumber, MemInUse);
#endif
 
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

  /* Particle Splitter. Split the particles into 13 (=1+12) children 
     particles */
  
  if (MetaData.FirstTimestepAfterRestart == TRUE &&
      ParticleSplitterIterations > 0)
    ParticleSplitter(LevelArray, 0, &MetaData);
 
  



#ifdef MEM_TRACE
  MemInUse = mused();
  fprintf(memtracePtr, "1st rebuild %8"ISYM"  %16"ISYM" \n", 
	  MetaData.CycleNumber, MemInUse);
#endif
 
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

  /* ====== MAIN LOOP ===== */

  bool FirstLoop = true;
  while (!Stop) {

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
    MemInUse = mused();
    fprintf(memtracePtr, "Top %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif

    /* Load balance the root grids if this isn't the initial call */

    if ((CheckpointRestart == FALSE) && (!FirstLoop))
      CommunicationLoadBalanceRootGrids(LevelArray, MetaData.TopGridRank, 
					MetaData.CycleNumber);

    /* Output level information to log file. */
 
    if (MyProcessorNumber == ROOT_PROCESSOR)
      LevelInfofptr = fopen("OutputLevelInformation.out", "a");

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
        dtProc = min(dtProc, Temp->GridData->ComputeTimeStep());
        Temp = Temp->NextGridThisLevel;
      }

      dt = RootGridCourantSafetyNumber*CommunicationMinValue(dtProc);

    dt = RootGridCourantSafetyNumber*CommunicationMinValue(dtProc);
    dt = min(MetaData.MaximumTopGridTimeStep, dt);

    if (debug) fprintf(stderr, "dt, Initialdt: %g %g \n", dt, Initialdt);
    if (Initialdt != 0) {
      
      dt = min(dt, Initialdt);
      if (debug) fprintf(stderr, "dt, Initialdt: %g %g \n", dt, Initialdt);
#ifdef TRANSFER
        dtPhoton = dt;
#endif
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
      printf("TopGrid dt = %"ESYM"     time = %"GOUTSYM"    cycle = %"ISYM,
	     dt, MetaData.Time, MetaData.CycleNumber);

      if (ComovingCoordinates) {
	FLOAT a, dadt;
	CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
	printf("    z = %"GOUTSYM, (1 + InitialRedshift)/a - 1);
      }
      printf("\n");
    }
    //}
 
    /* Inline halo finder */

    FOF(&MetaData, LevelArray, WroteData);

    /* Evolve the top grid (and hence the entire hierarchy). */

#ifdef USE_MPI 
    CommunicationBarrier();
    tlev0 = MPI_Wtime();
#endif
 
    if (HydroMethod == PPM_DirectEuler || HydroMethod == Zeus_Hydro || 
	HydroMethod == PPM_LagrangeRemap || HydroMethod == HydroMethodUndefined ||
	HydroMethod < 0) {
      if (EvolveLevel(&MetaData, LevelArray, 0, dt, Exterior) == FAIL) {
        if (NumberOfProcessors == 1) {
          fprintf(stderr, "Error in EvolveLevel.\n");
          fprintf(stderr, "--> Dumping data (output number %d).\n",
                  MetaData.DataDumpNumber);
	Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
		     &TopGrid, MetaData, Exterior);
        }
        return FAIL;
      }
    } else {
      if (HydroMethod == HD_RK || HydroMethod == MHD_RK)
	if (EvolveLevel_RK2(&MetaData, LevelArray, 0, dt, Exterior, dt) == FAIL) {
	  if (NumberOfProcessors == 1) {
	    fprintf(stderr, "Error in EvolveLevel_RK2.\n");
	    fprintf(stderr, "--> Dumping data (output number %d).\n",
		    MetaData.DataDumpNumber);
	    Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
			       &TopGrid, MetaData, Exterior);
	  }
        return FAIL;
      }
    }



#ifdef USE_MPI 
    CommunicationBarrier();
    tlev1 = MPI_Wtime();
#endif
 
    /* Rebuild the grids from level 0. */

#ifdef USE_MPI
    treb0 = MPI_Wtime();
#endif

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Pre loop rebuild %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif
 
    if (ProblemType != 25)
      if (RebuildHierarchy(&MetaData, LevelArray, 0) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	ENZO_FAIL("");
      }

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Post loop rebuild %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif

#ifdef USE_MPI
    treb1 = MPI_Wtime();
#endif
 
    /* Add time and check stopping criteria (steps #21 & #22)
       (note the topgrid is also keeping its own time but this statement will
       keep the two in synch). */
 
    MetaData.Time += dt;
    MetaData.CycleNumber++;
    MetaData.LastCycleCPUTime = ReturnWallTime() - LastCPUTime;
    LastCPUTime = ReturnWallTime();
	
    if (MetaData.Time >= MetaData.StopTime) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on top grid time limit.\n");
      Stop = TRUE;
    }
    if (MetaData.CycleNumber >= MetaData.StopCycle) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on top grid cycle limit.\n");
      Stop = TRUE;
    }
    if (ReturnWallTime() - MetaData.StartCPUTime >= MetaData.StopCPUTime) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("Stopping on CPU time limit.\n");
      Stop = TRUE;
    }
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
 
    /* Check for time-actions. */
 
    CheckForTimeAction(LevelArray, MetaData);
 
    /* Check for output. */
 
    CheckForOutput(&TopGrid, MetaData, Exterior, WroteData);

    /* Check for resubmission */
    
    if (!Restart)
      CheckForResubmit(MetaData, Stop);

    /* If stopping, inline halo finder one more time */

    if (Stop && !Restart)
      FOF(&MetaData, LevelArray, TRUE);

    /* Try to cut down on memory fragmentation. */
 
#ifdef REDUCE_FRAGMENTATION
 
    if (WroteData && !Stop)
      ReduceFragmentation(TopGrid, MetaData, Exterior, LevelArray);
 
#endif /* REDUCE_FRAGMENTATION */

#ifdef USE_LCAPERF
    lcaperf.stop("EL");
    if (((lcaperf_iter+1) % LCAPERF_DUMP_FREQUENCY)==0) lcaperf.end("EL");
#endif

#ifdef MEM_TRACE
    MemInUse = mused();
    fprintf(memtracePtr, "Bot %8"ISYM"  %16"ISYM" \n", MetaData.CycleNumber, MemInUse);
#endif

  for ( i = 0; i < MAX_NUMBER_OF_TASKS; i++ ) {
    TaskMemory[i] = -1;
  }

#ifdef MEM_TRACE

  MPI_Datatype DataTypeInt = (sizeof(Eint64) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  MPI_Arg ThisTask;
  MPI_Arg TaskCount;
  MPI_Arg Count = 1;
  MPI_Arg stat;

  stat = MPI_Comm_size(MPI_COMM_WORLD, &TaskCount);
  stat = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  stat = MPI_Allgather(&MemInUse, Count, DataTypeInt, TaskMemory, Count, DataTypeInt, MPI_COMM_WORLD);

/*
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
    if (WroteData) {
      if (MemInUse > MemoryLimit) {
        if (MyProcessorNumber == ROOT_PROCESSOR)
          printf("Stopping due to memory limit.\n");
        Stop = TRUE;
      }
    }
#endif

    FirstLoop = false;
 
  } // ===== end of main loop ====
 
#ifdef USE_LCAPERF
  if (((lcaperf_iter+1) % LCAPERF_DUMP_FREQUENCY)!=0) lcaperf.end("EL");
  lcaperf.attribute ("timestep",0, LCAPERF_NULL);
#endif

#ifdef USE_MPI
  MetaData.CPUTime = MPI_Wtime() - MetaData.CPUTime;
#endif
 
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
      !WroteData)
    //#ifdef USE_HDF5_GROUPS
    if (Group_WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
			   &TopGrid, MetaData, Exterior, -666) == FAIL)
      ENZO_FAIL("Error in Group_WriteAllData.");
// #else
//     if (WriteAllData(MetaData.DataDumpName, MetaData.DataDumpNumber,
// 		     &TopGrid, MetaData, Exterior, -666) == FAIL) {
//       fprintf(stderr, "Error in WriteAllData.\n");
//       ENZO_FAIL("");
//     }
// #endif
 
  /* Write a file to indicate that we're finished. */

  FILE *Exit_fptr;
  if (!Restart && MyProcessorNumber == ROOT_PROCESSOR) {
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
