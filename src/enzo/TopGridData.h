#ifndef __TopGridData_h_
#define __TopGridData_h_
/***********************************************************************
/
/  TOP GRID DATA STRUCTURE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       February 29th, 2008
/
/  PURPOSE:
/
************************************************************************/
#include "AMRH5writer.h"
struct TopGridData
{

  /* Counters for the TopGrid. */

  int   CycleNumber;         // Number of top grid timestep performed
  int   SubcycleNumber;         // Number of top grid tximestep performed
  FLOAT Time;                // Current problem time
  double CPUTime;            // Current CPU time used
  double StartCPUTime;
  double LastCycleCPUTime;    // CPU time used in the last cycle
  int ResubmitOn;             // Resubmit job after StopCPUTime

  /* Script names for resubmission to queues and restarting to reduce
     memory fragmentation. */

  char *ResubmitCommand;      // Script name for job resubmission

  /* Stopping criteria for TopGrid. */

  FLOAT StopTime;            // time to stop at
  int   StopCycle;           // timestep number to stop at
  int   StopSteps;           // stop after N steps (for heap fragmentation fix)
  float StopCPUTime;         // Maximum CPU time to be used

  /* Time step limit */
  float MaximumTopGridTimeStep; // limit the topgrid time step to be smaller than this

  /* Parameters governing when output is done. */

  float TimeLastRestartDump;  // CPU time of the last restart dump (seconds)
  float dtRestartDump;        // CPU time between restart dumps (0 = never)

  FLOAT TimeLastDataDump;     // Problem time of the last data dump
  FLOAT dtDataDump;           // Problem time between data dumps (0 = never)

  FLOAT TimeLastHistoryDump; // Problem time of the last history (small) dump
  FLOAT dtHistoryDump;       // Problem time between history dumps (0 = never)

  FLOAT TimeLastMovieDump;   // Problem time of the last movie dump
  FLOAT dtMovieDump;          // Problem time between movie dumps (0 = never)

  FLOAT TimeLastTracerParticleDump;   // Problem time of last tracer part dump
  FLOAT dtTracerParticleDump;         // Problem time between dumps (0 = never)

  FLOAT TimeLastInterpolatedDataDump;     // Problem time of the last interpolated data dump
  FLOAT dtInterpolatedDataDump;           // Problem time between interpolated data dumps (0 = never)

  FLOAT NewMovieLeftEdge[MAX_DIMENSION];  // region for seq. movie output
  FLOAT NewMovieRightEdge[MAX_DIMENSION];

  int CycleLastRestartDump;  // Cycle of the last restart dump (seconds)
  int CycleSkipRestartDump;  // Cycles between restart dumps (0 = never)

  int CycleLastDataDump;     // Cycle of the last data dump
  int CycleSkipDataDump;     // Number of cycles between data dumps (0 = never)

  int SubcycleLastDataDump;     // SubCycle of the last data dump
  int SubcycleSkipDataDump;     // Number of subcycles between data dumps (0 = never)

  int CycleLastHistoryDump;  // Cycle of the last history (small) dump
  int CycleSkipHistoryDump;  // Cycles between history dumps (0 = never)
  int CycleSkipGlobalDataDump;//AK Cycles between global data dumps (0 = never)

  int OutputFirstTimeAtLevel; // Outputs when a new level is generated
  int StopFirstTimeAtLevel;   // Stops when this level is first reached

  int NumberOfOutputsBeforeExit; // Number of datadumps to write before exiting (0 = never)
  int OutputsLeftBeforeExit;     // Number of datadumps actually left before exiting

  int WroteData;              // Flag if data dump written this iteration.

  /* Parameters governing output names. */

  int RestartDumpNumber;        // number appended to end of restart dump name
  int DataDumpNumber;           // number appended to end of data dump name
  int HistoryDumpNumber;        // number appended to end of history dump
  int MovieDumpNumber;          // number appended to end of movie dump name
  int TracerParticleDumpNumber; // number of dump

  char *RestartDumpName;         // restart dump base name
  char *DataDumpName;            // data dump base name
  char *HistoryDumpName;         // history dump base name
  char *MovieDumpName;           // movie dump base name
  char *TracerParticleDumpName;  // movie dump name
  char *RedshiftDumpName;        // redshift dump base name

  char *RestartDumpDir;         // restart dump directory name
  char *DataDumpDir;            // data dump directory name
  char *HistoryDumpDir;         // history dump directory name
  char *MovieDumpDir;           // movie dump directory name
  char *TracerParticleDumpDir;  // tracer particle dump directory name
  char *RedshiftDumpDir;        // redshift dump directory name
  char * ExtraDumpDir;          // Another directory name.  For the extra dumps.
  char * ExtraDumpName;          // Another directory name.  For the extra dumps.

  char *LocalDir;               // local disk directory name
  char *GlobalDir;              // global disk directory name

  char *MetaDataIdentifier;     // A name (string) that will be persisted between datasets
  char *SimulationUUID;         // Unique identifier for the simulation
  char *RestartDatasetUUID;     // Identifier of the dataset restarting from
  char *InitialConditionsUUID;  // Identifier of the initial conditions used

  /* TopGrid Parameters governing hierarchy */

  int StaticHierarchy;     // TRUE for static mesh refinement

  /* Some grid defining data
     These are here out of convenience, the real ones are in the grids. */

  int TopGridRank;
  int TopGridDims[MAX_DIMENSION];
  boundary_type  LeftFaceBoundaryCondition[MAX_DIMENSION],
                RightFaceBoundaryCondition[MAX_DIMENSION];
  char *BoundaryConditionName;

  /* Gravity data -- used only for top grid potential field solve */

  gravity_boundary_type GravityBoundary;

#ifdef TRANSFER
  /* Implicit solver data */
  char *RadHydroParameterFname;
  FLOAT FLDTime;
  float dtFLD;
#endif

  /* Particle and Particle boundary data. (real one in ExternalBoundary). */

  boundary_type ParticleBoundaryType;
  PINT          NumberOfParticles;

  /* Hydro Parameters.  
     These are here out of convenience, the real ones are in the grids. */

  float  CourantSafetyNumber;                       // Hydro parameter
  int    PPMFlatteningParameter;                    // PPM parameter
  int    PPMDiffusionParameter;                     // PPM parameter
  int    PPMSteepeningParameter;                    // PPM parameter

  AMRHDF5Writer AmiraGrid;
  int FirstTimestepAfterRestart;
  int MovieTimestepCounter;
  float GlobalMaximumkphIfront;

};

#endif
