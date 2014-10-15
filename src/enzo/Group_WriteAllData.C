/***********************************************************************
/
/  WRITE OUT ALL THE DATA (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       February, 2004
/              October, 2004
/              Direct or indirect SRB driver
/              Local or Global file system, or cwd
/              July , 2006
/              Assemble group file in-core
/  modified2:  Robert Harkness
/  date:       April 2008
/  modified3:  Michael Kuhlen, October 2010, HDF5 hierarchy
/
/  PURPOSE:
/
************************************************************************/
#define SYSCALL
 
// This function writes out the data hierarchy (TopGrid), the External
//   Boundary (Exterior), the TopGridData, and the global_data.
#include "preincludes.h"
#include "EnzoTiming.h" 

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <hdf5.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "h5utilities.h"

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
#include "CosmologyParameters.h"
#include "CommunicationUtilities.h"
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif

void my_exit(int status);
 
// HDF5 function prototypes
 
// function prototypes

int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

int SysMkdir(char *startdir, char *directory);
 
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int WriteDataCubes(HierarchyEntry *TopGrid, int TDdims[], char *gridbasename, int &GridID, FLOAT WriteTime);
int Group_WriteDataHierarchy(FILE *fptr, TopGridData &MetaData, HierarchyEntry *TopGrid,
		       char *gridbasename, int &GridID, FLOAT WriteTime, hid_t file_id,
               int CheckpointDump = FALSE);
int WriteHDF5HierarchyFile(char *base_name, HierarchyEntry *TopGrid, TopGridData MetaData, LevelHierarchyEntry *LevelArray[]);
int WriteMemoryMap(FILE *fptr, HierarchyEntry *TopGrid,
		   char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteConfigure(FILE *optr);
int WriteTaskMap(FILE *fptr, HierarchyEntry *TopGrid,
		 char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteParameterFile(FILE *fptr, TopGridData &MetaData, char *Filename);
int WriteStarParticleData(FILE *fptr, TopGridData &MetaData);
int WriteRadiationData(FILE *fptr);
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int CommunicationCombineGrids(HierarchyEntry *OldHierarchy,
			      HierarchyEntry **NewHierarchyPointer,
			      FLOAT WriteTime, int CheckpointDump);
void DeleteGridHierarchy(HierarchyEntry *GridEntry);
void ContinueExecution(void);
int CreateSmoothedDarkMatterFields(TopGridData &MetaData, HierarchyEntry *TopGrid);
 
int CreateGriddedStarParticleFields(TopGridData &MetaData, HierarchyEntry *TopGrid); 

#ifndef FAST_SIB
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          int level, TopGridData *MetaData,
                          ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#endif 


#ifdef TRANSFER
extern char RTSuffix[];
#endif
extern char BCSuffix[];
extern char GridSuffix[];
extern char HierarchySuffix[];
extern char hdfsuffix[];
extern char RadiationSuffix[];
extern char TaskMapSuffix[];
extern char MemoryMapSuffix[];
extern char ForcingSuffix[]; // WS
extern char ConfigureSuffix[];

char CPUSuffix[]       = ".cpu";
char BHierarchySuffix[] = ".harrays";
 
extern char LastFileNameWritten[MAX_LINE_LENGTH];
 
 
 
 
int Group_WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, 
#ifdef TRANSFER
		 ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1, 
		 int CheckpointDump = FALSE)
{

  TIMER_START("Group_WriteAllData");

  char id[MAX_CYCLE_TAG_SIZE], *cptr, name[MAX_LINE_LENGTH];
  char dumpdirname[MAX_LINE_LENGTH];
  char dumpdirroot[MAX_LINE_LENGTH];
  char unixcommand[MAX_LINE_LENGTH];
  char gridbasename[MAX_LINE_LENGTH];
  char hierarchyname[MAX_LINE_LENGTH];
  char bhierarchyname[MAX_LINE_LENGTH];
  char radiationname[MAX_LINE_LENGTH];
  char taskmapname[MAX_LINE_LENGTH];
  char memorymapname[MAX_LINE_LENGTH];
  char configurename[MAX_LINE_LENGTH];
  char groupfilename[MAX_LINE_LENGTH];
  char forcingname[MAX_LINE_LENGTH]; // WS
 
  int unixresult;
  int status;
  int local, global;
  int file_status;
  int ii, pe, nn;
  double twrite0, twrite1;

  char pid[MAX_TASK_TAG_SIZE];
 
  FILE *fptr;
  FILE *sptr;
  FILE *gptr;
  FILE *tptr;
  FILE *mptr;
  FILE *optr;
 
  hid_t       file_id;
  hid_t       file_acc_template;
  size_t      memory_increment; // in bytes
  hbool_t     dump_flag;
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int GridID = 1;
  int GridJD = 1;
  int GridKD = 1;
  int GridLD = 1;
 
  /* If this is an interpolated time step, then temporary replace  the time
     in MetaData.  Note:  Modified 6 Feb 2006 to fix interpolated  data outputs. */

  FLOAT SavedTime = MetaData.Time;
  MetaData.Time = ((WriteTime < 0) || (CheckpointDump == TRUE)) ? MetaData.Time : WriteTime;

  // Global or local filesystem?
 
  local = 0;
  global = 0;
 
  if (MetaData.LocalDir != NULL)
  {
     local = 1;
     strcpy(dumpdirroot, MetaData.LocalDir);
     // fprintf(stdout, "XXXX local dir: %s\n", MetaData.LocalDir);
     // Create on node - locking?
  }
 
  if (MetaData.GlobalDir != NULL)
  {
     global = 1;
     strcpy(dumpdirroot, MetaData.GlobalDir);
     // fprintf(stdout, "XXXX global dir: %s\n", MetaData.GlobalDir);
     // Create on task 0 only
  }
  
  if (( local == 1) && (global == 1))
    fprintf(stdout, "Local AND Global !!\n");
 
  // Create main name
 
  if (ComovingCoordinates && (cptr = strstr(name, "RRRR"))) {
    FLOAT a, dadt;
    CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
    sprintf(cptr, "%"CYCLE_TAG_FORMAT""ISYM, nint(100*((1 + InitialRedshift)/a - 1)));
  } else {
 
    sprintf(id, "%"CYCLE_TAG_FORMAT""ISYM, filenumber);
    sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);

    /******************** Extra Output ********************/

    if ( (cptr = strstr(basename, MetaData.ExtraDumpName)) ) {
 
      if (MetaData.ExtraDumpDir != NULL)
      {
        if (MetaData.LocalDir != NULL) {
          // Local fs
          strcpy(dumpdirname, MetaData.LocalDir);
          strcat(dumpdirname, "/");
          strcat(dumpdirname, MetaData.ExtraDumpDir);
          strcat(dumpdirname, id);
  
          // Create once per node...
#ifdef USE_NODE_LOCAL
            strcat(dumpdirname, "/mpi");
            strcat(dumpdirname, pid);
#endif /* USE_NODE_LOCAL */
 
          strcpy(name, dumpdirname);
          strcat(name, "/");
          strcat(name, basename);
        } // if LocalDir
 
        else
 
        {
          if (MetaData.GlobalDir != NULL) {
            // Global fs
            strcpy(dumpdirname, MetaData.GlobalDir);
            strcat(dumpdirname, "/");
            strcat(dumpdirname, MetaData.ExtraDumpDir);
            strcat(dumpdirname, id);
            // Do mkdir on cpu 0 only
            strcpy(name, dumpdirname);
            strcat(name, "/");
            strcat(name, basename);
          } // if GlobalDir
 
          else
 
          {
            // No local or global specified
            strcpy(name, basename);
          } // else GlobalDir
 
        } // else LocalDir
 
      } // if DataDumpDir
 
      else
 
      {
        strcpy(name, basename);
      } // else DataDumpDir
 
      if (debug) fprintf(stdout, "Extra dump: %s\n", name);

    } // if ExtraDumpName
    /******************** CYCLE / DT BASED OUTPUTS ********************/
 
    if ( (cptr = strstr(basename, MetaData.DataDumpName)) ) {
 
      if (MetaData.DataDumpDir != NULL)
      {
        if (MetaData.LocalDir != NULL) {
          // Local fs
          strcpy(dumpdirname, MetaData.LocalDir);
          strcat(dumpdirname, "/");
          strcat(dumpdirname, MetaData.DataDumpDir);
          strcat(dumpdirname, id);
  
          // Create once per node...
#ifdef USE_NODE_LOCAL
            strcat(dumpdirname, "/mpi");
            strcat(dumpdirname, pid);
#endif /* USE_NODE_LOCAL */
 
          strcpy(name, dumpdirname);
          strcat(name, "/");
          strcat(name, basename);
        } // if LocalDir
 
        else
 
        {
          if (MetaData.GlobalDir != NULL) {
            // Global fs
            strcpy(dumpdirname, MetaData.GlobalDir);
            strcat(dumpdirname, "/");
            strcat(dumpdirname, MetaData.DataDumpDir);
            strcat(dumpdirname, id);
            // Do mkdir on cpu 0 only
            strcpy(name, dumpdirname);
            strcat(name, "/");
            strcat(name, basename);
          } // if GlobalDir
 
          else
 
          {
            // No local or global specified
            strcpy(name, basename);
          } // else GlobalDir
 
        } // else LocalDir
 
      } // if DataDumpDir
 
      else
 
      {
        strcpy(name, basename);
      } // else DataDumpDir
 
      if (debug) fprintf(stdout, "DATA dump: %s\n", name);
 
    } // if DataDumpName

    /******************** REDSHIFT BASED OUTPUTS ********************/
 
    if ( (cptr = strstr(basename, MetaData.RedshiftDumpName)) ) {
 
      if (MetaData.RedshiftDumpDir != NULL)
      {
        if (MetaData.LocalDir != NULL) {
          // Local fs
          strcpy(dumpdirname, MetaData.LocalDir);
          strcat(dumpdirname, "/");
          strcat(dumpdirname, MetaData.RedshiftDumpDir);
          strcat(dumpdirname, id);
  
          // Create once per node...
#ifdef USE_NODE_LOCAL
          strcat(dumpdirname, "/mpi");
          strcat(dumpdirname, pid);
#endif /* USE_NODE_LOCAL */
 
          strcpy(name, dumpdirname);
          strcat(name, "/");
          strcat(name, basename);
        } // if LocalDir
 
        else
 
        {
          if (MetaData.GlobalDir != NULL) {
            // Global fs
            strcpy(dumpdirname, MetaData.GlobalDir);
            strcat(dumpdirname, "/");
            strcat(dumpdirname, MetaData.RedshiftDumpDir);
            strcat(dumpdirname, id);
            // Do mkdir on cpu 0 only
            strcpy(name, dumpdirname);
            strcat(name, "/");
            strcat(name, basename);
          } // if GlobalDir
 
          else
 
          {
            // No local or global specified
            strcpy(name, basename);
          } // else GlobalDir
 
        } // else LocalDir
 
      } // if RedshiftDumpDir
 
      else
 
      {
        strcpy(name, basename);
      } // else RedshiftDumpDir

      if (debug)
	fprintf(stdout, "REDSHIFT dump: %s\n", name);
 
    } // if RedshiftDumpName

    /******************** RESTART BASED OUTPUTS ********************/

    if ( (cptr = strstr(basename, MetaData.RestartDumpName)) ) {
 
      if (MetaData.RestartDumpDir != NULL)
      {
        if (MetaData.LocalDir != NULL) {
          // Local fs
          strcpy(dumpdirname, MetaData.LocalDir);
          strcat(dumpdirname, "/");
          strcat(dumpdirname, MetaData.RestartDumpDir);
          strcat(dumpdirname, id);
  
          // Create once per node...
#ifdef USE_NODE_LOCAL
          strcat(dumpdirname, "/mpi");
          strcat(dumpdirname, pid);
#endif /* USE_NODE_LOCAL */
 
          strcpy(name, dumpdirname);
          strcat(name, "/");
          strcat(name, basename);
        } // if LocalDir
 
        else
 
        {
          if (MetaData.GlobalDir != NULL) {
            // Global fs
            strcpy(dumpdirname, MetaData.GlobalDir);
            strcat(dumpdirname, "/");
            strcat(dumpdirname, MetaData.RestartDumpDir);
            strcat(dumpdirname, id);
            // Do mkdir on cpu 0 only
            strcpy(name, dumpdirname);
            strcat(name, "/");
            strcat(name, basename);
          } // if GlobalDir
 
          else
 
          {
            // No local or global specified
            strcpy(name, basename);
          } // else GlobalDir
 
        } // else LocalDir
 
      } // if RedshiftDumpDir
 
      else
 
      {
        strcpy(name, basename);
      } // else RedshiftDumpDir

      if (debug)
	fprintf(stdout, "RESTART dump: %s\n", name);
 
    } // if RestartDumpName
 
    if (filenumber >= 0)
      strcat(name, id);
  }
 
  strcpy(LastFileNameWritten, name);
 
  strcpy(groupfilename, name);
  strcat(groupfilename, CPUSuffix);
  strcat(groupfilename, pid);
 
  if (debug)
    fprintf(stdout, "WriteAllData: writing group file %s\n", groupfilename);
 
//  Synchronization point for directory creation
 
  CommunicationBarrier();
#ifdef USE_MPI
  twrite0 = MPI_Wtime();
#endif
 
//  Get cwd
//  Generate command
//  Execute system call
 
    if ( local )
    {

      MPI_Arg mpi_size;
      MPI_Arg mpi_rank;

#ifdef USE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else
      mpi_rank = 0;
      mpi_size = 1;
#endif

      pe = mpi_rank;
      nn = mpi_size;

 
      for ( ii = 0; ii < nn; ii++ )
      {

        CommunicationBarrier();
        if( pe == ii )
        {
          if ( (cptr = strstr(basename, MetaData.ExtraDumpName)) ) {
            if (MetaData.ExtraDumpDir != NULL) {
#ifdef SYSCALL
              unixresult = SysMkdir("", dumpdirname);
              if (debug) fprintf(stdout, "Extra dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
              strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
              unixresult = system(unixcommand);
              if (debug) fprintf(stdout, "DATA dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
            }

          } // ENDIF extradump
 
          if ( (cptr = strstr(basename, MetaData.DataDumpName)) ) {
            if (MetaData.DataDumpDir != NULL) {
#ifdef SYSCALL
              unixresult = SysMkdir("", dumpdirname);
              if (debug) fprintf(stdout, "DATA dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
              strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
              unixresult = system(unixcommand);
              if (debug) fprintf(stdout, "DATA dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
            }
          } // ENDIF datadump
 
          if ( (cptr = strstr(basename, MetaData.RedshiftDumpName)) ) {
            if (MetaData.RedshiftDumpDir != NULL) {
#ifdef SYSCALL
              unixresult = SysMkdir("", dumpdirname);
              fprintf(stdout, "REDSHIFT dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
              strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
              unixresult = system(unixcommand);
              fprintf(stdout, "REDSHIFT dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
            }
          } // ENDIF redshift

          if ( (cptr = strstr(basename, MetaData.RestartDumpName)) ) {
            if (MetaData.RestartDumpDir != NULL) {
#ifdef SYSCALL
              unixresult = SysMkdir("", dumpdirname);
              fprintf(stdout, "RESTART dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
              strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
              unixresult = system(unixcommand);
              fprintf(stdout, "RESTART dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
            }
          } // ENDIF restart

        } // ENDIF pe == ii
      } // ENDFOR ii
    } // ENDIF local
 
    if ( global )
    {
      if ( MyProcessorNumber == ROOT_PROCESSOR )
      {
        if ( (cptr = strstr(basename, MetaData.ExtraDumpName)) ) {
          if (MetaData.ExtraDumpDir != NULL) {
#ifdef SYSCALL
            unixresult = SysMkdir("", dumpdirname);
            if (debug) fprintf(stdout, "Extra dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
            strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
            unixresult = system(unixcommand);
            if (debug) fprintf(stdout, "Extra dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
          }
        } // ENDIF datadump
 
        if ( (cptr = strstr(basename, MetaData.DataDumpName)) ) {
          if (MetaData.DataDumpDir != NULL) {
#ifdef SYSCALL
            unixresult = SysMkdir("", dumpdirname);
            if (debug) fprintf(stdout, "DATA dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
            strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
            unixresult = system(unixcommand);
            if (debug) fprintf(stdout, "DATA dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
          }
        } // ENDIF datadump
 
        if ( (cptr = strstr(basename, MetaData.RedshiftDumpName)) ) {
          if (MetaData.RedshiftDumpDir != NULL) {
#ifdef SYSCALL
            unixresult = SysMkdir("", dumpdirname);
            fprintf(stdout, "REDSHIFT dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
            strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
            unixresult = system(unixcommand);
            fprintf(stdout, "REDSHIFT dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
          }
        } // ENDIF redshift

        if ( (cptr = strstr(basename, MetaData.RestartDumpName)) ) {
          if (MetaData.RestartDumpDir != NULL) {
#ifdef SYSCALL
            unixresult = SysMkdir("", dumpdirname);
            fprintf(stdout, "RESTART dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
            strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
            unixresult = system(unixcommand);
            fprintf(stdout, "RESTART dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
          }
        } // ENDIF restart
 
      }
    }


  
//  fprintf(stdout, "Sync point ok\n");
 
  CommunicationBarrier();
 
//  Start I/O timing
 
#ifdef USE_HDF5_OUTPUT_BUFFERING

  memory_increment = 1024*1024;
  dump_flag = 1;

  file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
    if( file_acc_template == h5_error ){my_exit(EXIT_FAILURE);}

  h5_status = H5Pset_fapl_core(file_acc_template, memory_increment, dump_flag);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

  file_id = H5Fcreate(groupfilename, H5F_ACC_TRUNC, H5P_DEFAULT, file_acc_template);
    if( file_id == h5_error ){my_exit(EXIT_FAILURE);}

#else

  file_id = H5Fcreate(groupfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //  h5_status = H5Fclose(file_id);
  //    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#endif

  // WS: Output forcing spectrum
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if (DrivenFlowProfile) {
      strcpy(forcingname, name);
      strcat(forcingname, ForcingSuffix);
      if (debug)
          printf("Group_WriteAllData: writing file %s.\n", forcingname);
      if (Forcing.WriteSpectrum(forcingname) == FAIL) {
          fprintf(stderr, "Error in WriteSpectrum.\n");
          return FAIL;
      }
    }
  } 

 
  // Set MetaData.BoundaryConditionName
 
  if (MetaData.BoundaryConditionName != NULL)
    delete [] MetaData.BoundaryConditionName;
  MetaData.BoundaryConditionName = new char[MAX_LINE_LENGTH];
  strcpy(MetaData.BoundaryConditionName, name);
  strcat(MetaData.BoundaryConditionName, BCSuffix);

  /* We set our global variable CheckpointRestart to TRUE here, so that it gets
     output in the parameter file. */

  CheckpointRestart = CheckpointDump;
 
#ifdef TRANSFER
  if (ImplicitProblem && MyProcessorNumber == ROOT_PROCESSOR) {
    // Output ImplicitSolver module parameter file

    //    Reset MetaData.RadHydroParameterFname
    if (MetaData.RadHydroParameterFname != NULL)
      delete MetaData.RadHydroParameterFname;
    MetaData.RadHydroParameterFname = new char[MAX_LINE_LENGTH];
    strcpy(MetaData.RadHydroParameterFname, name);
    strcat(MetaData.RadHydroParameterFname, RTSuffix);
    
    // Open RT module parameter file
    if ((fptr = fopen(MetaData.RadHydroParameterFname, "w")) == NULL) {
      fprintf(stderr, "Error opening RT module parameter file: %s\n",
	      MetaData.RadHydroParameterFname);
      return FAIL;
    }
    
    // Write RT module parameters to file
    if (ImplicitSolver->WriteParameters(fptr) == FAIL) {
      fprintf(stderr, "Error in ImplicitSolver::WriteParameters\n");
      return FAIL;
    }
    fclose(fptr);
  }
#endif		 

  // Output TopGrid data
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((fptr = fopen(name, "w")) == NULL) 
      ENZO_VFAIL("Error opening output file %s\n", name)
    if (CheckpointDump == TRUE) {
      fprintf(fptr, "# WARNING! This is a checkpoint dump! Lots of data!\n");
    }
    else if (WriteTime >= 0) {
      fprintf(fptr, "# WARNING! Interpolated output: level = %"ISYM"\n",
	      MetaData.OutputFirstTimeAtLevel-1);
    }
    if (WriteParameterFile(fptr, MetaData, name) == FAIL)
      ENZO_FAIL("Error in WriteParameterFile");
    fclose(fptr);
  
  }

 
  // Output Boundary condition info
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((fptr = fopen(MetaData.BoundaryConditionName, "w")) == NULL) 
      ENZO_VFAIL("Error opening boundary condition file: %s\n",
	      MetaData.BoundaryConditionName)
    strcat(MetaData.BoundaryConditionName, hdfsuffix);
    if (Exterior->WriteExternalBoundary(fptr, MetaData.BoundaryConditionName)
	== FAIL)
      ENZO_FAIL("Error in WriteExternalBoundary");
    fclose(fptr);
 
  }
 
  // Create hierarchy name and grid base name
 
  strcpy(hierarchyname, name);
  strcat(hierarchyname, HierarchySuffix);
 
  strcpy(gridbasename, name);
//  strcat(gridbasename, GridSuffix);

  strcpy(taskmapname, name);
  strcat(taskmapname, TaskMapSuffix);
  strcat(taskmapname, pid);

  strcpy(memorymapname, name);
  strcat(memorymapname, MemoryMapSuffix);

  strcpy(configurename, name);
  strcat(configurename, ConfigureSuffix);

 
  /* If we're writing interpolated dark matter fields, create them now. */

  CreateSmoothedDarkMatterFields(MetaData, TopGrid);
 
  /* Combine the top level grids into a single grid for output
     (TempTopGrid is the top of an entirely new hierarchy). */
 
  int level;
  HierarchyEntry *TempTopGrid;
  CommunicationCombineGrids(TopGrid, &TempTopGrid, WriteTime, CheckpointDump);
 
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;

#ifndef FAST_SIB
    if (VelAnyl==1 || BAnyl==1) {
      AddLevel(LevelArray, TempTopGrid, 0);
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	HierarchyEntry **Grids;
	int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
	if (LevelArray[level] != NULL)
	  SetBoundaryConditions(Grids, NumberOfGrids, level, &MetaData, 
				Exterior, LevelArray[level]);
      }
    }
#endif
  
  // Output Data Hierarchy

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    if( HierarchyFileOutputFormat % 2 == 0 ) {
      AddLevel(LevelArray, TempTopGrid, 0);
      WriteHDF5HierarchyFile(name, TempTopGrid, MetaData, LevelArray);      
    }
    
    if (HierarchyFileOutputFormat > 0)
      if ((fptr = fopen(hierarchyname, "w")) == NULL) 
	ENZO_VFAIL("Error opening hierarchy file %s\n", hierarchyname);
  }

  /* Clean-up LevelArray */

  if (LevelArray[0] != NULL) {
    LevelHierarchyEntry *Temp;
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      while (LevelArray[level] != NULL) {
	Temp = LevelArray[level]->NextGridThisLevel;
	delete LevelArray[level];
	LevelArray[level] = Temp;
      } // ENDWHILE
  } // ENDIF

  if (Group_WriteDataHierarchy(fptr, MetaData, TempTopGrid,
            gridbasename, GridID, WriteTime, file_id, CheckpointDump) == FAIL)
    ENZO_FAIL("Error in Group_WriteDataHierarchy");

    hid_t metadata_group = H5Gcreate(file_id, "Metadata", 0);
    if(metadata_group == h5_error)ENZO_FAIL("Error writing metadata!");
    writeArrayAttribute(metadata_group, HDF5_INT, MAX_DEPTH_OF_HIERARCHY,
                        "LevelCycleCount", LevelCycleCount);
  if(CheckpointDump == TRUE){
    // Write our supplemental (global) data
    FLOAT dtThisLevelCopy[MAX_DEPTH_OF_HIERARCHY];
    FLOAT dtThisLevelSoFarCopy[MAX_DEPTH_OF_HIERARCHY];
    for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
        dtThisLevelCopy[level] = dtThisLevel[level];
        dtThisLevelSoFarCopy[level] = dtThisLevelSoFar[level];
    }
    writeArrayAttribute(metadata_group, HDF5_PREC, MAX_DEPTH_OF_HIERARCHY,
                        "dtThisLevel", dtThisLevelCopy);
    writeArrayAttribute(metadata_group, HDF5_PREC, MAX_DEPTH_OF_HIERARCHY,
                        "dtThisLevelSoFar", dtThisLevelSoFarCopy);
    writeScalarAttribute(metadata_group, HDF5_PREC, "Time", &MetaData.Time);
  }
    H5Gclose(metadata_group);

  // At this point all the grid data has been written

  h5_status = H5Fclose(file_id);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#ifdef USE_HDF5_OUTPUT_BUFFERING

  h5_status = H5Pclose(file_acc_template);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#endif


  if (MyProcessorNumber == ROOT_PROCESSOR)
    if ((mptr = fopen(memorymapname, "w")) == NULL) 
      ENZO_VFAIL("Error opening memory map file %s\n", memorymapname)

  if (WriteMemoryMap(mptr, TempTopGrid, gridbasename, GridKD, WriteTime) == FAIL)
    ENZO_FAIL("Error in WriteMemoryMap");

  // Output configure

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((optr = fopen(configurename, "w")) == NULL) {
      fprintf(stdout, "Error opening configure file %s\n", configurename);
      fprintf(stdout, "Not crucial but worrysome. Will continue.\n" );
      //      ENZO_VFAIL("Error opening configure file %s\n", configurename)
    }

    WriteConfigure(optr);

    fclose(optr);
  }

  // Output task map

#ifdef TASKMAP
  if ((tptr = fopen(taskmapname, "w")) == NULL)
    ENZO_VFAIL("Error opening task map file %s\n", taskmapname)

  if (WriteTaskMap(tptr, TempTopGrid, gridbasename, GridLD, WriteTime) == FAIL)
    ENZO_FAIL("Error in WriteTaskMap");
#endif
 
  int TGdims[3];
 
  TGdims[0] = MetaData.TopGridDims[0];
  TGdims[1] = MetaData.TopGridDims[1];
  TGdims[2] = MetaData.TopGridDims[2];
 
  //  fprintf(stdout, "TGdims  %"ISYM"  %"ISYM"  %"ISYM"\n", TGdims[0], TGdims[1], TGdims[2]);
 
  if (CubeDumpEnabled == 1)
    if (WriteDataCubes(TempTopGrid, TGdims, name, GridJD, WriteTime) == FAIL)
      ENZO_FAIL("Error in WriteDataCubes");

 
  // Clean up combined top level grid, and first two levels of hierarchy
 
  if (TempTopGrid != TopGrid) {
    if (TempTopGrid->NextGridNextLevel != NULL)
      DeleteGridHierarchy(TempTopGrid->NextGridNextLevel);
    delete TempTopGrid->GridData;
    delete TempTopGrid;
  }
 
  // Output star particle data 
 
  if (WriteStarParticleData(fptr, MetaData) == FAIL)
    ENZO_FAIL("Error in WriteStarParticleData");
 
  // Output MBH particle data
  
  if (MBHParticleIO == TRUE && MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *MBHfptr;

    if ((MBHfptr = fopen(MBHParticleIOFilename, "a")) == NULL) {
      ENZO_VFAIL("Error opening file %s\n", MBHParticleIOFilename)
    }
    
    // printing order: time, regular star count, MBH id, MBH mass, MBH angular momentum
    for (int i = 0; i < G_TotalNumberOfStars; i++) { 
      fprintf(MBHfptr, " %"FSYM"  %"ISYM"  %"ISYM"  %lf  %"FSYM"  %"FSYM"  %"FSYM"  %lf\n", 
	      MetaData.Time, NumberOfStarParticles, (int)(MBHParticleIOTemp[i][0]), 
	      MBHParticleIOTemp[i][1], (float)(MBHParticleIOTemp[i][2]), 
	      (float)(MBHParticleIOTemp[i][3]), (float)(MBHParticleIOTemp[i][4]),
	      MBHParticleIOTemp[i][5]);
    }
    
    fclose(MBHfptr);
    
  }

  // Create radiation name and write radiation data
 
  if (((RadiationFieldType >= 10 && RadiationFieldType <= 11) ||
       RadiationData.RadiationShield == TRUE) &&
      MyProcessorNumber == ROOT_PROCESSOR) {
 
    FILE *Radfptr;
 
    strcpy(radiationname, name);
    strcat(radiationname, RadiationSuffix);
 
    if ((Radfptr = fopen(radiationname, "w")) == NULL) {
      ENZO_VFAIL("Error opening radiation file %s\n", radiationname)
    }
    if (WriteRadiationData(Radfptr) == FAIL)
      ENZO_FAIL("Error in WriteRadiationData");
 
    fclose(Radfptr);
 
  }
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if (HierarchyFileOutputFormat  > 0)
      fclose(fptr);
    fclose(mptr);
  }

#ifdef TASKMAP
  fclose(tptr);
#endif

  // Replace the time in metadata with the saved value (above)
 
  MetaData.Time = SavedTime;

  /* This should always be false outside of writing and the restart process */
    
  CheckpointRestart = FALSE;

  CommunicationBarrier();
// if (debug)
    //  fprintf(stdout, "WriteAllData: finished writing data\n");
 
//  Stop I/O timing
 
//  Synchronization point for SRB
 
  CommunicationBarrier();
 
  ContinueExecution();
 
  CommunicationBarrier();
#ifdef USE_MPI
  twrite1 = MPI_Wtime();
#endif

  if ( MyProcessorNumber == ROOT_PROCESSOR ){
    sptr = fopen("OutputLog", "a");
    fprintf(sptr, "DATASET WRITTEN %s %8"ISYM" %18.16"GSYM" %18.8"FSYM" %18.8"FSYM"\n", 
	    name, MetaData.CycleNumber, MetaData.Time, twrite0, (twrite1-twrite0));
    fclose(sptr);
  }

  TIMER_STOP("Group_WriteAllData"); 
  return SUCCESS;
}

/* 
void DeleteGridHierarchy(HierarchyEntry *GridEntry)
{
  if (GridEntry->NextGridThisLevel != NULL)

     DeleteGridHierarchy(GridEntry->NextGridThisLevel);
 
  delete GridEntry;
 
  return;
}
*/
