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
/ modified2:   Robert HArkness
/ date:        April 2008
/
/  PURPOSE:
/
************************************************************************/
 
// This function writes out the data hierarchy (TopGrid), the External
//   Boundary (Exterior), the TopGridData, and the global_data.
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
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
#include "CosmologyParameters.h"
#include "CommunicationUtilities.h"
void my_exit(int status);
 
// function prototypes
 
int SysMkdir(char *startdir, char *directory);
 
int WriteDataCubes(HierarchyEntry *TopGrid, int TDdims[], char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteDataHierarchy(FILE *fptr, TopGridData &MetaData, HierarchyEntry *TopGrid,
		       char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteMemoryMap(FILE *fptr, HierarchyEntry *TopGrid,
		   char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteConfigure(FILE *optr);
int WriteTaskMap(FILE *fptr, HierarchyEntry *TopGrid,
		 char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteParameterFile(FILE *fptr, TopGridData &MetaData);
int WriteStarParticleData(FILE *fptr);
int WriteRadiationData(FILE *fptr);
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int CommunicationCombineGrids(HierarchyEntry *OldHierarchy,
			      HierarchyEntry **NewHierarchyPointer,
			      FLOAT WriteTime, int RestartDump = FALSE);
void DeleteGridHierarchy(HierarchyEntry *GridEntry);
void ContinueExecution(void);
int CreateSmoothedDarkMatterFields(TopGridData &MetaData, HierarchyEntry *TopGrid);
 
 
char BCSuffix[]        = ".boundary";
char GridSuffix[]      = ".grid";
char HierarchySuffix[] = ".hierarchy";
char hdfsuffix[]       = ".hdf";
char RadiationSuffix[] = ".radiation";
char TaskMapSuffix[]   = ".taskmap";
char MemoryMapSuffix[] = ".memorymap";
char ConfigureSuffix[] = ".configure";
 
extern char LastFileNameWritten[MAX_LINE_LENGTH];
 
 
 
 
int WriteAllData(char *basename, int filenumber,
		 HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior, FLOAT WriteTime = -1)
{
 
 
  char id[MAX_CYCLE_TAG_SIZE], *cptr, name[MAX_LINE_LENGTH];
  char dumpdirname[MAX_LINE_LENGTH];
  char dumpdirroot[MAX_LINE_LENGTH];
  char unixcommand[MAX_LINE_LENGTH];
  char gridbasename[MAX_LINE_LENGTH];
  char hierarchyname[MAX_LINE_LENGTH];
  char radiationname[MAX_LINE_LENGTH];
  char taskmapname[MAX_LINE_LENGTH];
  char memorymapname[MAX_LINE_LENGTH];
  char configurename[MAX_LINE_LENGTH];
 
  int unixresult;
  int status;
  int local, global;
  int file_status;
  int ii, pe, nn;
 
#ifdef USE_MPI
  double io_start, io_stop;
  double dc_start, dc_stop;
  double ttenter, ttexit;
  double iot1a, iot1b, iot2a, iot2b, iot3a, iot3b, iot4a, iot4b;
  char io_logfile[MAX_NAME_LENGTH];
  FILE *xptr;
#endif /* USE_MPI */
  char pid[MAX_TASK_TAG_SIZE];
 
  FILE *fptr;
  FILE *sptr;
  FILE *tptr;
  FILE *mptr;
  FILE *optr;
 
  int GridID = 1;
  int GridJD = 1;
  int GridKD = 1;
  int GridLD = 1;
 
#ifdef USE_MPI
  ttenter = MPI_Wtime();
#endif /* USE_MPI */

  /* If this is an interpolated time step, then temporary replace  the time
     in MetaData.  Note:  Modified 6 Feb 2006 to fix interpolated  data outputs. */

  FLOAT SavedTime = MetaData.Time;
  MetaData.Time = (WriteTime < 0) ? MetaData.Time : WriteTime;

  /* If we're writing interpolated dark matter fields, create them now. */

  CreateSmoothedDarkMatterFields(MetaData, TopGrid);

  // Global or local filesystem?
 
  local = 0;
  global = 0;
 
  strcpy(dumpdirname, " ");
 
  if (MetaData.LocalDir != NULL)
  {
     local = 1;
     strcpy(dumpdirroot, MetaData.LocalDir);
     // fprintf(stderr, "XXXX local dir: %s\n", MetaData.LocalDir);
     // Create on node - locking?
  }
 
  if (MetaData.GlobalDir != NULL)
  {
     global = 1;
     strcpy(dumpdirroot, MetaData.GlobalDir);
     // fprintf(stderr, "XXXX global dir: %s\n", MetaData.GlobalDir);
     // Create on task 0 only
  }
 
  if (( local == 1) && (global == 1))
    fprintf(stderr, "Local AND Global !!\n");
 
  // Create main name
 
  if (ComovingCoordinates && (cptr = strstr(name, "RRRR"))) {
    FLOAT a, dadt;
    CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
    sprintf(cptr, "%"CYCLE_TAG_FORMAT""ISYM, nint(100*((1 + InitialRedshift)/a - 1)));
  } else {
 
    sprintf(id, "%"CYCLE_TAG_FORMAT""ISYM, filenumber);
 
#ifdef USE_MPI
    sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
#endif /* USE_MPI */
 
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
 
      if (debug) fprintf(stderr, "DATA dump: %s\n", name);
 
    } // if DataDumpName
 
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
 
      fprintf(stderr, "REDSHIFT dump: %s\n", name);
 
    } // if RedshiftDumpName
 
    if (filenumber >= 0)
      strcat(name, id);
  }
 
  strcpy(LastFileNameWritten, name);
 
//  if (debug)
    fprintf(stderr, "WriteAllData: writing file %s\n", name);
 
//  Synchronization point for directory creation
 
#ifdef USE_MPI
  iot1a = MPI_Wtime();
  CommunicationBarrier();
  dc_start = MPI_Wtime();
  iot1b = MPI_Wtime();
#endif /* USE_MPI */
 
//  Get cwd
//  Generate command
//  Execute system call
 
    if ( local )
    {

      MPI_Arg mpi_rank;
      MPI_Arg mpi_size;

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
 
          if ( (cptr = strstr(basename, MetaData.DataDumpName)) ) {
            if (MetaData.DataDumpDir != NULL) {
#ifdef SYSCALL
              unixresult = SysMkdir("", dumpdirname);
              if (debug) fprintf(stderr, "DATA dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
              strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
              unixresult = system(unixcommand);
              if (debug) fprintf(stderr, "DATA dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
            }
          }
 
          if ( (cptr = strstr(basename, MetaData.RedshiftDumpName)) ) {
            if (MetaData.RedshiftDumpDir != NULL) {
#ifdef SYSCALL
              unixresult = SysMkdir("", dumpdirname);
              fprintf(stderr, "REDSHIFT dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
              strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
              unixresult = system(unixcommand);
              fprintf(stderr, "REDSHIFT dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
            }
          }
        }
      }
    }
 
    if ( global )
    {
      if ( MyProcessorNumber == ROOT_PROCESSOR )
      {
 
        if ( (cptr = strstr(basename, MetaData.DataDumpName)) ) {
          if (MetaData.DataDumpDir != NULL) {
#ifdef SYSCALL
            unixresult = SysMkdir("", dumpdirname);
            if (debug) fprintf(stderr, "DATA dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
            strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
            unixresult = system(unixcommand);
            if (debug) fprintf(stderr, "DATA dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
          }
        }
 
        if ( (cptr = strstr(basename, MetaData.RedshiftDumpName)) ) {
          if (MetaData.RedshiftDumpDir != NULL) {
#ifdef SYSCALL
            unixresult = SysMkdir("", dumpdirname);
            fprintf(stderr, "REDSHIFT dump: dumpdirname=(%s) == unixresult=%"ISYM"\n", dumpdirname, unixresult);
#else
            strcat(strcpy(unixcommand, "mkdir -p "), dumpdirname);
            unixresult = system(unixcommand);
            fprintf(stderr, "REDSHIFT dump: %s == %"ISYM"\n", unixcommand, unixresult);
#endif
          }
        }
 
      }
    }
 
//  fprintf(stderr, "Sync point ok\n");
 
#ifdef USE_MPI
  iot2a = MPI_Wtime();
  CommunicationBarrier();
  dc_stop = MPI_Wtime();
  iot2b = MPI_Wtime();
#endif /* USE_MPI */
 
//  Start I/O timing
 
#ifdef USE_MPI
  io_start = MPI_Wtime();
#endif /* USE_MPI */
 
  // Set MetaData.BoundaryConditionName
 
  if (MetaData.BoundaryConditionName != NULL)
    delete MetaData.BoundaryConditionName;
  MetaData.BoundaryConditionName = new char[MAX_LINE_LENGTH];
  strcpy(MetaData.BoundaryConditionName, name);
  strcat(MetaData.BoundaryConditionName, BCSuffix);
 
  // Output TopGrid data
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((fptr = fopen(name, "w")) == NULL) {
      fprintf(stderr, "Error opening output file %s\n", name);
      ENZO_FAIL("");
    }
    if (WriteTime >= 0)
      fprintf(fptr, "# WARNING! Interpolated output: level = %"ISYM"\n",
	      MetaData.OutputFirstTimeAtLevel-1);
    if (WriteParameterFile(fptr, MetaData) == FAIL) {
      fprintf(stderr, "Error in WriteParameterFile\n");
      ENZO_FAIL("");
    }
    fclose(fptr);
  
  }
 
  // Output Boundary condition info
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((fptr = fopen(MetaData.BoundaryConditionName, "w")) == NULL) {
      fprintf(stderr, "Error opening boundary condition file: %s\n",
	      MetaData.BoundaryConditionName);
      ENZO_FAIL("");
    }
    strcat(MetaData.BoundaryConditionName, hdfsuffix);
    if (Exterior->WriteExternalBoundary(fptr, MetaData.BoundaryConditionName)
	== FAIL) {
      fprintf(stderr, "Error in WriteExternalBoundary\n");
      ENZO_FAIL("");
    }
    fclose(fptr);
  
  }
 
  // Create hierarchy name and grid base name
 
  strcpy(hierarchyname, name);
  strcat(hierarchyname, HierarchySuffix);
 
  strcpy(gridbasename, name);
  strcat(gridbasename, GridSuffix);

  strcpy(taskmapname, name);
  strcat(taskmapname, TaskMapSuffix);
  strcat(taskmapname, pid);

  strcpy(memorymapname, name);
  strcat(memorymapname, MemoryMapSuffix);

  strcpy(configurename, name);
  strcat(configurename, ConfigureSuffix);
 
  /* Combine the top level grids into a single grid for output
     (TempTopGrid is the top of an entirely new hierarchy). */
 
  HierarchyEntry *TempTopGrid;
  CommunicationCombineGrids(TopGrid, &TempTopGrid, WriteTime);
 
  // Output Data Hierarchy
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    if ((fptr = fopen(hierarchyname, "w")) == NULL) {
      fprintf(stderr, "Error opening hierarchy file %s\n", hierarchyname);
      ENZO_FAIL("");
    }
 
  if (WriteDataHierarchy(fptr, MetaData, TempTopGrid, gridbasename, GridID, WriteTime) == FAIL) {
    fprintf(stderr, "Error in WriteDataHierarchy\n");
    ENZO_FAIL("");
  }

  // Output StarParticle data (actually just number of stars)
 
  if (WriteStarParticleData(fptr) == FAIL) {
    fprintf(stderr, "Error in WriteStarParticleData\n");
    ENZO_FAIL("");
  }
 
  // Output memory map

  if (MyProcessorNumber == ROOT_PROCESSOR)
    if ((mptr = fopen(memorymapname, "w")) == NULL) {
      fprintf(stderr, "Error opening memory map file %s\n", memorymapname);
      ENZO_FAIL("");
    }

  if (WriteMemoryMap(mptr, TempTopGrid, gridbasename, GridKD, WriteTime) == FAIL) {
    fprintf(stderr, "Error in WriteMemoryMap\n");
    ENZO_FAIL("");
  }

  // Output configure

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((optr = fopen(configurename, "w")) == NULL) {
      fprintf(stderr, "Error opening configure file %s\n", configurename);
      ENZO_FAIL("");
    }

    WriteConfigure(optr);

    fclose(optr);
  }

  // Output task map

  if ((tptr = fopen(taskmapname, "w")) == NULL) {
    fprintf(stderr, "Error opening task map file %s\n", taskmapname);
    ENZO_FAIL("");
  }

  if (WriteTaskMap(tptr, TempTopGrid, gridbasename, GridLD, WriteTime) == FAIL) {
    fprintf(stderr, "Error in WriteTaskMap\n");
    ENZO_FAIL("");
  }
 
  int TGdims[3];
 
  TGdims[0] = MetaData.TopGridDims[0];
  TGdims[1] = MetaData.TopGridDims[1];
  TGdims[2] = MetaData.TopGridDims[2];
 
  //  fprintf(stderr, "TGdims  %"ISYM"  %"ISYM"  %"ISYM"\n", TGdims[0], TGdims[1], TGdims[2]);
 
  if (CubeDumpEnabled == 1) {
    if (WriteDataCubes(TempTopGrid, TGdims, name, GridJD, WriteTime) == FAIL) {
      fprintf(stderr, "Error in WriteDataCubes\n");
      ENZO_FAIL("");
    }
  }
 
  // Clean up combined top level grid, and first two levels of hierarchy
 
  if (TempTopGrid != TopGrid) {
    if (TempTopGrid->NextGridNextLevel != NULL)
      DeleteGridHierarchy(TempTopGrid->NextGridNextLevel);
    delete TempTopGrid->GridData;
    delete TempTopGrid;
  }
 
  // Create radiation name and write radiation data
 
  if (RadiationFieldType >= 10 && RadiationFieldType <= 11 &&
      MyProcessorNumber == ROOT_PROCESSOR) {
 
    FILE *Radfptr;
 
    strcpy(radiationname, name);
    strcat(radiationname, RadiationSuffix);
 
    if ((Radfptr = fopen(radiationname, "w")) == NULL) {
      fprintf(stderr, "Error opening radiation file %s\n", radiationname);
      ENZO_FAIL("");
    }
    if (WriteRadiationData(Radfptr) == FAIL) {
      fprintf(stderr, "Error in WriteRadiationData\n");
      ENZO_FAIL("");
    }
 
    fclose(Radfptr);
 
  }
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fclose(fptr);
    fclose(mptr);
  }

  fclose(tptr);
 
  // Replace the time in metadata with the saved value (above)
 
  MetaData.Time = SavedTime;

  CommunicationBarrier();
// if (debug)
    //  fprintf(stderr, "WriteAllData: finished writing data\n");
 
//  Stop I/O timing
 
#ifdef USE_MPI
  io_stop = MPI_Wtime();
#endif /* USE_MPI */
 
#ifdef USE_MPI
  iot3a = MPI_Wtime();
  CommunicationBarrier();
  iot3b = MPI_Wtime();
#endif /* USE_MPI */
 
  ContinueExecution();
 
#ifdef USE_MPI
  iot4a = MPI_Wtime();
  CommunicationBarrier();
  iot4b = MPI_Wtime();
#endif /* USE_MPI */

  if ( MyProcessorNumber == ROOT_PROCESSOR ){
    sptr = fopen("OutputLog", "a");
    fprintf(sptr, "DATASET WRITTEN %s \n", name);
    fclose(sptr);
  }
 
#ifdef USE_MPI
  ttexit = MPI_Wtime();
#endif /* USE_MPI */
 
#ifdef USE_MPI
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
  strcpy(io_logfile, "IO_perf.");
  strcat(io_logfile, pid);
  xptr = fopen(io_logfile, "a");
  fprintf(xptr, "IO %12.4e  %s\n", (io_stop-io_start), name);
  fprintf(xptr, "DC %12.4e  %s\n", (dc_stop-dc_start), name);
  fprintf(xptr, "XX %12.4e  %12.4e  %12.4e  %12.4e\n",
                (iot1b-iot1a), (iot2b-iot2a), (iot3b-iot3a), (iot4b-iot4a));
  fprintf(xptr, "TT %12.4e\n", (ttexit-ttenter));
  fclose(xptr);
#endif /* USE_MPI */

  //  fprintf(stderr,"Safe exit from WriteAllData\n");
 
  return SUCCESS;
}
 
void DeleteGridHierarchy(HierarchyEntry *GridEntry)
{
  if (GridEntry->NextGridThisLevel != NULL)
     DeleteGridHierarchy(GridEntry->NextGridThisLevel);
 
  delete GridEntry;
 

 
  return;
}
