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
/  modified2:  Robert HArkness
/  date:       April 2008
/  modified3:  Michael Kuhlen, October 2010, HDF5 hierarchy
/
/  PURPOSE:
/
************************************************************************/
 
// This function writes out the data hierarchy (TopGrid), the External
//   Boundary (Exterior), the TopGridData, and the global_data.
#include "preincludes.h"

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
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif
void my_exit(int status);
 
// function prototypes
 
int SysMkdir(char *startdir, char *directory);

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int WriteDataCubes(HierarchyEntry *TopGrid, int TDdims[], char *gridbasename, int &GridID, FLOAT WriteTime);
int WriteDataHierarchy(FILE *fptr, TopGridData &MetaData, HierarchyEntry *TopGrid,
		       char *gridbasename, int &GridID, FLOAT WriteTime);
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
			      FLOAT WriteTime, int RestartDump = FALSE);
void DeleteGridHierarchy(HierarchyEntry *GridEntry);
void ContinueExecution(void);
int CreateSmoothedDarkMatterFields(TopGridData &MetaData, HierarchyEntry *TopGrid);
 
 
#ifdef TRANSFER
char RTSuffix[]        = ".rtmodule";
#endif
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
		 ExternalBoundary *Exterior, 
#ifdef TRANSFER
		 ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1)
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
 
  /* If this is an interpolated time step, then temporary replace  the time
     in MetaData.  Note:  Modified 6 Feb 2006 to fix interpolated  data outputs. */

  FLOAT SavedTime = MetaData.Time;
  MetaData.Time = (WriteTime < 0) ? MetaData.Time : WriteTime;

  /* If we're writing interpolated dark matter fields, create them now. */

  CreateSmoothedDarkMatterFields(MetaData, TopGrid);

  /* Combine the top level grids into a single grid for output
     (TempTopGrid is the top of an entirely new hierarchy). */
 
  HierarchyEntry *TempTopGrid;
  CommunicationCombineGrids(TopGrid, &TempTopGrid, WriteTime);

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
 
  CommunicationBarrier();
 
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
 
  CommunicationBarrier();
 
  // Set MetaData.BoundaryConditionName
 
  if (MetaData.BoundaryConditionName != NULL)
    delete MetaData.BoundaryConditionName;
  MetaData.BoundaryConditionName = new char[MAX_LINE_LENGTH];
  strcpy(MetaData.BoundaryConditionName, name);
  strcat(MetaData.BoundaryConditionName, BCSuffix);

#ifdef TRANSFER
  if (ImplicitProblem) {
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
    if ((fptr = fopen(name, "w")) == NULL) {
      ENZO_VFAIL("Error opening output file %s\n", name)
    }
    if (WriteTime >= 0)
      fprintf(fptr, "# WARNING! Interpolated output: level = %"ISYM"\n",
	      MetaData.OutputFirstTimeAtLevel-1);
    if (WriteParameterFile(fptr, MetaData, name) == FAIL) {
      ENZO_FAIL("Error in WriteParameterFile\n");
    }
    fclose(fptr);
  
  }
 
  // Output Boundary condition info
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((fptr = fopen(MetaData.BoundaryConditionName, "w")) == NULL) {
      ENZO_VFAIL("Error opening boundary condition file: %s\n",
	      MetaData.BoundaryConditionName)
    }
    strcat(MetaData.BoundaryConditionName, hdfsuffix);
    if (Exterior->WriteExternalBoundary(fptr, MetaData.BoundaryConditionName)
	== FAIL) {
      ENZO_FAIL("Error in WriteExternalBoundary\n");
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
 
 
  // Output Data Hierarchy

  int level;
  LevelHierarchyEntry *Temp;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    if (HierarchyFileOutputFormat % 2 == 0) {
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
	LevelArray[level] = NULL;
      AddLevel(LevelArray, TempTopGrid, 0);
      WriteHDF5HierarchyFile(name, TempTopGrid, MetaData, LevelArray);

      // Delete LevelArray linked list
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
	while (LevelArray[level] != NULL) {
	  Temp = LevelArray[level]->NextGridThisLevel;
	  delete LevelArray[level];
	  LevelArray[level] = Temp;
	} // ENDWHILE
    } // ENDIF HierarchyFileOutputFormat

    if (HierarchyFileOutputFormat > 0)
      if ((fptr = fopen(hierarchyname, "w")) == NULL)
	ENZO_VFAIL("Error opening hierarchy file %s\n", hierarchyname);
  }
 
  if (WriteDataHierarchy(fptr, MetaData, TempTopGrid, gridbasename, GridID, WriteTime) == FAIL) {
    ENZO_FAIL("Error in WriteDataHierarchy\n");
  }

  // Output StarParticle data (actually just number of stars)
 
  if (WriteStarParticleData(fptr, MetaData) == FAIL) {
    ENZO_FAIL("Error in WriteStarParticleData\n");
  }

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
   
  // Output memory map

  if (MyProcessorNumber == ROOT_PROCESSOR)
    if ((mptr = fopen(memorymapname, "w")) == NULL) {
      ENZO_VFAIL("Error opening memory map file %s\n", memorymapname)
    }

  if (WriteMemoryMap(mptr, TempTopGrid, gridbasename, GridKD, WriteTime) == FAIL) {
    ENZO_FAIL("Error in WriteMemoryMap\n");
  }

  // Output configure

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if ((optr = fopen(configurename, "w")) == NULL) {
      ENZO_VFAIL("Error opening configure file %s\n", configurename)
    }

    WriteConfigure(optr);

    fclose(optr);
  }

  // Output task map

#ifdef TASKMAP
  if ((tptr = fopen(taskmapname, "w")) == NULL) {
    ENZO_VFAIL("Error opening task map file %s\n", taskmapname)
  }

  if (WriteTaskMap(tptr, TempTopGrid, gridbasename, GridLD, WriteTime) == FAIL) {
    ENZO_FAIL("Error in WriteTaskMap\n");
  }
#endif
 
  int TGdims[3];
 
  TGdims[0] = MetaData.TopGridDims[0];
  TGdims[1] = MetaData.TopGridDims[1];
  TGdims[2] = MetaData.TopGridDims[2];
 
  //  fprintf(stderr, "TGdims  %"ISYM"  %"ISYM"  %"ISYM"\n", TGdims[0], TGdims[1], TGdims[2]);
 
  if (CubeDumpEnabled == 1) {
    if (WriteDataCubes(TempTopGrid, TGdims, name, GridJD, WriteTime) == FAIL) {
      ENZO_FAIL("Error in WriteDataCubes\n");
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
 
  if (((RadiationFieldType >= 10 && RadiationFieldType <= 11) ||
       RadiationData.RadiationShield == TRUE) &&
      MyProcessorNumber == ROOT_PROCESSOR) {
    
    FILE *Radfptr;
 
    strcpy(radiationname, name);
    strcat(radiationname, RadiationSuffix);
 
    if ((Radfptr = fopen(radiationname, "w")) == NULL) {
      ENZO_VFAIL("Error opening radiation file %s\n", radiationname)
    }
    if (WriteRadiationData(Radfptr) == FAIL) {
      ENZO_FAIL("Error in WriteRadiationData\n");
    }
 
    fclose(Radfptr);
 
  }
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if (HierarchyFileOutputFormat > 0)
      fclose(fptr);
    fclose(mptr);
  }

#ifdef TASKMAP
  fclose(tptr);
#endif
 
  // Replace the time in metadata with the saved value (above)
 
  MetaData.Time = SavedTime;

  CommunicationBarrier();
// if (debug)
    //  fprintf(stderr, "WriteAllData: finished writing data\n");
 
//  Stop I/O timing
 
  CommunicationBarrier();
 
  ContinueExecution();
 
  CommunicationBarrier();

  if ( MyProcessorNumber == ROOT_PROCESSOR ){
    sptr = fopen("OutputLogA", "a");
    fprintf(sptr, "DATASET WRITTEN %s %8"ISYM" %18.16e\n", name, MetaData.CycleNumber, MetaData.Time);
    fclose(sptr);
  }
 
  return SUCCESS;
}
 
void DeleteGridHierarchy(HierarchyEntry *GridEntry)
{
  if (GridEntry->NextGridThisLevel != NULL)

     DeleteGridHierarchy(GridEntry->NextGridThisLevel);
 
  delete GridEntry;
 

 
  return;
}
