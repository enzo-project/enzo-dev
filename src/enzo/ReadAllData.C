/***********************************************************************
/
/  READ ALL THE DATA (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness
/              January, 2006
/  modified2:  Robert harkness
/              April, 2006
/              I/O performance on input
/  modified3:  Robert Harkness
/              April 2008
/
/  PURPOSE:
/
************************************************************************/
 
// This function reads in the data hierarchy (TopGrid), the External
//   Boundary (Exterior), the TopGridData, and the global_data.
 
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
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
#include "StarParticleData.h"
#include "CommunicationUtilities.h"
void my_exit(int status);
 
/* function prototypes */
 
int ReadDataHierarchy(FILE *fptr, HierarchyEntry *TopGrid, int GridID,
		      HierarchyEntry *ParentGrid);
int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt);
int ReadStarParticleData(FILE *fptr);
int ReadRadiationData(FILE *fptr);
int AssignGridToTaskMap(Eint64 *GridIndex, Eint64 *Mem, int Ngrids);
 
extern char RadiationSuffix[];
extern char HierarchySuffix[];
extern char hdfsuffix[];
extern char TaskMapSuffix[];
extern char MemoryMapSuffix[]; 


int ReadAllData(char *name, HierarchyEntry *TopGrid, TopGridData &MetaData,
		 ExternalBoundary *Exterior)
 
{
 
  /* declarations */
 
  char hierarchyname[MAX_LINE_LENGTH], radiationname[MAX_LINE_LENGTH];
  // Code shrapnel. See comments below. --Rick
  // char taskmapname[MAX_LINE_LENGTH];
  char memorymapname[MAX_LINE_LENGTH];

  FILE *fptr;
  FILE *tptr;
  FILE *mptr;

  int GridID = 1;
  int GridKD = 1;

  float dummy;

  // store the original parameter file name, in case we need it later
  strcpy(PrevParameterFileName, name);

#ifdef USE_MPI
  double io_start, io_stop;
  char io_logfile[MAX_NAME_LENGTH];
  char pid[MAX_TASK_TAG_SIZE];
  FILE *xptr;
#endif /* USE_MPI */

//  Start I/O timing

#ifdef USE_MPI
  CommunicationBarrier();
  io_start = MPI_Wtime();
#endif /* USE_MPI */
 
  /* Read TopGrid data. */
 
  if ((fptr = fopen(name, "r")) == NULL) {
    fprintf(stderr, "Error opening input file %s.\n", name);
    ENZO_FAIL("");
  }
  if (ReadParameterFile(fptr, MetaData, &dummy) == FAIL) {
    fprintf(stderr, "Error in ReadParameterFile.\n");
    ENZO_FAIL("");
  }
 
  /* Close main file. */
 
  fclose(fptr);
 
  /* Read Boundary condition info. */
 
  if ((fptr = fopen(MetaData.BoundaryConditionName, "r")) == NULL) {
    fprintf(stderr, "Error opening boundary condition file: %s\n",
	    MetaData.BoundaryConditionName);
    ENZO_FAIL("");
  }

#ifdef USE_HDF4
  if (Exterior->ReadExternalBoundaryHDF4(fptr) == FAIL) {  
    fprintf(stderr, "Error in ReadExternalBoundary (%s).\n",           
            MetaData.BoundaryConditionName);                  
    return FAIL;                                                                 
  }
#else
  if(LoadGridDataAtStart){    
    if (Exterior->ReadExternalBoundary(fptr) == FAIL) {
      fprintf(stderr, "Error in ReadExternalBoundary (%s).\n",
	      MetaData.BoundaryConditionName);
      ENZO_FAIL("");
    }
  }else{
    if (Exterior->ReadExternalBoundary(fptr, TRUE, FALSE) == FAIL) {
      fprintf(stderr, "Error in ReadExternalBoundary (%s).\n",
	      MetaData.BoundaryConditionName);
      ENZO_FAIL("");
    }
  }
#endif

  strcat(MetaData.BoundaryConditionName, hdfsuffix);
  fclose(fptr);

  /* Create the memory map name */

  strcpy(memorymapname, name);
  strcat(memorymapname, MemoryMapSuffix);

#ifdef USE_MPI
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
#endif
 
  // Code shrapnel. taskmapname is used in the commented block below.
  // --Rick
  //   strcpy(taskmapname, name);
  //   strcat(taskmapname, TaskMapSuffix);
  //   strcat(taskmapname, pid);

  /* Read the memory map */

  CommunicationBarrier();

#ifdef TASKMAP
  if ((mptr = fopen(memorymapname, "r")) == NULL) {
    fprintf(stderr, "Error opening MemoryMap file %s.\n", memorymapname);
    ENZO_FAIL("");
  }

  Eint64 GridIndex[MAX_NUMBER_OF_TASKS], OldPN, Mem[MAX_NUMBER_OF_TASKS];
  int i, ntask;
  i = 0;

  while( fscanf(mptr, "Grid %8lld  PN %8lld  Memory %16lld\n", &GridIndex[i], &OldPN, &Mem[i]) != EOF) {
    // fprintf(stderr, "Grid %8lld  PN %8lld  Memory %16lld\n", GridIndex[i], OldPN, Mem[i]);
    i++;
  }
  ntask = i;

  if (AssignGridToTaskMap(GridIndex, Mem, ntask) == FAIL) {
    fprintf(stderr, "Error in AssignGridToTaskMap.\n");
    ENZO_FAIL("");
  }

  fclose(mptr);
#endif


  // This is code shrapnel. It's only purpose is to print out the number of grids
  // that are expected to be opened. It also crashes if the taskmap file isn't found.
  // The variable ngrids isn't used any where below.
  // --Rick
  
  //   /* Count the number of grids to read */

  // #ifdef SINGLE_HDF5_OPEN_ON_INPUT
  //   if ((tptr = fopen(taskmapname, "r")) == NULL) {
  //     fprintf(stderr, "Error opening TaskMap file %s.\n", taskmapname);
  //     ENZO_FAIL("");
  //   }
  
  //   Eint64 OldPN;
  //   int OldGR, ngrids;
  //   int i = 0;
  
  //   while( fscanf(tptr, "Grid %8lld  PN %8lld\n", &OldGR, &OldPN) != EOF) {
  //     // fprintf(stderr, "Grid %8lld  PN %8lld\n", OldGR, OldPN);
  //     i++;
  //   }
  //   ngrids = i;
  
  //   fprintf(stderr, "CPU %"ISYM" owns %"ISYM" grids\n", MyProcessorNumber, ngrids);
  
  //   fclose(tptr);
  // #endif

  /* Create hierarchy name and grid base name. */

  strcpy(hierarchyname, name);
  strcat(hierarchyname, HierarchySuffix);
 
  /* Read Data Hierarchy. */
 
  if ((fptr = fopen(hierarchyname, "r")) == NULL) {
    fprintf(stderr, "Error opening hierarchy file %s.\n", hierarchyname);
    ENZO_FAIL("");
  }
  GridID = 1;
  if (ReadDataHierarchy(fptr, TopGrid, GridID, NULL) == FAIL) {
    fprintf(stderr, "Error in ReadDataHierarchy (%s).\n", hierarchyname);
    ENZO_FAIL("");
  }

  /* Read StarParticle data. */
 
  if (ReadStarParticleData(fptr) == FAIL) {
    fprintf(stderr, "Error in ReadStarParticleData.\n");
    ENZO_FAIL("");
  }
 
  /* Create radiation name and read radiation data. */
 
  if (RadiationFieldType >= 10 && RadiationFieldType <= 11) {
    FILE *Radfptr;
    strcpy(radiationname, name);
    strcat(radiationname, RadiationSuffix);
    if ((Radfptr = fopen(radiationname, "r")) == NULL) {
      fprintf(stderr, "Error opening radiation file %s.\n", name);
      ENZO_FAIL("");
    }
    if (ReadRadiationData(Radfptr) == FAIL) {
      fprintf(stderr, "Error in ReadRadiationData.\n");
      ENZO_FAIL("");
    }
    fclose(Radfptr);
  }
 
  fclose(fptr);

  /* If we added new particle attributes, unset flag so we don't carry
     this parameter to later data. */

  AddParticleAttributes = FALSE;

//  Stop I/O timing

#ifdef USE_MPI
  io_stop = MPI_Wtime();
#endif /* USE_MPI */

#ifdef USE_MPI
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
  strcpy(io_logfile, "IN_perf.");
  strcat(io_logfile, pid);
  xptr = fopen(io_logfile, "a");
  fprintf(xptr, "IN %12.4e  %s\n", (io_stop-io_start), name);
  fclose(xptr);
#endif /* USE_MPI */
 
  return SUCCESS;
}
