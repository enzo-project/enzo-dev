/***********************************************************************
/
/  READ ALL THE DATA (DATA & RESTART DUMP)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness
/              January, 2006
/  modified2:  Robert Harkness
/              April, 2006
/              I/O performance on input
/  modified3:  Robert Harkness
/              January, 2007
/              Read group file in-core
/  modified4:  Robert Harkness
/              April 2008
/
/  PURPOSE:
/
************************************************************************/
 
// This function reads in the data hierarchy (TopGrid), the External
//   Boundary (Exterior), the TopGridData, and the global_data.

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <hdf5.h> 
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
#include "CommunicationUtilities.h"
void my_exit(int status);

// HDF5 function prototypes


 
/* function prototypes */
 
int Group_ReadDataHierarchy(FILE *fptr, HierarchyEntry *TopGrid, int GridID,
			    HierarchyEntry *ParentGrid, hid_t file_id);
int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt);
int ReadStarParticleData(FILE *fptr);
int ReadRadiationData(FILE *fptr);
int AssignGridToTaskMap(Eint64 *GridIndex, Eint64 *Mem, int Ngrids);
 
extern char RadiationSuffix[];
extern char HierarchySuffix[];
extern char hdfsuffix[];
extern char TaskMapSuffix[];
extern char MemoryMapSuffix[];
extern char CPUSuffix[];

 
int Group_ReadAllData(char *name, HierarchyEntry *TopGrid, TopGridData &MetaData,
		      ExternalBoundary *Exterior)
 
{
 
  /* declarations */
 
  char hierarchyname[MAX_LINE_LENGTH], radiationname[MAX_LINE_LENGTH];
  // Code shrapnel. See comments below. --Rick
  //  char taskmapname[MAX_LINE_LENGTH];
  char memorymapname[MAX_LINE_LENGTH];
  char groupfilename[MAX_LINE_LENGTH];

  FILE *fptr;
  FILE *tptr;
  FILE *mptr;

  hid_t       file_id;
  hid_t       file_acc_template;
  size_t      memory_increment; // in bytes
  hbool_t     dump_flag;

  herr_t      h5_status;
  herr_t      h5_error = -1;

  int GridID = 1;
  int GridKD = 1;

  float dummy;

  // store the original parameter file name, in case we need it later
  strcpy(PrevParameterFileName, name);

#ifdef USE_MPI
  double io_start, io_stop;
  char io_logfile[MAX_NAME_LENGTH];
  FILE *xptr;
#endif /* USE_MPI */
  char pid[MAX_TASK_TAG_SIZE];

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


  // name is something like /dsgpfs/harkness/NewL7/Dumps/DD0156/DD0156
  // open the hdf file on this processor /dsgpfs/harkness/NewL7/Dumps/DD0156/DD0156.cpu0000, etc.
  // the task map should respect this, otherwise the map is scrambled on input

  // get the task number

  MPI_Arg mpi_rank;

#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#else
  mpi_rank = 0;
#endif

  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);

  strcpy(groupfilename, name);
  strcat(groupfilename, CPUSuffix);
  strcat(groupfilename, pid);

 
  /* Read Boundary condition info. */
 
  if ((fptr = fopen(MetaData.BoundaryConditionName, "r")) == NULL) {
    fprintf(stderr, "Error opening boundary condition file: %s\n",
	    MetaData.BoundaryConditionName);
    ENZO_FAIL("");
  }
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
  strcat(MetaData.BoundaryConditionName, hdfsuffix);
  fclose(fptr);

  /* Create the memory map name */

  strcpy(memorymapname, name);
  strcat(memorymapname, MemoryMapSuffix);
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  // Code shrapnel. taskmapname is used in the commented block below.
  // --Rick
  //   strcpy(taskmapname, name);
  //   strcat(taskmapname, TaskMapSuffix);
  //   strcat(taskmapname, pid);

  /* Read the memory map */

  CommunicationBarrier();

  // fprintf(stderr, "All at sync point - read MemoryMap and assign Grid->Task\n");

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


  if(LoadGridDataAtStart){
    // open of HDF5 file here with incore properties

    fprintf(stderr, "OPEN %s on processor %"ISYM"\n", groupfilename, MyProcessorNumber);

    CommunicationBarrier();

#ifdef SINGLE_HDF5_OPEN_ON_INPUT

#ifdef USE_HDF5_INPUT_BUFFERING

    memory_increment = 1024*1024;
    dump_flag = 0;

    file_acc_template = H5Pcreate (H5P_FILE_ACCESS);
    if( file_acc_template == h5_error ){my_exit(EXIT_FAILURE);}

    h5_status = H5Pset_fapl_core(file_acc_template, memory_increment, dump_flag);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

    file_id = H5Fopen(groupfilename, H5F_ACC_RDWR, file_acc_template);
    if( file_id == h5_error ){my_exit(EXIT_FAILURE);}

#else

    file_id = H5Fopen(groupfilename, H5F_ACC_RDWR, H5P_DEFAULT);
    if( file_id == h5_error ){my_exit(EXIT_FAILURE);}

#endif

#else

    file_id = h5_error;

#endif /* SINGLE OPEN */
  }

  GridID = 1;
  if (Group_ReadDataHierarchy(fptr, TopGrid, GridID, NULL, file_id) == FAIL) {
    fprintf(stderr, "Error in ReadDataHierarchy (%s).\n", hierarchyname);
    ENZO_FAIL("");
  }

  if(LoadGridDataAtStart){
    // can close HDF5 file here

#ifdef SINGLE_HDF5_OPEN_ON_INPUT

    h5_status = H5Fclose(file_id);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#ifdef USE_HDF5_INPUT_BUFFERING

    h5_status = H5Pclose(file_acc_template);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#endif

#endif /* SINGLE OPEN */
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
