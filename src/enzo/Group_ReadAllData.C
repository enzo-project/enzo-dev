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
/  modified5:  Michael Kuhlen, October 2010, HDF5 hierarchy
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
#include "CommunicationUtilities.h"
void my_exit(int status);

// HDF5 function prototypes


/* function prototypes */
 
int Group_ReadDataHierarchy(FILE *fptr, hid_t Hfile_id, HierarchyEntry *TopGrid, int GridID,
			    HierarchyEntry *ParentGrid, hid_t file_id,
			    int NumberOfRootGrids, int *RootGridProcessors,
			    bool ReadParticlesOnly=false, FILE *log_fptr=NULL);
int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt);
int ReadStarParticleData(FILE *fptr, hid_t Hfile_id, FILE *log_fptr);
int ReadRadiationData(FILE *fptr);
int AssignGridToTaskMap(Eint64 *GridIndex, Eint64 *Mem, int Ngrids);
int InitialLoadBalanceRootGrids(FILE *fptr, hid_t Hfile_id, int TopGridRank,
				int TopGridDim, int &NumberOfRootGrids,
				int* &RootProcessors);
 
extern char RadiationSuffix[];
extern char HierarchySuffix[];
extern char hdfsuffix[];
extern char TaskMapSuffix[];
extern char MemoryMapSuffix[];
extern char CPUSuffix[];


#define IO_LOG

// the following HDF5 helper routines are defined in
// Grid_ReadHierarchyInformationHDF5.C
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, int &Attribute, FILE *log_fptr);
int HDF5_ReadDataset(hid_t group_id, const char *DatasetName, int Dataset[], FILE *log_fptr);

int Group_ReadAllData(char *name, HierarchyEntry *TopGrid, TopGridData &MetaData,
		      ExternalBoundary *Exterior, float *Initialdt,
		      bool ReadParticlesOnly)
 
{
 
  /* declarations */
 
  char hierarchyname[MAX_LINE_LENGTH], radiationname[MAX_LINE_LENGTH];
  char HDF5hierarchyname[MAX_LINE_LENGTH];
  // Code shrapnel. See comments below. --Rick
  //  char taskmapname[MAX_LINE_LENGTH];
  char memorymapname[MAX_LINE_LENGTH];
  char groupfilename[MAX_LINE_LENGTH];
  char line[MAX_LINE_LENGTH];

  FILE *log_fptr;
  FILE *fptr;
  FILE *tptr;
  FILE *mptr;

  hid_t       file_id, Hfile_id;
  hid_t       attr_id, dset_id;
  hid_t       file_acc_template;
  size_t      memory_increment; // in bytes
  hbool_t     dump_flag;

  herr_t      h5_status;
  herr_t      h5_error = -1;

  int GridID = 1;
  int GridKD = 1;


  float dummy;
  int dummy_int;

  int *NumberOfSubgridCells;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  // store the original parameter file name, in case we need it later
  strcpy(PrevParameterFileName, name);

  char pid[MAX_TASK_TAG_SIZE];

  CommunicationBarrier();
 
  /* Read TopGrid data. */
 
  if ((fptr = fopen(name, "r")) == NULL) {
    ENZO_VFAIL("Error opening input file %s.\n", name)
  }
  if (ReadParameterFile(fptr, MetaData, Initialdt) == FAIL) {
        ENZO_FAIL("Error in ReadParameterFile.");
  }
 
  /* Close main file. */
 
  fclose(fptr);

  /* Set the number of particle attributes, if left unset. */

  if (NumberOfParticleAttributes == INT_UNDEFINED ||
      NumberOfParticleAttributes == 0) {
    if (StarParticleCreation || StarParticleFeedback) {
      NumberOfParticleAttributes = 3;
      if (StarMakerTypeIaSNe) NumberOfParticleAttributes++;
      AddParticleAttributes = TRUE;
    } else {
      NumberOfParticleAttributes = 0;
    }

  }

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
  int BRerr=0 ;
  if ((fptr = fopen(MetaData.BoundaryConditionName, "r")) == NULL) {
    fprintf(stderr, "Error opening boundary condition file: %s\n",
	    MetaData.BoundaryConditionName);
    return FAIL;
  }

  // Below, ENZO_FAIL is changed to "return FAIL" to deal with various data formats including HDF4, HDF5, packed-HDF5
  // because boundary should be the one that distinguishes these different data formats.
  // This will allow a graceful exit when the dataformat is not packed-HDF5.
  // - Ji-hoon Kim
  // Try to read external boundaries. If they don't fit grid data we'll set them later below
    if(LoadGridDataAtStart){    
      if (Exterior->ReadExternalBoundary(fptr) == FAIL) {
	fprintf(stderr, "Error in ReadExternalBoundary (%s).\n",
		MetaData.BoundaryConditionName);
	return FAIL;
      }
    }else{
      if (Exterior->ReadExternalBoundary(fptr, TRUE, FALSE) == FAIL) {
	fprintf(stderr, "Error in ReadExternalBoundary (%s).\n",
		MetaData.BoundaryConditionName);
	return FAIL;
      }
    }

  strcat(MetaData.BoundaryConditionName, hdfsuffix);
  if (fptr != NULL) fclose(fptr);

  /* Create the memory map name */

  strcpy(memorymapname, name);
  strcat(memorymapname, MemoryMapSuffix);
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  /* Read the memory map */

  CommunicationBarrier();

  // fprintf(stderr, "All at sync point - read MemoryMap and assign Grid->Task\n");

#ifdef TASKMAP
  if ((mptr = fopen(memorymapname, "r")) == NULL) {
    ENZO_VFAIL("Error opening MemoryMap file %s.\n", memorymapname)
  }

  Eint64 GridIndex[MAX_NUMBER_OF_TASKS], OldPN, Mem[MAX_NUMBER_OF_TASKS];
  int i, ntask;
  i = 0;

  while( fscanf(mptr, "Grid %8lld  PN %8lld  Memory %16lld\n", &GridIndex[i], &OldPN, &Mem[i]) != EOF) {
    i++;
  }
  ntask = i;

  if (AssignGridToTaskMap(GridIndex, Mem, ntask) == FAIL) {
        ENZO_FAIL("Error in AssignGridToTaskMap.");
  }

  fclose(mptr);
#endif

  /* Create hierarchy name and grid base name. */

  strcpy(hierarchyname, name);
  strcat(hierarchyname, HierarchySuffix);

  //  if (HierarchyFileInputFormat == 0) {
  if (HierarchyFileInputFormat % 2 == 0) {
    sprintf(HDF5hierarchyname,"%s.hdf5",hierarchyname);

    if (io_log) {
      char logname[MAX_LINE_LENGTH];
      sprintf(logname,"%s.in_log",HDF5hierarchyname);
      
      log_fptr = fopen(logname,"w");
    }
    
    Hfile_id = H5Fopen(HDF5hierarchyname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if( Hfile_id == h5_error )
      ENZO_VFAIL("Error opening HDF5 hierarchy file: %s",HDF5hierarchyname)

    // read NumberOfProcessors attribute
    HDF5_ReadAttribute(Hfile_id, "NumberOfProcessors", PreviousMaxTask, log_fptr);
    PreviousMaxTask--;
    
    // read TotalNumberOfGrids attribute
    HDF5_ReadAttribute(Hfile_id, "TotalNumberOfGrids", TotalNumberOfGrids, log_fptr);

    // read LevelLookupTable
    LevelLookupTable = new int[TotalNumberOfGrids];
    HDF5_ReadDataset(Hfile_id, "LevelLookupTable", LevelLookupTable, log_fptr);

//     if(MyProcessorNumber == ROOT_PROCESSOR)
//       for(int i=0;i<TotalNumberOfGrids;i++)
// 	fprintf(stderr,"LevelLookupTable[%d] = %d\n",i,LevelLookupTable[i]);

  } 

  if (HierarchyFileInputFormat == 1) {
    if ((fptr = fopen(hierarchyname, "r")) == NULL) {
      ENZO_VFAIL("Error opening hierarchy file %s.\n", hierarchyname)
    }

    /* scan data hierarchy for maximum task number */
    
    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL)
      if (sscanf(line, "Task = %"ISYM, &dummy_int) > 0)
	PreviousMaxTask = max(PreviousMaxTask, dummy_int);

    rewind(fptr);

  }


  /* If we're load balancing only within nodes, count level-1 cells in
     each level-0 grid and load balance the entire nodes. */

  if (ResetLoadBalancing)
    LoadBalancing = 1;

  int *RootGridProcessors = NULL, NumberOfRootGrids = 1;
  InitialLoadBalanceRootGrids(fptr, Hfile_id, MetaData.TopGridRank, MetaData.TopGridDims[0], NumberOfRootGrids, RootGridProcessors);

  /* Read Data Hierarchy. */

  if(LoadGridDataAtStart){
    // open of HDF5 file here with incore properties

    CommunicationBarrier();

#ifdef SINGLE_HDF5_OPEN_ON_INPUT

    fprintf(stderr, "OPEN %s on processor %"ISYM"\n", groupfilename, MyProcessorNumber);

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

    if (debug && HierarchyFileInputFormat == 1)
      fprintf(stdout, "OPEN data hierarchy %s\n", hierarchyname);
    if (debug && HierarchyFileInputFormat % 2 == 0)
      fprintf(stdout, "OPEN data hierarchy %s\n", HDF5hierarchyname);
    file_id = h5_error;

#endif /* SINGLE OPEN */
  }

  GridID = 1;
  if (Group_ReadDataHierarchy(fptr, Hfile_id, TopGrid, GridID, NULL, file_id,
			      NumberOfRootGrids, RootGridProcessors,
			      ReadParticlesOnly, log_fptr) == FAIL) {
    fprintf(stderr, "Error in ReadDataHierarchy (%s).\n", hierarchyname);
    return FAIL;
  }

//   printf("P%d: out of Group_RDH\n", MyProcessorNumber);
//   CommunicationBarrier();
  
  if(LoadGridDataAtStart){
    // can close HDF5 file here

    if(CheckpointRestart == TRUE) {
#ifndef SINGLE_HDF5_OPEN_ON_INPUT

    file_id = H5Fopen(groupfilename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(file_id == h5_error)ENZO_VFAIL("Could not open %s", groupfilename)

#endif
      // Now we load our metadata back in
      hid_t metadata_group;
    H5E_BEGIN_TRY{
      metadata_group = H5Gopen(file_id, "Metadata");
    }H5E_END_TRY
    if(metadata_group != h5_error) {
      readAttribute(metadata_group, HDF5_INT, "LevelCycleCount",
          LevelCycleCount, TRUE);
      if(CheckpointRestart == TRUE) { // We only need these in a checkpoint
        FLOAT dtThisLevelCopy[MAX_DEPTH_OF_HIERARCHY];
        FLOAT dtThisLevelSoFarCopy[MAX_DEPTH_OF_HIERARCHY];
        readAttribute(metadata_group, HDF5_PREC, "dtThisLevel",
            dtThisLevelCopy, TRUE);
        readAttribute(metadata_group, HDF5_PREC, "dtThisLevelSoFar",
            dtThisLevelSoFarCopy, TRUE);
        for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
            dtThisLevel[level] = dtThisLevelCopy[level];
            dtThisLevelSoFar[level] = dtThisLevelSoFarCopy[level];
        }
        readAttribute(metadata_group, HDF5_PREC, "Time",
            &MetaData.Time, TRUE);
      }
    } else if(CheckpointRestart == TRUE) {
      ENZO_FAIL("Couldn't open Metadata!");
    }
    H5Gclose(metadata_group);

#ifndef SINGLE_HDF5_OPEN_ON_INPUT
      H5Fclose(file_id);
#endif
    }
#ifdef SINGLE_HDF5_OPEN_ON_INPUT

    h5_status = H5Fclose(file_id);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#ifdef USE_HDF5_INPUT_BUFFERING

    h5_status = H5Pclose(file_acc_template);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

#endif

#endif /* SINGLE OPEN */
  }


  /* If there was trouble reading the boundary file attempt to sanely set them now. */
  /* We do this after all the grids have been read in case the number of baryon fields 
     etc. was modified in that stage. 
     This allows you to add extra fields etc. and then restart. Proceed with caution.
  */
  if (BRerr) {
    fprintf(stderr,"Setting External Boundaries.\n");
  float Dummy[MAX_DIMENSION];
  int dim;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Dummy[dim] = 0.0;
  fprintf(stderr, "Because of trouble in reading the boundary we reset it now.\n");
  Exterior->Prepare(TopGrid->GridData);
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    Exterior->InitializeExternalBoundaryFace(dim,
					     MetaData.LeftFaceBoundaryCondition[dim],
					     MetaData.RightFaceBoundaryCondition[dim],
					     Dummy, Dummy);
  TopGrid->GridData->SetExternalBoundaryValues(Exterior);
  }

  /* Read StarParticle data. */
 
  if (ReadStarParticleData(fptr, Hfile_id, log_fptr) == FAIL) {
        ENZO_FAIL("Error in ReadStarParticleData.");
  }
 
  /* Create radiation name and read radiation data. */
 
  if ((RadiationFieldType >= 10 && RadiationFieldType <= 11) || 
      (RadiationData.RadiationShield == TRUE && RadiationFieldType != 12)) {
    FILE *Radfptr;
    strcpy(radiationname, name);
    strcat(radiationname, RadiationSuffix);
    if ((Radfptr = fopen(radiationname, "r")) == NULL) {
      ENZO_VFAIL("Error opening radiation file %s.\n", name)
    }
    if (ReadRadiationData(Radfptr) == FAIL) {
            ENZO_FAIL("Error in ReadRadiationData.");
    }
    fclose(Radfptr);
  }
 
  //  if (HierarchyFileInputFormat == 0) {
  if (HierarchyFileInputFormat % 2 == 0) {
    h5_status = H5Fclose(Hfile_id);
    if (h5_status == h5_error)
      ENZO_FAIL("Error closing HDF5 hierarchy file.");    

    delete [] LevelLookupTable;
  }

  if (HierarchyFileInputFormat == 1)
    fclose(fptr);

  /* If we added new particle attributes, unset flag so we don't carry
     this parameter to later data. */

  AddParticleAttributes = FALSE;

  /* If we're reseting load balancing (i.e. the host processors), turn
     off the reset flag because we've already done this and don't want
     it to propagate to later datadumps. */

  if (ResetLoadBalancing)

    ResetLoadBalancing = FALSE;

  delete [] RootGridProcessors;

  if (HierarchyFileInputFormat % 2 == 0 && io_log) 
      fclose(log_fptr);
  

  return SUCCESS;
}
