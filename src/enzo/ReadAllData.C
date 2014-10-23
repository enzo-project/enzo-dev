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
/  modified4:  Michael Kuhlen, October 2010, HDF5 hierarchy
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
#include "CommunicationUtilities.h"
void my_exit(int status);
 
/* function prototypes */
 
int ReadDataHierarchy(FILE *fptr, hid_t Hfile_id, HierarchyEntry *TopGrid, int GridID, HierarchyEntry *ParentGrid, FILE *log_fptr);
int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt);
int ReadStarParticleData(FILE *fptr, hid_t Hfile_id, FILE *log_fptr);
int ReadRadiationData(FILE *fptr);
int AssignGridToTaskMap(Eint64 *GridIndex, Eint64 *Mem, int Ngrids);
 
extern char RadiationSuffix[];
extern char HierarchySuffix[];
extern char hdfsuffix[];
extern char TaskMapSuffix[];
extern char MemoryMapSuffix[]; 
 
//#define IO_LOG
#ifdef IO_LOG
int io_log = 1;
#else
int io_log = 0;
#endif

// the following HDF5 helper routines are defined in
// Grid_ReadHierarchyInformationHDF5.C
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, int &Attribute, FILE *log_fptr);
int HDF5_ReadDataset(hid_t group_id, const char *DatasetName, int Dataset[], FILE *log_fptr);

int ReadAllData(char *name, HierarchyEntry *TopGrid, TopGridData &MetaData,
		ExternalBoundary *Exterior, float *Initialdt)
 
{
 
  /* declarations */
 
  char pid[MAX_TASK_TAG_SIZE];
  char hierarchyname[MAX_LINE_LENGTH], radiationname[MAX_LINE_LENGTH];
  char HDF5hierarchyname[MAX_LINE_LENGTH];
  // Code shrapnel. See comments below. --Rick
  // char taskmapname[MAX_LINE_LENGTH];
  char memorymapname[MAX_LINE_LENGTH];

  FILE *log_fptr;
  FILE *fptr;
  FILE *tptr;
  FILE *mptr;

  hid_t Hfile_id;
  herr_t h5_status;
  herr_t h5_error = -1;

  int GridID = 1;
  int GridKD = 1;

  float dummy;

  // store the original parameter file name, in case we need it later
  strcpy(PrevParameterFileName, name);

  CommunicationBarrier();

  /* Read TopGrid data. */
 
  if ((fptr = fopen(name, "r")) == NULL) {
    ENZO_VFAIL("Error opening input file %s.\n", name)
  }
  if (ReadParameterFile(fptr, MetaData, Initialdt) == FAIL) {
        ENZO_FAIL("Error in ReadParameterFile.");
  }
 
  /* Close main file. */
  fprintf(stderr, "fclose: opening boundary condition file: %s\n", MetaData.BoundaryConditionName);
 
  fclose(fptr);
 
  /* Set the number of particle attributes, if left unset. */

  if (NumberOfParticleAttributes == INT_UNDEFINED ||
      NumberOfParticleAttributes == 0) {
    if (StarParticleCreation || StarParticleFeedback) {
      NumberOfParticleAttributes = 3;
      if (StarMakerTypeIaSNe) NumberOfParticleAttributes++;
      if (StarMakerTypeIISNeMetalField) NumberOfParticleAttributes++;
      AddParticleAttributes = TRUE;
    } else {
      NumberOfParticleAttributes = 0;
    }

  }

  /* Read Boundary condition info. */
  fprintf(stderr, "fopen: opening boundary condition file: %s\n", MetaData.BoundaryConditionName);
 
  int BRerr = 0 ;
  if ((fptr = fopen(MetaData.BoundaryConditionName, "r")) == NULL) {
    fprintf(stderr, "Error opening boundary condition file: %s\n",
	    MetaData.BoundaryConditionName);
    BRerr = 1;
  } 

  // Try to read external boundaries. If they don't fit grid data we'll set them later below
  int ReadHDF4B = 0;
#ifdef USE_HDF4
  if (Exterior->ReadExternalBoundaryHDF4(fptr) == FAIL) {  
    fprintf(stderr, "Error in ReadExternalBoundary using HDF4  (%s).\n",           
	    MetaData.BoundaryConditionName);                  
    fprintf(stderr, "Will try HDF5 instead.\n");
    BRerr = 1;
  } else ReadHDF4B = 1;
#endif // HDF4
  if (ReadHDF4B != 1) 
    if(LoadGridDataAtStart){    
      if (Exterior->ReadExternalBoundary(fptr) == FAIL) {
	fprintf(stderr, "Error in ReadExternalBoundary (%s).\n",
		MetaData.BoundaryConditionName);
	BRerr = 1;
      }
    }else{
      if (Exterior->ReadExternalBoundary(fptr, TRUE, FALSE) == FAIL) {
	fprintf(stderr, "Error in ReadExternalBoundary (%s).\n",
		MetaData.BoundaryConditionName);
	BRerr = 1;
      }
    }


  if (BRerr ==0) strcat(MetaData.BoundaryConditionName, hdfsuffix);
  if ((fptr != NULL) && (BRerr == 0)) fclose(fptr);

  /* Create the memory map name */

  strcpy(memorymapname, name);
  strcat(memorymapname, MemoryMapSuffix);

#ifdef USE_MPI
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
#endif
 
  /* Read the memory map */

  CommunicationBarrier();

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
  }

  /* Read Data Hierarchy. */ 

  GridID = 1;
  if (ReadDataHierarchy(fptr, Hfile_id, TopGrid, GridID, NULL, log_fptr) == FAIL) {
    ENZO_VFAIL("Error in ReadDataHierarchy (%s).\n", hierarchyname)
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
    
    SimpleConstantBoundary = TRUE;
    
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      if (MetaData.LeftFaceBoundaryCondition[dim] != periodic ||
	  MetaData.RightFaceBoundaryCondition[dim] != periodic) {
	SimpleConstantBoundary = FALSE;
      }
    }
    
    Exterior->Prepare(TopGrid->GridData);
    
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      Exterior->InitializeExternalBoundaryFace(dim,
					       MetaData.LeftFaceBoundaryCondition[dim],
					       MetaData.RightFaceBoundaryCondition[dim],
					       Dummy, Dummy);
      fprintf(stderr, " %i  %i \n", MetaData.LeftFaceBoundaryCondition[dim],
	      MetaData.RightFaceBoundaryCondition[dim]);
    }
    
    Exterior->InitializeExternalBoundaryParticles(
						  MetaData.ParticleBoundaryType);
    
    
    /*
    HierarchyEntry *NextGrid;
    NextGrid = TopGrid->NextGridThisLevel;
    while (NextGrid != NULL) {
      Exterior->Prepare(NextGrid->GridData);
      NextGrid->GridData->SetExternalBoundaryValues(Exterior);
      NextGrid = NextGrid->NextGridThisLevel;
    }
    */
    
  } /* Setting External Boundary */

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

  if (io_log)
    fclose(log_fptr);

  return SUCCESS;
}
