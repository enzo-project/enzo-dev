/***********************************************************************
/
/  GRID CLASS (READ IN GRID HIERARCHY INFORMATION FROM HDF5 FILE)
/
/  written by: Michael Kuhlen
/  date:       October, 2010
/
/  PURPOSE:
/
************************************************************************/
 
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

void my_exit(int status);
 

//#define IO_LOG

int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, int &Attribute, FILE *log_fptr);
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, float &Attribute, FILE *log_fptr);
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eint32 &Attribute, FILE *log_fptr);
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eint64 &Attribute, FILE *log_fptr);
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eflt32 &Attribute, FILE *log_fptr);
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eflt64 &Attribute, FILE *log_fptr);
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eflt128 &Attribute, FILE *log_fptr);
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, int Attribute[], int NumberOfElements, FILE *log_fptr);
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, char Attribute[], FILE *log_fptr);

int HDF5_ReadDataset(hid_t group_id, const char *DatasetName, int &Dataset, FILE *log_fptr);
int HDF5_ReadDataset(hid_t group_id, const char *DatasetName, int Dataset[], FILE *log_fptr);
int HDF5_ReadDataset(hid_t group_id, const char *DatasetName, FLOAT Dataset[], FILE *log_fptr);


int grid::ReadHierarchyInformationHDF5(hid_t Hfile_id, int GridID, int &Task, int &NextGridThisLevelID, int &NextGridNextLevelID, char DataFilename[], FILE *log_fptr) {
  char GroupName[MAX_LINE_LENGTH];

  int level;

  hid_t       group_id;

  herr_t      h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  
  // ***** Open grid group *****
  level = LevelLookupTable[GridID-1];
  sprintf(GroupName,"Level%"ISYM"/Grid%"GROUP_TAG_FORMAT""ISYM, level, GridID);
  
  group_id = H5Gopen(Hfile_id, GroupName);
  

  // ***** Read attributes *****

  HDF5_ReadAttribute(group_id, "Task", Task, log_fptr);

  //  if ( MyProcessorNumber == 0 )
  //    fprintf(stderr, "Reading Grid %"ISYM" assigned to Task %"ISYM"\n", TestGridID, Task);
  
  HDF5_ReadAttribute(group_id, "GridRank", GridRank, log_fptr);

  HDF5_ReadAttribute(group_id, "Time", Time, log_fptr);

  HDF5_ReadAttribute(group_id, "OldTime", OldTime, log_fptr);

  HDF5_ReadAttribute(group_id, "SubgridsAreStatic", SubgridsAreStatic, log_fptr);

  HDF5_ReadAttribute(group_id, "NumberOfBaryonFields", NumberOfBaryonFields, log_fptr);

  if (NumberOfBaryonFields > 0) {
    HDF5_ReadAttribute(group_id, "FieldType", FieldType, NumberOfBaryonFields, log_fptr);

    HDF5_ReadAttribute(group_id, "BaryonFileName", DataFilename, log_fptr);

    HDF5_ReadAttribute(group_id, "CourantSafetyNumber", CourantSafetyNumber, log_fptr);
    
    HDF5_ReadAttribute(group_id, "PPMFlatteningParameter", PPMFlatteningParameter, log_fptr);
    
    HDF5_ReadAttribute(group_id, "PPMDiffusionParameter", PPMDiffusionParameter, log_fptr);
    
    HDF5_ReadAttribute(group_id, "PPMSteepeningParameter", PPMSteepeningParameter, log_fptr);
  }

  if (SelfGravity) {
    // GravityBoundaryType
    HDF5_ReadAttribute(group_id, "GravityBoundaryType", GravityBoundaryType, log_fptr);
  }
  
  // NextGridThisLevelID
  HDF5_ReadAttribute(group_id, "NextGridThisLevelID", NextGridThisLevelID, log_fptr);
  
  // NextGridNextLevelID
  HDF5_ReadAttribute(group_id, "NextGridNextLevelID", NextGridNextLevelID, log_fptr);



  // ***** Read datasets *****

  HDF5_ReadDataset(group_id, "NumberOfParticles", NumberOfParticles, log_fptr);
  
  if (NumberOfParticles > 0) {
    HDF5_ReadAttribute(group_id, "ParticleFileName", DataFilename, log_fptr);
  }

  HDF5_ReadDataset(group_id, "GridDimension", GridDimension, log_fptr);

  HDF5_ReadDataset(group_id, "GridStartIndex", GridStartIndex, log_fptr);

  HDF5_ReadDataset(group_id, "GridEndIndex", GridEndIndex, log_fptr);

  HDF5_ReadDataset(group_id, "GridRightEdge", GridRightEdge, log_fptr);

  HDF5_ReadDataset(group_id, "GridLeftEdge", GridLeftEdge, log_fptr);


  // If HierarchyFile has different Ghostzones (which should be a parameter not a macro ...)
  // (useful in a restart with different hydro/mhd solvers) 
  int ghosts = NumberOfGhostZones;
  if (GridStartIndex[0] != ghosts)  {
    if (GridID < 2)
      fprintf(stderr,"Grid_ReadHierarchyInformationHDF5: Adjusting Ghostzones which in the HDF5 hierarchy file did not match the selected HydroMethod.\n");
    for (int dim=0; dim < GridRank; dim++) {
      GridDimension[dim]  = GridEndIndex[dim]-GridStartIndex[dim]+1+2*ghosts;
      GridStartIndex[dim] = ghosts;
      GridEndIndex[dim]   = GridStartIndex[dim]+GridDimension[dim]-1-2*ghosts;
      if (GridID < 2) fprintf(stderr, "dim: GridStart,GridEnd,GridDim:  %"ISYM": %"ISYM" %"ISYM" %"ISYM"\n",
			      dim, GridStartIndex[dim], GridEndIndex[dim], GridDimension[dim]);
    }
  }
  

  // ***** Close this grid *****
  h5_status = H5Gclose(group_id);
  if (io_log) fprintf(log_fptr, "H5Gclose: status = %"ISYM"\n", (int) h5_status);

    
  return SUCCESS;
}




// HDF5 utility routines (to read attributes and datasets)

// 32-bit int (Eint32)
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eint32 &Attribute, FILE *log_fptr) {

  hid_t attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  attr_id = H5Aopen_name(group_id, AttributeName);
  if (io_log) fprintf(log_fptr, "H5Aopen_name: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Aread(attr_id, HDF5_I4, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Aread: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// 64-bit int (Eint64)
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eint64 &Attribute, FILE *log_fptr) {

  hid_t attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  attr_id = H5Aopen_name(group_id, AttributeName);
  if (io_log) fprintf(log_fptr, "H5Aopen_name: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Aread(attr_id, HDF5_I8, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Aread: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// 32-bit float (Eflt32)
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eflt32 &Attribute, FILE *log_fptr) {

  hid_t attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif
  attr_id = H5Aopen_name(group_id, AttributeName);
  if (io_log) fprintf(log_fptr, "H5Aopen_name: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Aread(attr_id, HDF5_R4, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Aread: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// 64-bit float (Eflt64)
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eflt64 &Attribute, FILE *log_fptr) {

  hid_t attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  attr_id = H5Aopen_name(group_id, AttributeName);
  if (io_log) fprintf(log_fptr, "H5Aopen_name: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Aread(attr_id, HDF5_R8, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Aread: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// 128-bit float (Eflt128)
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, Eflt128 &Attribute, FILE *log_fptr) {

  hid_t attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  attr_id = H5Aopen_name(group_id, AttributeName);
  if (io_log) fprintf(log_fptr, "H5Aopen_name: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Aread(attr_id, HDF5_R16, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Aread: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// int vector
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, int Attribute[], int NumberOfElements, FILE *log_fptr) {

  hid_t attr_id;
  hsize_t dims[1];

  herr_t h5_status;

  int io_log = 1;
#ifdef IO_LOG
  io_log = 1;
#endif

  attr_id = H5Aopen_name(group_id, AttributeName);
  if (io_log) fprintf(log_fptr, "H5Aopen_name: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Aread(attr_id, HDF5_INT, Attribute);
  if (io_log) fprintf(log_fptr, "H5Aread: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// string (char vector)
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, char Attribute[], FILE *log_fptr) {

  hid_t attr_id, string_type_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  attr_id = H5Aopen_name(group_id, AttributeName);
  if (io_log) fprintf(log_fptr, "H5Aopen_name: attr_id = %"ISYM"\n", (int) attr_id);

  string_type_id = H5Aget_type(attr_id);
  
  h5_status = H5Aread(attr_id, string_type_id, Attribute);
  if (io_log) fprintf(log_fptr, "H5Aread: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}


// int
int HDF5_ReadDataset(hid_t group_id, const char *DatasetName, int &Dataset, FILE *log_fptr) {

  hid_t dset_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dset_id = H5Dopen(group_id, DatasetName);
  if (io_log) fprintf(log_fptr, "H5Dopen: dset_id = %"ISYM"\n", (int) dset_id);

  h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) &Dataset);
  if (io_log) fprintf(log_fptr, "H5Dread: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Dclose(dset_id);
  if (io_log) fprintf(log_fptr, "H5Dclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// int vector
int HDF5_ReadDataset(hid_t group_id, const char *DatasetName, int Dataset[], FILE *log_fptr) {

  hid_t dset_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dset_id = H5Dopen(group_id, DatasetName);
  if (io_log) fprintf(log_fptr, "H5Dopen: dset_id = %"ISYM"\n", (int) dset_id);

  h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) Dataset);
  if (io_log) fprintf(log_fptr, "H5Dread: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Dclose(dset_id);
  if (io_log) fprintf(log_fptr, "H5Dclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// FLOAT vector
int HDF5_ReadDataset(hid_t group_id, const char *DatasetName, FLOAT Dataset[], FILE *log_fptr) {

  hid_t dset_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dset_id = H5Dopen(group_id, DatasetName);
  if (io_log) fprintf(log_fptr, "H5Dopen: dset_id = %"ISYM"\n", (int) dset_id);

  h5_status = H5Dread(dset_id, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) Dataset);
  if (io_log) fprintf(log_fptr, "H5Dread: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Dclose(dset_id);
  if (io_log) fprintf(log_fptr, "H5Dclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}
