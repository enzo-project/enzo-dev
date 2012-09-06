/***********************************************************************
/
/  GRID CLASS (WRITE OUT GRID HIERARCHY INFORMATION IN HDF5)
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
 

// from HDF5 1.8.7+  (H5_VERSION_GE, H5_VERSION_LE)
/* macros for comparing the version */
#define HDF5_VERSION_GE(Maj,Min,Rel) \
       (((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR==Min) && (H5_VERS_RELEASE>=Rel)) || \
        ((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR>Min)) || \
        (H5_VERS_MAJOR>Maj))

#define HDF5_VERSION_LE(Maj,Min,Rel) \
       (((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR==Min) && (H5_VERS_RELEASE<=Rel)) || \
        ((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR<Min)) || \
        (H5_VERS_MAJOR<Maj))



// set this if you want to soft links to the daughter and parent
// grids, and external links to the data file
#define WITH_HDF5_LINKS


//#define IO_LOG

int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eint32 Attribute, FILE *log_fptr);
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eint64 Attribute, FILE *log_fptr);
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eflt32 Attribute, FILE *log_fptr);
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eflt64 Attribute, FILE *log_fptr);
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eflt128 Attribute, FILE *log_fptr);
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, int *Attribute, int NumberOfElements, FILE *log_fptr);
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, char *Attribute, FILE *log_fptr);

int HDF5_WriteDataset(hid_t group_id, const char *DatasetName, int Dataset, FILE *log_fptr);
int HDF5_WriteDataset(hid_t group_id, const char *DatasetName, int *Dataset, int NumberOfElements, FILE *log_fptr);
int HDF5_WriteDataset(hid_t group_id, const char *DatasetName, FLOAT *Dataset, int NumberOfElements, FILE *log_fptr);


int grid::WriteHierarchyInformationHDF5(char *base_name, hid_t level_group_id, int level, int ParentGridIDs[], int NumberOfDaughterGrids, int DaughterGridIDs[], int NextGridThisLevelID, int NextGridNextLevelID, FILE *log_fptr) {
  char GroupName[MAX_LINE_LENGTH];
  char BaryonFileName[MAX_LINE_LENGTH], *ParticleFileName = BaryonFileName;
  char LinkName[MAX_LINE_LENGTH], TargetName[MAX_LINE_LENGTH];

  int dim;
  int GridGlobalPosition[MAX_DIMENSION];
  
  hid_t       group_id, subgroup_id;

  herr_t      h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif


  sprintf(BaryonFileName,"%s.cpu%"TASK_TAG_FORMAT""ISYM, base_name,ProcessorNumber);


  // ***** Create Group For This Grid *****

  sprintf(GroupName,"Grid%"GROUP_TAG_FORMAT""ISYM, ID);

  if (io_log) fprintf(log_fptr, "Calling H5Gcreate with Name %s\n", GroupName);
  group_id = H5Gcreate(level_group_id, GroupName, 0);
  if (io_log) fprintf(log_fptr, "H5Gcreate: group_id = %"ISYM"\n", (int) group_id);

  // ***** Write Grid Attributes *****

  // Task
  HDF5_WriteAttribute(group_id, "Task", ProcessorNumber, log_fptr);

  // GridRank
  HDF5_WriteAttribute(group_id, "GridRank", GridRank, log_fptr);

  // Time
  HDF5_WriteAttribute(group_id, "Time", Time, log_fptr);

  // OldTime
  HDF5_WriteAttribute(group_id, "OldTime", OldTime, log_fptr);

  // SubgridsAreStatic
  HDF5_WriteAttribute(group_id, "SubgridsAreStatic", SubgridsAreStatic, log_fptr);

  // NumberOfBaryonFields
  HDF5_WriteAttribute(group_id, "NumberOfBaryonFields", NumberOfBaryonFields, log_fptr);

  if (NumberOfBaryonFields > 0) {
    // FieldType
    HDF5_WriteAttribute(group_id, "FieldType", FieldType, NumberOfBaryonFields, log_fptr);

    // BaryonFileName
    HDF5_WriteAttribute(group_id, "BaryonFileName", BaryonFileName, log_fptr);

    // CourantSafetyNumber
    HDF5_WriteAttribute(group_id, "CourantSafetyNumber", CourantSafetyNumber, log_fptr);
   
    // PPMFlatteningParameter
    HDF5_WriteAttribute(group_id, "PPMFlatteningParameter", PPMFlatteningParameter, log_fptr);

    // PPMDiffusionParameter
    HDF5_WriteAttribute(group_id, "PPMDiffusionParameter", PPMDiffusionParameter, log_fptr);

    // PPMSteepeningParameter
    HDF5_WriteAttribute(group_id, "PPMSteepeningParameter", PPMSteepeningParameter, log_fptr);
  }
  
  if (NumberOfParticles > 0) {
    // ParticleFileName
    HDF5_WriteAttribute(group_id, "ParticleFileName", ParticleFileName, log_fptr);
  }

  if (SelfGravity) {
    // GravityBoundaryType
    HDF5_WriteAttribute(group_id, "GravityBoundaryType", GravityBoundaryType, log_fptr);
  }
  
  // NumberOfDaughterGrids
  HDF5_WriteAttribute(group_id, "NumberOfDaughterGrids", NumberOfDaughterGrids, log_fptr);

  // NextGridThisLevelID
  HDF5_WriteAttribute(group_id, "NextGridThisLevelID", NextGridThisLevelID, log_fptr);
  
  // NextGridNextLevelID
  HDF5_WriteAttribute(group_id, "NextGridNextLevelID", NextGridNextLevelID, log_fptr);




  // ***** Write Grid Datasets *****

  // GridDimension
  HDF5_WriteDataset(group_id, "GridDimension", GridDimension, GridRank, log_fptr);

  // GridStartIndex
  HDF5_WriteDataset(group_id, "GridStartIndex", GridStartIndex, GridRank, log_fptr);

  // GridEndIndex
  HDF5_WriteDataset(group_id, "GridEndIndex", GridEndIndex, GridRank, log_fptr);

  // GridLeftEdge
  HDF5_WriteDataset(group_id, "GridLeftEdge", GridLeftEdge, GridRank, log_fptr);

  // GridRightEdge
  HDF5_WriteDataset(group_id, "GridRightEdge", GridRightEdge, GridRank, log_fptr);

  // GridGlobalPosition

  // This is LeftGridEdge[] expressed in integer indices of this
  // level, i.e. running from 0 to RootGridDimension[] *
  // RefinementFactors[]**level - 1. This may be useful for
  // re-calculating positions in long double precision (which is not
  // universally supported by HDF5) at runtime.

  for (int i=0; i<GridRank; i++)
    GridGlobalPosition[i] = (int) ( (GridLeftEdge[i] - DomainLeftEdge[i]) / CellWidth[i][0] );

  HDF5_WriteDataset(group_id, "GridGlobalPosition", GridGlobalPosition, GridRank, log_fptr);

  // NumberOfParticles
  HDF5_WriteDataset(group_id, "NumberOfParticles", NumberOfParticles, log_fptr);


  // ***** External Link to the grid data in the external .cpuXXXX file *****
  // (only supported under HDF5 1.8+)

#ifdef WITH_HDF5_LINKS

#if HDF5_VERSION_GE(1,8,0)
  sprintf(TargetName,"Grid%"GROUP_TAG_FORMAT""ISYM, ID);
  sprintf(LinkName,"GridData");
  if (io_log) fprintf(log_fptr,"H5Lcreate_external: %s:%s -> %s\n", BaryonFileName, TargetName, LinkName);
  h5_status = H5Lcreate_external(BaryonFileName, TargetName, group_id, LinkName, H5P_DEFAULT, H5P_DEFAULT);


  if (io_log) fprintf(log_fptr, "H5Lcreate_external: status = %"ISYM"\n", (int) h5_status);
#endif // HDF5_VERSION_GE(1,8,0)

#endif

  // ***** Links to daughter grids *****
  if (NumberOfDaughterGrids > 0) {

    // create DaughterGrids subgroup
    sprintf(GroupName,"DaughterGrids");
    if (io_log) fprintf(log_fptr, "Calling H5Gcreate with Name %s\n", GroupName);
    subgroup_id = H5Gcreate(group_id, GroupName, 0);
    if (io_log) fprintf(log_fptr, "H5Gcreate: subgroup_id = %"ISYM"\n", (int) subgroup_id);

    // write daughter grid id array
    HDF5_WriteAttribute(subgroup_id, "DaughterGridIDs", DaughterGridIDs, NumberOfDaughterGrids, log_fptr);
    
#ifdef WITH_HDF5_LINKS    
    // create soft links to daughter grids
    for(int i=0;i<NumberOfDaughterGrids;i++) {
      
      sprintf(LinkName,"DaughterGrid%"GRID_TAG_FORMAT""ISYM,i);
      sprintf(TargetName,"/Level%"ISYM"/Grid%"GROUP_TAG_FORMAT""ISYM,level+1,DaughterGridIDs[i]);

#if HDF5_VERSION_GE(1,8,0)
      if (io_log) fprintf(log_fptr,"H5Lcreate_soft: %s -> %s\n", TargetName, LinkName);
      h5_status = H5Lcreate_soft(TargetName, subgroup_id, LinkName, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Lcreate_soft: status = %"ISYM"\n", (int) h5_status);
#else
      if (io_log) fprintf(log_fptr,"H5Glink: %s -> %s\n", TargetName, LinkName);
      h5_status = H5Glink(subgroup_id, H5G_LINK_SOFT, TargetName, LinkName);
      if (io_log) fprintf(log_fptr, "H5Glink: status = %"ISYM"\n", (int) h5_status);
#endif // HDF5_VERSION_GE(1,8,0)

    }
#endif

    h5_status = H5Gclose(subgroup_id);
    if (io_log) fprintf(log_fptr, "H5Gclose: status = %"ISYM"\n", (int) h5_status);
    
  }

#ifdef WITH_HDF5_LINKS
  // ***** Links to parent grids *****
  if (level > 0) {

    // create ParentGrids subgroup
    sprintf(GroupName,"ParentGrids");
    if (io_log) fprintf(log_fptr, "Calling H5Gcreate with Name %s\n", GroupName);
    subgroup_id = H5Gcreate(group_id, GroupName, 0);
    if (io_log) fprintf(log_fptr, "H5Gcreate: subgroup_id = %"ISYM"\n", (int) subgroup_id);

    for(int i=0;i<level;i++) {

      sprintf(LinkName,"ParentGrid_Level%"ISYM,i);
      sprintf(TargetName,"/Level%"ISYM"/Grid%"GROUP_TAG_FORMAT""ISYM,i,ParentGridIDs[level-1-i]);

#if HDF5_VERSION_GE(1,8,0)

      if (io_log) fprintf(log_fptr,"H5Lcreate_soft: %s -> %s\n", TargetName, LinkName);
      h5_status = H5Lcreate_soft(TargetName, subgroup_id, LinkName, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Lcreate_soft: status = %"ISYM"\n", (int) h5_status);
#else
      if (io_log) fprintf(log_fptr,"H5Glink: %s -> %s\n", TargetName, LinkName);
      h5_status = H5Glink(subgroup_id, H5G_LINK_SOFT, TargetName, LinkName);
      if (io_log) fprintf(log_fptr, "H5Glink: status = %"ISYM"\n", (int) h5_status);
#endif // HDF5_VERSION_GE(1,8,0)

    }

    h5_status = H5Gclose(subgroup_id);
    if (io_log) fprintf(log_fptr, "H5Gclose: status = %"ISYM"\n", (int) h5_status);
    
  }
#endif


  // ***** Close this grid *****
  h5_status = H5Gclose(group_id);
  if (io_log) fprintf(log_fptr, "H5Gclose: status = %"ISYM"\n", (int) h5_status);

    
  return SUCCESS;
}




// HDF5 utility routines (to write attributes and datasets)

// int
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eint32 Attribute, FILE *log_fptr) {

  hid_t dspace_id, attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dspace_id = H5Screate(H5S_SCALAR);
  if (io_log) fprintf(log_fptr, "H5Screate: dspace_id = %"ISYM"\n", (int) dspace_id);

  attr_id = H5Acreate(group_id, AttributeName, HDF5_FILE_INT, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Acreate: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Awrite(attr_id,  HDF5_INT, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Awrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eint64 Attribute, FILE *log_fptr) {

  hid_t dspace_id, attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dspace_id = H5Screate(H5S_SCALAR);
  if (io_log) fprintf(log_fptr, "H5Screate: dspace_id = %"ISYM"\n", (int) dspace_id);

  attr_id = H5Acreate(group_id, AttributeName, HDF5_FILE_I8, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Acreate: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Awrite(attr_id,  HDF5_I8, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Awrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}


// 32-bit float (Eflt32)
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eflt32 Attribute, FILE *log_fptr) {

  hid_t dspace_id, attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dspace_id = H5Screate(H5S_SCALAR);
  if (io_log) fprintf(log_fptr, "H5Screate: dspace_id = %"ISYM"\n", (int) dspace_id);

  attr_id = H5Acreate(group_id, AttributeName, HDF5_R4, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Acreate: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Awrite(attr_id,  HDF5_R4, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Awrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}


// 64-bit float (Eflt64)
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eflt64 Attribute, FILE *log_fptr) {

  hid_t dspace_id, attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dspace_id = H5Screate(H5S_SCALAR);
  if (io_log) fprintf(log_fptr, "H5Screate: dspace_id = %"ISYM"\n", (int) dspace_id);

  attr_id = H5Acreate(group_id, AttributeName, HDF5_R8, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Acreate: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Awrite(attr_id,  HDF5_R8, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Awrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// 128-bit float (Eflt128)
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eflt128 Attribute, FILE *log_fptr) {

  hid_t dspace_id, attr_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dspace_id = H5Screate(H5S_SCALAR);
  if (io_log) fprintf(log_fptr, "H5Screate: dspace_id = %"ISYM"\n", (int) dspace_id);

  attr_id = H5Acreate(group_id, AttributeName, HDF5_R16, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Acreate: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Awrite(attr_id,  HDF5_R16, &Attribute);
  if (io_log) fprintf(log_fptr, "H5Awrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}


// int vector
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, int *Attribute, int NumberOfElements, FILE *log_fptr) {

  hid_t dspace_id, attr_id;

  hsize_t OutDims[MAX_DIMENSION];

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  OutDims[0] = NumberOfElements;
  dspace_id = H5Screate_simple(1, OutDims, NULL);
  if (io_log) fprintf(log_fptr, "H5Screate_simple: dspace_id = %"ISYM"\n", (int) dspace_id);
  
  attr_id = H5Acreate(group_id, AttributeName, HDF5_FILE_INT, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Acreate: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Awrite(attr_id,  HDF5_INT, Attribute);
  if (io_log) fprintf(log_fptr, "H5Awrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// string (char vector)
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, char *Attribute, FILE *log_fptr) {

  hid_t dspace_id, attr_id, string_type_id;

  hsize_t string_size;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dspace_id = H5Screate(H5S_SCALAR);   
  if (io_log) fprintf(log_fptr, "H5Screate: dspace_id = %"ISYM"\n", (int) dspace_id);

  string_type_id = H5Tcopy(H5T_C_S1);
  string_size = (hsize_t) strlen(Attribute);
  h5_status = H5Tset_size(string_type_id, string_size);

  attr_id = H5Acreate(group_id, AttributeName, string_type_id, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Acreate: attr_id = %"ISYM"\n", (int) attr_id);

  h5_status = H5Awrite(attr_id,  string_type_id, Attribute);
  if (io_log) fprintf(log_fptr, "H5Awrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Aclose(attr_id);
  if (io_log) fprintf(log_fptr, "H5Aclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Tclose(string_type_id);
  if (io_log) fprintf(log_fptr, "H5Tclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// int
int HDF5_WriteDataset(hid_t group_id, const char *DatasetName, int Dataset, FILE *log_fptr) {

  hid_t dspace_id, dset_id;

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  dspace_id = H5Screate(H5S_SCALAR);
  if (io_log) fprintf(log_fptr, "H5Screate: dspace_id = %"ISYM"\n", (int) dspace_id);

  dset_id = H5Dcreate(group_id, DatasetName, HDF5_FILE_INT, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Dcreate: dset_id = %"ISYM"\n", (int) dset_id);

  h5_status = H5Dwrite(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) (&Dataset));
  if (io_log) fprintf(log_fptr, "H5Dwrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Dclose(dset_id);
  if (io_log) fprintf(log_fptr, "H5Dclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// int vector
int HDF5_WriteDataset(hid_t group_id, const char *DatasetName, int *Dataset, int NumberOfElements, FILE *log_fptr) {

  hid_t dspace_id, dset_id;

  hsize_t OutDims[MAX_DIMENSION];

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  OutDims[0] = NumberOfElements;
  dspace_id = H5Screate_simple(1, OutDims, NULL);
  if (io_log) fprintf(log_fptr, "H5Screate_simple: dspace_id = %"ISYM"\n", (int) dspace_id);

  dset_id = H5Dcreate(group_id, DatasetName, HDF5_FILE_INT, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Dcreate: dset_id = %"ISYM"\n", (int) dset_id);

  h5_status = H5Dwrite(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) Dataset);
  if (io_log) fprintf(log_fptr, "H5Dwrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Dclose(dset_id);
  if (io_log) fprintf(log_fptr, "H5Dclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

// FLOAT vector
int HDF5_WriteDataset(hid_t group_id, const char *DatasetName, FLOAT *Dataset, int NumberOfElements, FILE *log_fptr) {

  hid_t dspace_id, dset_id;

  hsize_t OutDims[MAX_DIMENSION];

  herr_t h5_status;

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  OutDims[0] = NumberOfElements;
  dspace_id = H5Screate_simple(1, OutDims, NULL);
  if (io_log) fprintf(log_fptr, "H5Screate_simple: dspace_id = %"ISYM"\n", (int) dspace_id);

  dset_id = H5Dcreate(group_id, DatasetName, HDF5_FILE_PREC, dspace_id, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Dcreate: dset_id = %"ISYM"\n", (int) dset_id);

  h5_status = H5Dwrite(dset_id, HDF5_PREC, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) Dataset);
  if (io_log) fprintf(log_fptr, "H5Dwrite: status = %"ISYM"\n", (int) h5_status);

  h5_status = H5Dclose(dset_id);
  if (io_log) fprintf(log_fptr, "H5Dclose: status = %"ISYM"\n", (int) h5_status);
  
  h5_status = H5Sclose(dspace_id);
  if (io_log) fprintf(log_fptr, "H5Sclose: status = %"ISYM"\n", (int) h5_status);

  return SUCCESS;
}

