/***********************************************************************
/
/  WRITE OUT THE HIERARCHY INTO AN HDF5 FILE
/
/  written by: Michael Kuhlen
/  date:       October, 2010
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

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

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int SetGlobalGridID(int &GlobalID, HierarchyEntry *Grid);

// the following HDF5 helper routines are defined in
// Grid_WriteHierarchyInformationHDF5.C
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, int Attribute, FILE *log_fptr);
#ifdef SMALL_INTS
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, Eint64 Attribute, FILE *log_fptr);
#endif
int HDF5_WriteAttribute(hid_t group_id, const char *AttributeName, FLOAT Attribute, FILE *log_fptr);
int HDF5_WriteDataset(hid_t group_id, const char *DatasetName, int *Dataset, int NumberOfElements, FILE *log_fptr);

// #define IO_LOG


int WriteHDF5HierarchyFile(char *base_name, HierarchyEntry *TopGrid, TopGridData MetaData, LevelHierarchyEntry *LevelArray[]) {
  char FileName[MAX_LINE_LENGTH];
  char GroupName[MAX_LINE_LENGTH];
  char DatasetName[MAX_LINE_LENGTH];
  char AttributeName[MAX_LINE_LENGTH];

  int level, i;
  int FinestLevel=0;
  int *ParentGridIDs, *DaughterGridIDs;
  int *LevelLookupTable;
  int *OriginalIDs;

  FILE *log_fptr;

  LevelHierarchyEntry *LTemp;
  HierarchyEntry *HTemp, *HTemp2;

  hid_t       file_id, group_id, dset_id, attr_id;
  hid_t       dspace_id;
  hid_t       int_file_type_id, int_mem_type_id;  

  hsize_t     OutDims[MAX_DIMENSION];

  herr_t      h5_status;
  herr_t      h5_error = -1;  

  int io_log = 0;
#ifdef IO_LOG
  io_log = 1;
#endif

  int jj = sizeof(int); 
  switch(jj)  {
    
  case 4:
    int_mem_type_id = HDF5_I4;
    int_file_type_id = HDF5_FILE_I4;
    break;
    
  case 8:
    int_mem_type_id = HDF5_I8;
    int_file_type_id = HDF5_FILE_I8;
    break;
    
  default:
    printf("INCORRECT int DEFINITION\n");    
  }

  if ( MyProcessorNumber != ROOT_PROCESSOR )
    return SUCCESS;

  sprintf(FileName,"%s.hierarchy.hdf5",base_name);

  char logname[MAX_LINE_LENGTH];
  sprintf(logname,"%s.log",FileName);
  if (io_log) log_fptr = fopen(logname, "w");

  // open the file
  if (io_log) fprintf(log_fptr, "Calling H5Fcreate with Name = %s\n", FileName);
  file_id = H5Fcreate(FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);

  
  // Calculate CurrentRedshift and add as attribute (if Cosmology)
  FLOAT CurrentRedshift, a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
  CurrentRedshift = (1 + InitialRedshift)/a - 1;

  HDF5_WriteAttribute(file_id, "Redshift", CurrentRedshift, log_fptr);


  // Add NumberOfProcessors attribute to the root group
  HDF5_WriteAttribute(file_id, "NumberOfProcessors", NumberOfProcessors, log_fptr);

  if (StarParticleCreation) {
    // Add NumberOfStarParticles attribute to the root group
    HDF5_WriteAttribute(file_id, "NumberOfStarParticles", NumberOfStarParticles, log_fptr);
    
    // Add NumberOfOtherParticles attribute to the root group
    HDF5_WriteAttribute(file_id, "NumberOfOtherParticles", NumberOfOtherParticles, log_fptr);
  }
  
  // find finest level with grids
  for(level = MAX_DEPTH_OF_HIERARCHY-1; LevelArray[level]==NULL; level--)
    FinestLevel = level;
  FinestLevel--;
  if (io_log) fprintf(log_fptr, "FinestLevel = %"ISYM"\n",FinestLevel);


  int NumberOfGrids[FinestLevel+1];
  int TotalNumberOfGrids = 0;
  for(level=0; level<=FinestLevel; level++) {
    NumberOfGrids[level] = 0;
    for(LTemp=LevelArray[level]; LTemp; LTemp=LTemp->NextGridThisLevel, NumberOfGrids[level]++);
    TotalNumberOfGrids += NumberOfGrids[level];
  }


  // Add TotalNumberOfGrids attribute to the root group
  HDF5_WriteAttribute(file_id, "TotalNumberOfGrids", TotalNumberOfGrids, log_fptr);


  // We need global GridID's (i.e. running from 0 to
  // TotalNumberOfGrids), but what's currently stored in Grid::ID is
  // the grid's ID on this level (i.e. running from 0 to
  // NumberOfGridsThisLevel). Here I temporarily change Grid:ID to
  // refer to the global GridID, but this is restored to the level
  // GridID at the end of this routine.

  // First, save original grid ids.
  OriginalIDs = new int[TotalNumberOfGrids];
  for(level=0, i=0; level<=FinestLevel; level++)
    for(LTemp=LevelArray[level]; LTemp; LTemp=LTemp->NextGridThisLevel)
      OriginalIDs[i++] = LTemp->GridData->GetGridID();
  
  // Now, set grid ids to the global ones.
  int GlobalID=1;
  SetGlobalGridID(GlobalID, TopGrid);


  // create all level groups  
  for(level=0; level<=FinestLevel; level++) {
    sprintf(GroupName,"/Level%"ISYM,level);
    
    if (io_log) fprintf(log_fptr, "Calling H5Gcreate with Name %s\n", GroupName);
    group_id = H5Gcreate(file_id, GroupName, 0);
    if (io_log) fprintf(log_fptr, "H5Gcreate: %"ISYM"\n", (int) group_id);
    

    // create all grid groups and and write its datasets
    for (LTemp = LevelArray[level]; LTemp; LTemp=LTemp->NextGridThisLevel) {

      // find ParentGridIs
      if (level > 0) {
	ParentGridIDs = new int[level];
	
	int counter=0;
	for (HTemp = LTemp->GridHierarchyEntry->ParentGrid; HTemp; HTemp = HTemp->ParentGrid) {
	  ParentGridIDs[counter++] = HTemp->GridData->GetGridID();
	}
      }

      // find DaughterGridIDs
      int NumberOfDaughterGrids=0;
      for (HTemp2 = LTemp->GridHierarchyEntry->NextGridNextLevel; HTemp2; HTemp2 = HTemp2->NextGridThisLevel)
	NumberOfDaughterGrids++;
      
      if (NumberOfDaughterGrids>0) {
	DaughterGridIDs = new int[NumberOfDaughterGrids];
	
	int counter=0;
	for (HTemp2 = LTemp->GridHierarchyEntry->NextGridNextLevel; HTemp2; HTemp2 = HTemp2->NextGridThisLevel)
	  DaughterGridIDs[counter++] = HTemp2->GridData->GetGridID();
      }
      
      
      int NextGridThisLevelID = 0;
      if (LTemp->GridHierarchyEntry->NextGridThisLevel) 
	NextGridThisLevelID = LTemp->GridHierarchyEntry->NextGridThisLevel->GridData->GetGridID();

      int NextGridNextLevelID = 0;
      if (LTemp->GridHierarchyEntry->NextGridNextLevel) 
	NextGridNextLevelID = LTemp->GridHierarchyEntry->NextGridNextLevel->GridData->GetGridID();

      LTemp->GridData->WriteHierarchyInformationHDF5(base_name, group_id, level, ParentGridIDs, NumberOfDaughterGrids, DaughterGridIDs, NextGridThisLevelID, NextGridNextLevelID, log_fptr);
      
      // fprintf(stderr,"level=%"ISYM" grid=%"ISYM"\n",level,LTemp->GridData->GetGridID());
      // fprintf(stderr,"  Parents:");
      // if(level==0) fprintf(stderr," none\n");
      // else {
      // 	for(int i=0;i<level;i++)
      // 	  fprintf(stderr," %"ISYM, ParentGridIDs[i]);
      // 	fprintf(stderr,"\n");	
      // }
      // fprintf(stderr,"  Daughters: %"ISYM" (",NumberOfDaughterGrids);
      // for(i=0;i<NumberOfDaughterGrids;i++)
      // 	fprintf(stderr," %"ISYM, DaughterGridIDs[i]);
      // fprintf(stderr,")\n");	
      // fflush(stderr);


      if (level > 0) delete [] ParentGridIDs;
      if (NumberOfDaughterGrids > 0) delete [] DaughterGridIDs;
	
    } // loop over grids this level

    //    fprintf(stderr,"  Number of grids this level = %"ISYM"\n",NumberOfGrids[level]);

    // Add NumberOfGrids Attribute
    HDF5_WriteAttribute(group_id, "NumberOfGrids", NumberOfGrids[level], log_fptr);

    // Close this group
    h5_status = H5Gclose(group_id);
    if (io_log) fprintf(log_fptr, "H5Gclose: %"ISYM"\n", (int) h5_status);
    
  } // loop over levels
  //  fprintf(stderr,"Total number of grids = %"ISYM"\n",TotalNumberOfGrids);


  // Create a LevelLookupTable, which is just an int array containing
  // for each grid the level that it's on.

  // are GridIDs always sequential, meaning the largest GridID is
  // equal to TotalNumberOfGrids? (I hope so...)
  LevelLookupTable = new int[TotalNumberOfGrids];
  for(level=0; level<=FinestLevel; level++) {
    for (LTemp = LevelArray[level]; LTemp; LTemp=LTemp->NextGridThisLevel) {
      LevelLookupTable[LTemp->GridData->GetGridID()-1] = level;
      //      fprintf(stderr,"LevelLookup[%d] = %d\n",LTemp->GridData->GetGridID()-1,level);
    }
  }

  // for(i=0;i<TotalNumberOfGrids;i++)
  //   fprintf(stderr,"%d %d\n",i,LevelLookupTable[i]);
  // fflush(stderr);

  // Write LevelLookupTable Dataset
  HDF5_WriteDataset(file_id, "LevelLookupTable", LevelLookupTable, TotalNumberOfGrids, log_fptr);


  // Close the file
  h5_status = H5Fclose(file_id);

  if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", (int) h5_status);


  if (io_log) fclose(log_fptr);

  // Restore Grid::ID to the original grid id.
  for(level=0, i=0; level<=FinestLevel; level++)
    for(LTemp=LevelArray[level]; LTemp; LTemp=LTemp->NextGridThisLevel)
      LTemp->GridData->SetGridID(OriginalIDs[i++]);

  delete [] LevelLookupTable;
  delete [] OriginalIDs;

  return SUCCESS;
}


int SetGlobalGridID(int &GlobalID, HierarchyEntry *Grid) {
  Grid->GridData->SetGridID(GlobalID++);

  if(Grid->NextGridThisLevel)
    SetGlobalGridID(GlobalID, Grid->NextGridThisLevel);

  if(Grid->NextGridNextLevel)
    SetGlobalGridID(GlobalID, Grid->NextGridNextLevel);
 
  return SUCCESS;
}
