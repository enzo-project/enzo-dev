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
		ExternalBoundary *Exterior, float *Initialdt)
 
{
 
  /* declarations */
 
  char pid[MAX_TASK_TAG_SIZE];
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
 
  /* Read Data Hierarchy. */
 
  if ((fptr = fopen(hierarchyname, "r")) == NULL) {
    ENZO_VFAIL("Error opening hierarchy file %s.\n", hierarchyname)
  }
  GridID = 1;
  if (ReadDataHierarchy(fptr, TopGrid, GridID, NULL) == FAIL) {
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
 
  if (ReadStarParticleData(fptr) == FAIL) {
        ENZO_FAIL("Error in ReadStarParticleData.");
  }
 
  /* Create radiation name and read radiation data. */
 
  if ((RadiationFieldType >= 10 && RadiationFieldType <= 11) || 
      RadiationData.RadiationShield == TRUE) {
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
 
  fclose(fptr);

  /* If we added new particle attributes, unset flag so we don't carry
     this parameter to later data. */

  AddParticleAttributes = FALSE;

  return SUCCESS;
}
