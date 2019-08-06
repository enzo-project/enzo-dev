/***********************************************************************
/
/  Convert star objects to actve particles
/
/  written by: John Regan
/  date:       May, 2018
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#include "preincludes.h"

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
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif
#include "CosmologyParameters.h"

#define DEBUG 0
void CollectParticleTypes(char **active_particle_types, int numparticles);

int ActiveParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
               int NumberOfGrids, LevelHierarchyEntry *LevelArray[],
               int level, int NumberOfNewActiveParticles[]);
#ifdef TRANSFER
int RadiativeTransferInitialize(char *ParameterFile, 
				HierarchyEntry &TopGrid, 
				TopGridData &MetaData,
				ExternalBoundary &Exterior, 
				ImplicitProblemABC* &ImplicitSolver,
				LevelHierarchyEntry *LevelArray[]);
#endif

int Group_WriteAllData(char *basename, int filenumber,
		       HierarchyEntry *TopGrid, TopGridData &MetaData,
		       ExternalBoundary *Exterior, 
#ifdef TRANSFER
		       ImplicitProblemABC *ImplicitSolver,
#endif
		       FLOAT WriteTime = -1, 
		       int CheckpointDump = FALSE);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

#ifdef TRANSFER
int RadiativeTransferReadParameters(FILE *fptr);
#endif
void my_exit(int status);
int ConvertParticles2ActiveParticles(char *ParameterFile,
			  LevelHierarchyEntry *LevelArray[], 
			  HierarchyEntry *TopGrid,
			  TopGridData &MetaData,
			  ExternalBoundary &Exterior
#ifdef TRANSFER
			  , ImplicitProblemABC *ImplicitSolver
#endif
			  )
{

  
  const float When = 0.5;
  LevelHierarchyEntry *Temp;
  HierarchyEntry **Grids;
  int numtypes = 0;
  
  int maxlevel = 0;
  while(LevelArray[maxlevel] != NULL)
    maxlevel++;
  /* Determine the parameter name prefix */

  char *cptr;
  int DumpNumber;
  char DumpName[MAX_LINE_LENGTH];
  char NumberString[MAX_LINE_LENGTH];
  if ( (cptr = strstr(ParameterFile, MetaData.DataDumpName)) ) {
    strcpy(DumpName, MetaData.DataDumpName);
  }
  else if ( (cptr = strstr(ParameterFile, MetaData.RestartDumpName)) ) {
    strcpy(DumpName, MetaData.RestartDumpName);
  }
  else if ( (cptr = strstr(ParameterFile, MetaData.RedshiftDumpName)) ) {
    strcpy(DumpName, MetaData.RedshiftDumpName);
  }
  else
    ENZO_FAIL("Cannot determine output type.");

  /* Extract output number */
  
  strncpy(NumberString, ParameterFile+2, 4);
  NumberString[4] = '\0';
  DumpNumber = atoi(NumberString);

  strcat(MetaData.GlobalDir, "/ActiveParticles");
  int Masterarray[MAX_ACTIVE_PARTICLE_TYPES], RMasterarray[MAX_ACTIVE_PARTICLE_TYPES];
  
  char **active_particle_types;
  active_particle_types = new char*[MAX_ACTIVE_PARTICLE_TYPES];
  for(int i=0; i<MAX_ACTIVE_PARTICLE_TYPES; i++) {
    active_particle_types[i] = new char[MAX_LINE_LENGTH];
    strcpy(active_particle_types[i], "");
    Masterarray[i] = 0;
    RMasterarray[i]=0;
  }
  int active_particles = 0, global_active_particles = 0;
  if(MyProcessorNumber == ROOT_PROCESSOR)
    printf("%s: You have chosen to read in a file and convert all the particles " \
	   "found in the output to active particles.\n", __FUNCTION__);


  
  /* Initialize the radiative transfer */

#ifdef TRANSFER
  if (RadiativeTransferInitialize(ParameterFile, *TopGrid, MetaData, Exterior, 
				  ImplicitSolver, LevelArray) == FAIL) {
    fprintf(stderr, "Error in RadiativeTransferInitialize.\n");
    my_exit(EXIT_FAILURE);
  }
#endif
 
  /* 
   * Determine what type of star objects exist and populate the 
   * ActiveParticleType array
   */
  for (int level = 0; level < maxlevel; level++) {
    int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
    for (int grid1 = 0; grid1 < NumberOfGrids; grid1++) {
       numtypes += Grids[grid1]->GridData->DetermineActiveParticleTypes(active_particle_types);
    }
  }
  for(int i = 0; i < MAX_ACTIVE_PARTICLE_TYPES; i++)
    {
      if(strlen(active_particle_types[i]) > 0) {
	active_particles++;
	Masterarray[i] = 1;
      }
    }
#ifdef USE_MPI
  MPI_Allreduce(&numtypes, &global_active_particles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(Masterarray, RMasterarray, MAX_ACTIVE_PARTICLE_TYPES, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if(MyProcessorNumber == ROOT_PROCESSOR)
    printf("%s: Number of particles found = %d\n", __FUNCTION__, global_active_particles);
#endif

#if DEBUG
  /* Now lets collect up the active particle types */
  //CollectParticleTypes(active_particle_types, global_active_particles);
  //if(MyProcessorNumber == ROOT_PROCESSOR)    {
  // for(int i = 0; i < global_active_particles; i++)
  //   fprintf(stdout, "P%d: particle type[%d] =  %s\n", MyProcessorNumber, 
  //	      i, active_particle_types[i]);
  //}
#endif
  /* Copy the active particle names into the arrays mapped from the particle numbers */
  for(int i = 0; i < MAX_ACTIVE_PARTICLE_TYPES; i++) {
    if(RMasterarray[i]) {
      if(i == PARTICLE_TYPE_SINGLE_STAR)
	strcpy(active_particle_types[i], "PopIII");
      if(i == PARTICLE_TYPE_CLUSTER)
	strcpy(active_particle_types[i], "CenOstriker");
    }
  }
 
  // Enable the active particles that were selected.
  for (int i = 0;i < MAX_ACTIVE_PARTICLE_TYPES; i++) {
    if(strlen(active_particle_types[i]) > 0)
      EnableActiveParticleType(active_particle_types[i]);
  }
 
  int TotalNumberOfNewParticles = 0, GlobalTotalNumberOfNewParticles=0;

  for (int level = 0; level < maxlevel; level++) {
    int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
    int *NumberOfNewActiveParticles = new int[NumberOfGrids]();
    for (int grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->ActiveParticleHandler_Convert(Grids[grid1]->NextGridThisLevel, level,
							    grid1, NumberOfNewActiveParticles[grid1]);
     
      TotalNumberOfNewParticles += NumberOfNewActiveParticles[grid1];
      
    }
    ActiveParticleFinalize(Grids, &MetaData, NumberOfGrids, LevelArray,
                          level, NumberOfNewActiveParticles);

  }

#ifdef USE_MPI
  MPI_Allreduce(&TotalNumberOfNewParticles, &GlobalTotalNumberOfNewParticles, 1, MPI_INT, MPI_SUM, 
		MPI_COMM_WORLD);
  if(GlobalTotalNumberOfNewParticles && MyProcessorNumber == ROOT_PROCESSOR)
    printf("%s: TotalNumberOfNewParticles = %d\n", __FUNCTION__, GlobalTotalNumberOfNewParticles);
#endif
  /* 
   * The new active particle types have now been created and the star objects destroyed.
   * Write the particles and active particles back to disk
   */
   Group_WriteAllData(DumpName, DumpNumber, TopGrid, MetaData, &Exterior
#ifdef TRANSFER
		       , ImplicitSolver
#endif
               );
  
  return SUCCESS;
}

void CollectParticleTypes(char **active_particle_types, int global_active_particles)
{
#ifdef USE_MPI
  Eint32 my_gsize = 0;
  char **ap_types;
  char **ap_types2;
  MPI_Comm_size(MPI_COMM_WORLD, &my_gsize);
  ap_types = new char*[global_active_particles];
  ap_types2 = new char*[global_active_particles*my_gsize];
  int j = 0;
  for(int i=0; i<global_active_particles; i++) {
    ap_types[j] = new char[MAX_LINE_LENGTH];
    strcpy(ap_types[j], "");
    j = j + 1;
  }
 
  j = 0;
  for(int i=0; i<global_active_particles*my_gsize; i++) {
    ap_types2[j] = new char[MAX_LINE_LENGTH];
    strcpy(ap_types2[j], "");
    j = j + 1;
  }
  
  j = 0;
  for(int i = 0; i <MAX_ACTIVE_PARTICLE_TYPES; i++)
    if(strlen(active_particle_types[i]) > 0) {
      strcpy(ap_types[j], active_particle_types[i]);
      j = j + 1;
    }
 
  for(int i=0; i<global_active_particles*my_gsize; i++) {
    if(MyProcessorNumber == 75) {
      printf("P%d: ap_types[%d] = %s\n", MyProcessorNumber, i, ap_types[i]);
    }
  }
  return;
  for(int i=0; i<global_active_particles; i++) {
    int offset = global_active_particles*MyProcessorNumber;
    MPI_Allgather(ap_types[i], MAX_LINE_LENGTH, MPI_CHAR,
		  ap_types2[i], MAX_LINE_LENGTH, MPI_CHAR, MPI_COMM_WORLD);
  }
  if(MyProcessorNumber == ROOT_PROCESSOR) {
    int k = 0;
    for(int i=0; i<global_active_particles*my_gsize; i++) {
      for(int j=0; j < i; j++) {
	if(strcmp(active_particle_types[j], ap_types2[i]) == 0)
	  continue;
      }
      printf("Populating array with %s\n", ap_types2[i]);
      strcpy(active_particle_types[k++], ap_types2[i]);
    // fprintf(stdout, "P%d: Enabling particle type[%d] =  %s\n", MyProcessorNumber, i, active_particle_types[i]); fflush(stdout);
    //MPI_Barrier(MPI_COMM_WORLD);
    }
  }
 
  return;
#endif
}
