/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A TURBULENCE SIMULATION)
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/  modified1:  Robert Harkness
/  date:       April 2008
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */
 
#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

 
//#include "performance.h"


 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "flowdefs.h"
#include "error.h"
void my_exit(int status);
 
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif
 
// HDF5 function prototypes
 

 
// function prototypes
 
int ReadFile(char *name, int Rank, int Dim[], int StartIndex[],
                  int EndIndex[], int BufferOffset[], float *buffer,
                  inits_type **tempbuffer, int Part, int Npart);
 
int grid::TurbulenceSimulationInitializeGrid(
                          float TurbulenceSimulationInitialDensity,
                          float TurbulenceSimulationInitialTemperature,
                          char *TurbulenceSimulationDensityName,
                          char *TurbulenceSimulationTotalEnergyName,
                          char *TurbulenceSimulationGasEnergyName,
                          char *TurbulenceSimulationVelocityNames[],
                          char *TurbulenceSimulationRandomForcingNames[],
                          int   TurbulenceSimulationSubgridsAreStatic,
                          int   TotalRefinement)
{
 
  /* declarations */
 
  int idim, dim, i, j, vel;
  int DeNum;
 
  int ExtraField[2];
 
  inits_type *tempbuffer = NULL;
 
  FILE *log_fptr;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
  /*
  if ( NumberOfProcessors > 64 )
    if (ParallelRootGridIO != TRUE) {
      ENZO_FAIL("ParallelRootGridIO MUST be set for > 64 cpus!\n");
    }
  */
 
  char pid[MAX_TASK_TAG_SIZE];
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  char *logname = new char[MAX_NAME_LENGTH];
  strcpy(logname, "TSlog.");
  strcat(logname,pid);
 
  if (io_log) {
    log_fptr = fopen(logname, "a");
    fprintf(log_fptr, "\n");
    fprintf(log_fptr, "TSIG ParallelRootGridIO = %"ISYM"\n", ParallelRootGridIO);
    fprintf(log_fptr, "Processor %"ISYM", Target processor %"ISYM"\n",
	    MyProcessorNumber, ProcessorNumber);
    fprintf(log_fptr, "TotalRefinement = %"ISYM"\n", TotalRefinement);
  }
 
  /* Determine if the data should be loaded in or not. */
 
  int ReadData = TRUE, Offset[] = {0,0,0};
 
  if (ParallelRootGridIO == TRUE && TotalRefinement == 1)
    ReadData = FALSE;
 
  if (io_log) fprintf(log_fptr, "ReadData = %"ISYM"\n", ReadData);
 
  /* Calculate buffer Offset (same as Grid unless doing ParallelRootGridIO
     (TotalRefinement = -1 if used as a signal that we should really load
     in the data regardless of the value of ParallelRootGridIO). */
 
  if (ParallelRootGridIO == TRUE && TotalRefinement == -1)
    for (dim = 0; dim < GridRank; dim++)
      Offset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/
			 CellWidth[dim][0]);
 
  /*----------------------------------------------------*/
  /* Create baryon fields. */
 
  NumberOfBaryonFields = 0;
  if (TurbulenceSimulationVelocityNames[0] != NULL) {
    FieldType[NumberOfBaryonFields++] = Density;
    FieldType[NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;
    FieldType[NumberOfBaryonFields++] = Velocity1;
    vel = NumberOfBaryonFields - 1;
    if (GridRank > 1)
      FieldType[NumberOfBaryonFields++] = Velocity2;
    if (GridRank > 2)
      FieldType[NumberOfBaryonFields++] = Velocity3;
  }
 
  /* Set the subgrid static flag. */
 
  SubgridsAreStatic = TurbulenceSimulationSubgridsAreStatic;
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber == MyProcessorNumber) {
 
  /* Skip following if NumberOfBaryonFields == 0. */
 
  if (NumberOfBaryonFields > 0) {
 
  /* Determine the size of the fields. */
 
  int size = 1;
 
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Allocate space for the fields. */
 
  if (ReadData == TRUE) {
    if (io_log) fprintf(log_fptr, "Allocate %"ISYM" fields, %"ISYM" floats per field\n",
			NumberOfBaryonFields, size);
    for (dim = 0; dim < GridRank; dim++)
      if (io_log) fprintf(log_fptr, "  Field dim %"ISYM" size %"ISYM"\n", dim, GridDimension[dim]);
    printf("Allocating %"ISYM" baryon fields of size %"ISYM"\n", NumberOfBaryonFields, size);
    for (int field = 0; field < NumberOfBaryonFields; field++)
      BaryonField[field] = new float[size];
  }
 
  /* If set, allocate space for RandomForcing fields. */
 
  if (RandomForcing == TRUE && ReadData == TRUE)
      for (int dim = 0; dim < GridRank; dim++)
        if (RandomForcingField[dim] == NULL)
          RandomForcingField[dim] = new float[size];
 
  /* Read the density field. */
 
  if (TurbulenceSimulationDensityName != NULL && ReadData)
    if (ReadFile(TurbulenceSimulationDensityName, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[0],
              &tempbuffer, 0, 1) == FAIL) {
      ENZO_FAIL("Error reading density field.\n");
    }
 
  /* Read the total energy field. */
 
  if (TurbulenceSimulationTotalEnergyName != NULL && ReadData)
    if (ReadFile(TurbulenceSimulationTotalEnergyName, GridRank,
		    GridDimension, GridStartIndex, GridEndIndex, Offset,
		    BaryonField[1], &tempbuffer, 0, 1) == FAIL) {
      ENZO_FAIL("Error reading total energy field.\n");
    }
 
  /* Read the gas energy field. */
 
  if (TurbulenceSimulationGasEnergyName != NULL && DualEnergyFormalism &&
      ReadData)
    if (ReadFile(TurbulenceSimulationGasEnergyName, GridRank, GridDimension,
		     GridStartIndex, GridEndIndex, Offset, BaryonField[2],
		     &tempbuffer, 0, 1) == FAIL) {
      ENZO_FAIL("Error reading gas energy field.\n");
    }
 
  /* Read the velocity fields. */
 
  if (TurbulenceSimulationVelocityNames[0] != NULL && ReadData)
    for (dim = 0; dim < GridRank; dim++)
      if (ReadFile(TurbulenceSimulationVelocityNames[dim], GridRank,
		   GridDimension, GridStartIndex, GridEndIndex, Offset,
		   BaryonField[vel+dim], &tempbuffer, 0, 1) == FAIL) {
	//	   BaryonField[vel+dim], &tempbuffer, dim, 3) == FAIL) {
	ENZO_VFAIL("Error reading velocity field %"ISYM".\n", dim)
      }
 
  /* Get RandomForcing data */
 
  if (RandomForcing == TRUE && TurbulenceSimulationRandomForcingNames[0] != NULL) {
 
  /* Read the random forcing fields. */
 
    if (ReadData)
      for (dim = 0; dim < GridRank; dim++)
	if (ReadFile(TurbulenceSimulationRandomForcingNames[dim], GridRank,
			  GridDimension, GridStartIndex, GridEndIndex, Offset,
			  RandomForcingField[dim], &tempbuffer, 0, 1) == FAIL) {
	  //		  RandomForcingField[dim], &tempbuffer, dim, 3) == FAIL) {
	  ENZO_VFAIL("Error reading RandomForcing field %"ISYM".\n", dim)
      }
  }
 
  else {
 
  /* OR: copy random forcing fields from initial velocities. */
 
    if (ReadData == TRUE && RandomForcing == TRUE )
    for (dim = 0; dim < GridRank; dim++)
      for (i = 0; i < size; i++)
        RandomForcingField[dim][i] = BaryonField[vel+dim][i];
  }
 
   /* If they were not read in above, set the total & gas energy fields now. */
 
  if (TurbulenceSimulationVelocityNames[0] != NULL && ReadData) {
    if (TurbulenceSimulationDensityName == NULL)
      for (i = 0; i < size; i++)
        BaryonField[0][i] = TurbulenceSimulationInitialDensity;
 
    if (TurbulenceSimulationTotalEnergyName == NULL)
      for (i = 0; i < size; i++)
        BaryonField[1][i] = TurbulenceSimulationInitialTemperature/(Gamma-1.);
 
    if (TurbulenceSimulationGasEnergyName == NULL && DualEnergyFormalism)
      for (i = 0; i < size; i++)
        BaryonField[2][i] = BaryonField[1][i];
 
    if (TurbulenceSimulationTotalEnergyName == NULL &&
        HydroMethod != Zeus_Hydro)
      for (dim = 0; dim < GridRank; dim++)
        for (i = 0; i < size; i++)
          BaryonField[1][i] +=
            0.5 * BaryonField[vel+dim][i] * BaryonField[vel+dim][i];
  }
 
  } // end: if (NumberOfBaryonFields > 0)
 
  } // end: if (ProcessorNumber == MyProcessorNumber)
 
  OldTime = Time;
 
  if (io_log) fclose(log_fptr);

 
  return SUCCESS;
}
