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
              float TurbulenceSimulationInitialDensityPerturbationAmplitude,
                          float TurbulenceSimulationInitialTemperature,
                          float TurbulenceSimulationInitialPressure,
                          float TurbulenceSimulationInitialMagneticField[],
                          char *TurbulenceSimulationMagneticNames[],
                          char *TurbulenceSimulationDensityName,
                          char *TurbulenceSimulationTotalEnergyName,
                          char *TurbulenceSimulationGasPressureName,
                          char *TurbulenceSimulationGasEnergyName,
                          char *TurbulenceSimulationVelocityNames[],
                          char *TurbulenceSimulationRandomForcingNames[],
                          int   TurbulenceSimulationSubgridsAreStatic,
                          int   TotalRefinement)
{
 
  /* declarations */
 
  int idim, dim, i, j, vel, ibx;
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
    FieldType[NumberOfBaryonFields++] = Density;
    if( EquationOfState == 0 )
      FieldType[NumberOfBaryonFields++] = TotalEnergy;
    if (DualEnergyFormalism)
      FieldType[NumberOfBaryonFields++] = InternalEnergy;
    FieldType[NumberOfBaryonFields++] = Velocity1;
    vel = NumberOfBaryonFields - 1;
    if (GridRank > 1)
      FieldType[NumberOfBaryonFields++] = Velocity2;
    if (GridRank > 2)
      FieldType[NumberOfBaryonFields++] = Velocity3;
    if (UseMHD) {
        ibx = NumberOfBaryonFields;
      FieldType[NumberOfBaryonFields++] = Bfield1;
      FieldType[NumberOfBaryonFields++] = Bfield2;
      FieldType[NumberOfBaryonFields++] = Bfield3;
    }
    if(HydroMethod == MHD_RK){
      FieldType[NumberOfBaryonFields++] = PhiField;
    }
   
    int idrivex, idrivey, idrivez;
    if (UseDrivingField && (HydroMethod == HD_RK || HydroMethod == MHD_RK)) {
      idrivex = NumberOfBaryonFields;
      idrivey = idrivex + 1;
      idrivez = idrivex + 2;
      FieldType[NumberOfBaryonFields++] = DrivingField1;
      FieldType[NumberOfBaryonFields++] = DrivingField2;
      FieldType[NumberOfBaryonFields++] = DrivingField3;
    }

  if( WritePotential )
      FieldType[NumberOfBaryonFields++] = GravPotential; 
 
  /* Set the subgrid static flag. */
 
  SubgridsAreStatic = TurbulenceSimulationSubgridsAreStatic;
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber == MyProcessorNumber) {
 
  /* Skip following if NumberOfBaryonFields == 0. */
 
  if (NumberOfBaryonFields > 0) {

      int DensNum = -1, GENum = -1, Vel1Num = -1, 
                         Vel2Num=-1, Vel3Num=-1, TENum=-1,
                         B1Num=-1, B2Num=-1, B3Num=-1, PhiNum=-1;
      IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, 
                           Vel2Num, Vel3Num, TENum,
                           B1Num, B2Num, B3Num, PhiNum);

 
  /* Determine the size of the fields. */
 
  int size = 1;
 
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Allocate space for the fields. */
 
  if (ReadData == TRUE) {
    this->AllocateGrids();
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
  if (TurbulenceSimulationTotalEnergyName != NULL && ReadData
      && EquationOfState == 0
      )

    if (ReadFile(TurbulenceSimulationTotalEnergyName, GridRank,
            GridDimension, GridStartIndex, GridEndIndex, Offset,
            BaryonField[TENum], &tempbuffer, 0, 1) == FAIL) {
      ENZO_FAIL("Error reading total energy field.\n");
    }
 
  /* Read the gas pressure field. */
 
  if (TurbulenceSimulationGasPressureName != NULL && ReadData)
    if (ReadFile(TurbulenceSimulationGasPressureName, GridRank,
            GridDimension, GridStartIndex, GridEndIndex, Offset,
            BaryonField[TENum], &tempbuffer, 0, 1) == FAIL) {
      ENZO_FAIL("Error reading total energy field.\n");
    }
 
  /* Read the gas energy field. */
 
  if (TurbulenceSimulationGasEnergyName != NULL && DualEnergyFormalism &&
      ReadData)
    if (ReadFile(TurbulenceSimulationGasEnergyName, GridRank, GridDimension,
             GridStartIndex, GridEndIndex, Offset, BaryonField[GENum],
             &tempbuffer, 0, 1) == FAIL) {
      ENZO_FAIL("Error reading gas energy field.\n");
    }
 
  /* Read the velocity fields. */
   
  int part=0, npart = 0;
  if (TurbulenceSimulationVelocityNames[0] != NULL && ReadData)
    for (dim = 0; dim < GridRank; dim++){
  if( strcmp( TurbulenceSimulationVelocityNames[0] , TurbulenceSimulationVelocityNames[1] ) == 0 ){
    part = dim;
    npart = 3;
  }else{
    part = 0;
    npart = 1;
  }
      if (ReadFile(TurbulenceSimulationVelocityNames[dim], GridRank,
           GridDimension, GridStartIndex, GridEndIndex, Offset,
           BaryonField[Vel1Num+dim], &tempbuffer, part, npart) == FAIL) {
    ENZO_VFAIL("Error reading velocity field %"ISYM".\n", dim)
      }
      }
  /* Read the magnetic fields. */
 
  if (TurbulenceSimulationMagneticNames[0] != NULL && ReadData && UseMHD ){

    for (dim = 0; dim < GridRank; dim++){
  if( strcmp( TurbulenceSimulationMagneticNames[0] , TurbulenceSimulationMagneticNames[1] ) == 0 ){
    part = dim;
    npart = 3;
  }else{
    part = 0;
    npart = 1;
  }
  if( UseMHDCT == TRUE ){
      if (ReadFile(TurbulenceSimulationMagneticNames[dim], GridRank,
                        MagneticDims[dim], MHDStartIndex[dim], MHDEndIndex[dim], Offset,
                        MagneticField[dim], &tempbuffer, part, npart) == FAIL) {
    ENZO_VFAIL("Error reading magnetic field %"ISYM".\n", dim)
      }
  }//MHDCT
  if( UseMHD ){
      if (ReadFile(TurbulenceSimulationMagneticNames[dim], GridRank,
           GridDimension, GridStartIndex, GridEndIndex, Offset,
           BaryonField[B1Num+dim], &tempbuffer, part, npart) == FAIL) {
    ENZO_VFAIL("Error reading magnetic field %"ISYM".\n", dim)
      }
      }//UseMHD

    }//dim (field)
  if( UseMHDCT ){
    this->CenterMagneticField();
      }
  }
  /* Get RandomForcing data */
 
  if (RandomForcing == TRUE && TurbulenceSimulationRandomForcingNames[0] != NULL) {
 
  /* Read the random forcing fields. */
 
    if (ReadData)
      for (dim = 0; dim < GridRank; dim++)
    if (ReadFile(TurbulenceSimulationRandomForcingNames[dim], GridRank,
              GridDimension, GridStartIndex, GridEndIndex, Offset,
              RandomForcingField[dim], &tempbuffer, 0, 1) == FAIL) {
      ENZO_VFAIL("Error reading RandomForcing field %"ISYM".\n", dim)
      }
  } else {
 
  /* OR: copy random forcing fields from initial velocities. */
    if( ReadData ){ 
        if (TurbulenceSimulationVelocityNames[0] == NULL) {
          //OR compute a brand new one!
          if( this->ComputeRandomForcingFields(0) == FAIL ){
            fprintf(stderr," Error in ComputeRandomForcingFields\n"); return FAIL;}
          
        }
      }//read
 
      if (ReadData == TRUE && RandomForcing == TRUE){
        for (dim = 0; dim < GridRank; dim++)
          for (i = 0; i < size; i++)
            RandomForcingField[dim][i] = BaryonField[vel+dim][i];
      }
  }
 
   /* If they were not read in above, set the total & gas energy fields now. */
 
  if (ReadData) {
    if (TurbulenceSimulationDensityName == NULL)
      for (i = 0; i < size; i++){
        BaryonField[DensNum][i] = TurbulenceSimulationInitialDensity;
      }
            

      for (i = 0; i < size; i++){
        BaryonField[DensNum][i] *= 1 + TurbulenceSimulationInitialDensityPerturbationAmplitude*
            ( (float)rand()/(float)(RAND_MAX)   - 0.5 );
      }

    if( UseMHD ){
        for( dim = 0; dim < 3; dim++){
            if( TurbulenceSimulationMagneticNames[ dim ] == NULL)
                for (i = 0; i < size; i++)
                     BaryonField[B1Num + dim ][i] =TurbulenceSimulationInitialMagneticField[dim];
        }
        for( i=0; i< size; i++)
            BaryonField[PhiNum][i]= 0;
    }
    if( UseMHDCT == TRUE && TurbulenceSimulationMagneticNames[0] == NULL ){
      for(int field=0;field<3;field++){
        MagneticField[field] = new float[MagneticSize[field]];
        for(i=0;i<MagneticSize[field];i++)
          MagneticField[field][i]=TurbulenceSimulationInitialMagneticField[field];
      }//field
    }
    if( EquationOfState == 0 ) {

      if (TurbulenceSimulationTotalEnergyName == NULL 
          && TurbulenceSimulationGasPressureName != NULL){
          for(i=0;i<size;i++)
            BaryonField[1][i] /= (BaryonField[0][i]*(Gamma -1));
      }
    
    if (TurbulenceSimulationTotalEnergyName == NULL && TurbulenceSimulationGasPressureName == NULL){
        if( TurbulenceSimulationInitialTemperature != FLOAT_UNDEFINED ){
            for (i = 0; i < size; i++)
              BaryonField[1][i] = TurbulenceSimulationInitialTemperature/(Gamma-1.);
        }
        if( TurbulenceSimulationInitialPressure != FLOAT_UNDEFINED ){
            for( i=0;i<size;i++){
                BaryonField[1][i] = TurbulenceSimulationInitialPressure/(Gamma - 1.)/
                    BaryonField[0][i];
            }
        }
    }
 
      if (TurbulenceSimulationGasEnergyName == NULL && DualEnergyFormalism){
      for (i = 0; i < size; i++)
        BaryonField[2][i] = BaryonField[1][i];
      }
 
    if (TurbulenceSimulationTotalEnergyName == NULL &&
        HydroMethod != Zeus_Hydro){

      for (dim = 0; dim < GridRank; dim++)
        for (i = 0; i < size; i++){
          BaryonField[1][i] +=
            0.5 * BaryonField[vel+dim][i] * BaryonField[vel+dim][i];
        }
      if(UseMHD){
          for(i=0;i<size;i++){
              BaryonField[TENum][i] += (0.5*(BaryonField[B1Num][i]*BaryonField[B1Num][i]+
                                            BaryonField[B2Num][i]*BaryonField[B2Num][i]+
                                            BaryonField[B3Num][i]*BaryonField[B3Num][i])/
                                        BaryonField[DensNum][i]);

          }
      }
    }
    }//EquationOfState  
  }//ReadData
 
  } // end: if (NumberOfBaryonFields > 0)
 
  } // end: if (ProcessorNumber == MyProcessorNumber)
 
  OldTime = Time;
 
  if (io_log) fclose(log_fptr);

 
  return SUCCESS;
}
