/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COSMOLOGY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness, March 2004
/              Changed default initialization of Total_Energy
/              The divisors DEFAULT_MU and Gamma-1.0 were missing,
/              leading to an effectively incorrect initial T.
/              64-bit Inits
/  modified2:  Robert Harkness, March 2006
/              Use MPICH V1.x Dims_create code
/  modified3:  Robert Harkness, February 2007
/              Read pre-decomposed grids
/  modified4:  Robert Harkness, April 2008
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include "mpi.h"
 
#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "fortran.def"
void my_exit(int status);
 
#ifdef PROTO  // Remove troublesome HDF PROTO declaration
#undef PROTO
#endif
 
#define DEFAULT_MU 0.6  // for temperature field
 
// HDF5 function prototypes
 

 
// Function prototypes
 
int ReadFile(char *name, int Rank, int Dims[], int StartIndex[],
                  int EndIndex[], int BufferOffset[], float *buffer,
                  inits_type **tempbuffer, int Part, int Npart);

int ReadGridFile(char *name, int Rank, int Dims[], int StartIndex[],
                      int EndIndex[], int BufferOffset[], float *buffer,
                      inits_type **tempbuffer, int Part, int Npart);

int ReadIntFile(char *name, int Rank, int Dims[], int StartIndex[],
                     int EndIndex[], int BufferOffset[], int *buffer,
                     int **tempbuffer, int Part, int Npart);
 
int ReadAttr(char *Fname, int *Rank, int Dims[], int *NSeg, int *LSeg, FILE *log_fptr);
 
void fcol(float *x, int n, int m, FILE *log_fptr);
void pcol(FLOAT *x, int n, int m, FILE *log_fptr);
void icol(int *x, int n, int m, FILE *log_fptr);
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);

int Enzo_Dims_create(int nnodes, int ndims, int *dims); 
 
 

 
int grid::CosmologySimulationInitializeGrid(
			  float CosmologySimulationOmegaBaryonNow,
			  float CosmologySimulationOmegaCDMNow,
			  float CosmologySimulationInitialTemperature,
			  char *CosmologySimulationDensityName,
			  char *CosmologySimulationTotalEnergyName,
			  char *CosmologySimulationGasEnergyName,
			  char *CosmologySimulationVelocityNames[],
			  char *CosmologySimulationParticlePositionName,
			  char *CosmologySimulationParticleVelocityName,
			  char *CosmologySimulationParticleMassName,
			  char *CosmologySimulationParticleTypeName,
			  int   CosmologySimulationSubgridsAreStatic,
			  int   TotalRefinement,
			  float CosmologySimulationInitialFractionHII,
			  float CosmologySimulationInitialFractionHeII,
			  float CosmologySimulationInitialFractionHeIII,
			  float CosmologySimulationInitialFractionHM,
			  float CosmologySimulationInitialFractionH2I,
			  float CosmologySimulationInitialFractionH2II,
			  int   UseMetallicityField,
			  int  &CurrentParticleNumber)
{
 
 
  int idim, dim, i, j, vel;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, MetalNum;
 
  int ExtraField[2];
 
  inits_type *tempbuffer = NULL;
  int *int_tempbuffer = NULL;
 
  hid_t       file_id, dset_id, attr_id, type_id, int_type_id;
  hid_t       mem_dsp_id, file_dsp_id, attr_dsp_id;
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  hsize_t     xfer_size;
  hsize_t     Slab_Dims[4];
  int         Slab_Rank;
 
  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     slab_stride[4], slab_count[4], slab_block[4];
  hsize_t     attr_count;
 
  hssize_t    mem_offset;
  hssize_t    slab_offset[4];
 
  int NSeg, LSeg, Part;
 
  int component_rank_attr;
  int component_size_attr;
  int field_rank_attr;
  int field_dims_attr[3];
 
  int enzo_layout[3];
  int A, B, C, I, J, K, M;
  int NOP, RNK;
  int mpi_layout[3];
  int fln, ext, lex;

  char *Extension;
  char pid[MAX_TASK_TAG_SIZE];

  MPI_Arg mpi_size;

  FILE *log_fptr;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
  if ( NumberOfProcessors > 64 )
  {
     if ( (ParallelRootGridIO != TRUE) || (ParallelParticleIO != TRUE) )
     {
       ENZO_FAIL("ParallelRootGridIO and ParallelParticleIO MUST be set for > 64 cpus!\n");
     }
  }
 
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  char *logname = new char[MAX_NAME_LENGTH];
  strcpy(logname, "CSlog.");
  strcat(logname, pid);
 
  if (io_log) log_fptr = fopen(logname, "a");
 
  delete [] logname;
 
  if (io_log) {
    fprintf(log_fptr, "\n");
    fprintf(log_fptr, "CSIG ParallelRootGridIO = %"ISYM"\n", ParallelRootGridIO);
    fprintf(log_fptr, "Processor %"ISYM", Target processor %"ISYM"\n", MyProcessorNumber, ProcessorNumber);
    fprintf(log_fptr, "TotalRefinement = %"ISYM"\n", TotalRefinement);
  }
 
  int ii = sizeof(inits_type);
 
  //  fprintf(stderr, "CSIG size of inits_type is %"ISYM"\n", ii);
 
  switch(ii)
  {
 
    case 4:
      type_id = HDF5_R4;
      break;
 
    case 8:
      type_id = HDF5_R8;
      break;
 
    default:
      type_id = HDF5_R4;
 
  }
 
  int_type_id = HDF5_INT;
 
  // Determine if the data should be loaded in or not
 
  int ReadData = TRUE, Offset[] = {0,0,0};
 
  if (ParallelRootGridIO == TRUE && TotalRefinement == 1)
    ReadData = FALSE;
 
  if (io_log) fprintf(log_fptr, "ReadData = %"ISYM"\n", ReadData);
 
  // Calculate buffer Offset (same as Grid unless doing ParallelRootGridIO
  // (TotalRefinement = -1 if used as a signal that we should really load
  // in the data regardless of the value of ParallelRootGridIO)
 
  if (ParallelRootGridIO == TRUE && TotalRefinement == -1)
    for (dim = 0; dim < GridRank; dim++)
      Offset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/CellWidth[dim][0]);
 
  // if (io_log) fprintf(log_fptr, "CSIG Offsets %"ISYM" %"ISYM" %"ISYM"\n", Offset[0], Offset[1], Offset[2]);
 
  // Create baryon fields (unless they are not needed)
 
  NumberOfBaryonFields = 0;

  if (CosmologySimulationDensityName != NULL) {

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
    if (MultiSpecies) {
      FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
      FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
      FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
      FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
      FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
      FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
      if (MultiSpecies > 1) {
	FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
	FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
	FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
      }
      if (MultiSpecies > 2) {
	FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
	FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
	FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
      }
    }
    if (UseMetallicityField) {
      FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity;
      FieldType[ExtraField[0] = NumberOfBaryonFields++] = ExtraType0;
      FieldType[ExtraField[1] = NumberOfBaryonFields++] = ExtraType1;
    }
    if (WritePotential)
      FieldType[NumberOfBaryonFields++] = GravPotential;
  }
 
  // Set the subgrid static flag
 
  SubgridsAreStatic = CosmologySimulationSubgridsAreStatic;




  // Return if this doesn't concern us
 
  if (ProcessorNumber == MyProcessorNumber) {


 

    // Skip following if NumberOfBaryonFields == 0
 
    if (NumberOfBaryonFields > 0) {
 
    // Get the cosmology units so we can convert temperature later
 
      float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1, 
	VelocityUnits=1;

      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		   &TimeUnits, &VelocityUnits, 
		   InitialTimeInCodeUnits) == FAIL) {
        ENZO_FAIL("Error in GetUnits.\n");
      }
 
      // Determine the size of the fields
 
      int size = 1;
 
      for (dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];
 
      // Allocate space for the fields
 
      if (ReadData == TRUE)
      {
        if (io_log) fprintf(log_fptr, "Allocate %"ISYM" fields, %"ISYM" floats per field\n", NumberOfBaryonFields, size);
        for (dim=0; dim < GridRank; dim++)
        {
          if (io_log) fprintf(log_fptr, "  Field dim %"ISYM" size %"ISYM"\n", dim, GridDimension[dim]);
        }
      }
 
      if (ReadData == TRUE)
        printf("Allocating %"ISYM" baryon fields of size %"ISYM"\n", NumberOfBaryonFields, size);
 
      if (ReadData == TRUE)
        for (int field = 0; field < NumberOfBaryonFields; field++)
          BaryonField[field] = new float[size];

      int Old_Style = FALSE;

      if (Old_Style)
      { 
 
        // Read the density field
 
        if (CosmologySimulationDensityName != NULL && ReadData) {
          if (ReadFile(CosmologySimulationDensityName, GridRank, GridDimension,
                        GridStartIndex, GridEndIndex, Offset, BaryonField[0],
                        &tempbuffer, 0, 1) == FAIL) {
            ENZO_FAIL("Error reading density field.\n");
          }
        }
 
        // Read the total energy field
 
        if (CosmologySimulationTotalEnergyName != NULL && ReadData) {
          if (ReadFile(CosmologySimulationTotalEnergyName, GridRank,
		    GridDimension, GridStartIndex, GridEndIndex, Offset,
		    BaryonField[1], &tempbuffer, 0, 1) == FAIL) {
            ENZO_FAIL("Error reading total energy field.\n");
          }
        }
 
        // Read the gas energy field
 
        if (CosmologySimulationGasEnergyName != NULL && DualEnergyFormalism && ReadData) {
          if (ReadFile(CosmologySimulationGasEnergyName, GridRank, GridDimension,
		 GridStartIndex, GridEndIndex, Offset, BaryonField[2],
		 &tempbuffer, 0, 1) == FAIL) {
            ENZO_FAIL("Error reading gas energy field.\n");
          }
        }
 
        // Read the velocity fields
 
        if (CosmologySimulationVelocityNames[0] != NULL && ReadData) {
          for (dim = 0; dim < GridRank; dim++)
            if (ReadFile(CosmologySimulationVelocityNames[dim], GridRank,
		   GridDimension, GridStartIndex, GridEndIndex, Offset,
		   BaryonField[vel+dim], &tempbuffer, dim, 3) == FAIL) {
              ENZO_VFAIL("Error reading velocity field %"ISYM".\n", dim)
            }
        }

      } else {

        // read with new stuff

 
        //  Generate subgrid to MPI task map for filenames

 
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        NOP = mpi_size;
        RNK = 3;

        mpi_layout[0] = 0;
        mpi_layout[1] = 0;
        mpi_layout[2] = 0;

        Enzo_Dims_create(NOP, RNK, mpi_layout);

        for (idim = 0; idim < GridRank; idim++)
        {
          enzo_layout[idim] = mpi_layout[GridRank-1-idim];
        }

        if (MyProcessorNumber == ROOT_PROCESSOR) {
          fprintf(stderr, "ENZO_layout %"ISYM" x %"ISYM" x %"ISYM"\n", enzo_layout[0], enzo_layout[1], enzo_layout[2]);
        }

        A = enzo_layout[0];
        B = enzo_layout[1];
        C = enzo_layout[2];
        I = nint( GridLeftEdge[0]/((DomainRightEdge[0]-DomainLeftEdge[0])/((float) A)));
        J = nint( GridLeftEdge[1]/((DomainRightEdge[1]-DomainLeftEdge[1])/((float) B)));
        K = nint( GridLeftEdge[2]/((DomainRightEdge[2]-DomainLeftEdge[2])/((float) C)));
//      M = ((I*B*C) + J*C) + K;
        M = ((K*B*A) + J*A) + I;
 
//      if (io_log) fprintf(log_fptr, "ABC %"ISYM" %"ISYM" %"ISYM";  IJK %"ISYM" %"ISYM" %"ISYM";  M = %"ISYM"\n", A,B,C,I,J,K,M);
 
        sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, M);

 
        lex = 0;
        fln = strlen(CosmologySimulationDensityName);
        ext = strcspn(CosmologySimulationDensityName, ".");
 
        if ( fln-ext > 0 )
        {
          lex = strlen(strstr(CosmologySimulationDensityName, "."));
          Extension = new char[lex+1];
          strcpy(Extension, strstr(CosmologySimulationDensityName, "."));
          lex = strlen(Extension);
        }
 
        char *GD = new char[MAX_NAME_LENGTH];
        strcpy(GD, "GD");
        strcat(GD, pid);
        if (lex > 0)
          strcat(GD, Extension);
 
        char *GV = new char[MAX_NAME_LENGTH];
        strcpy(GV, "GV");
        strcat(GV, pid);
        if (lex > 0)
          strcat(GV, Extension);
 
        if (io_log) fprintf(log_fptr, "Proc %"ISYM" reads density grid %s\n", MyProcessorNumber, GD);
        if (io_log) fprintf(log_fptr, "Proc %"ISYM" reads velocity grid %s\n", MyProcessorNumber, GV);
 
        if ( fln-ext > 0 )
          delete [] Extension;

        // Read the density field
 
        if (CosmologySimulationDensityName != NULL && ReadData) {
          if (ReadGridFile(GD, GridRank, GridDimension,
                        GridStartIndex, GridEndIndex, Offset, BaryonField[0],
                        &tempbuffer, 0, 1) == FAIL) {
            ENZO_FAIL("Error reading density field.\n");
          }
        }
 
        // Read the velocity fields
 
        if (CosmologySimulationVelocityNames[0] != NULL && ReadData) {
          for (dim = 0; dim < GridRank; dim++)
            if (ReadGridFile(GV, GridRank, GridDimension,
                          GridStartIndex, GridEndIndex, Offset, BaryonField[vel+dim],
                          &tempbuffer, dim, 3) == FAIL) {
              ENZO_VFAIL("Error reading velocity field %"ISYM".\n", dim)
            }
        }


        // end of read with new stuff

      }


      // If using multi-species, set the fields
 
      if (MultiSpecies && ReadData) {
        for (i = 0; i < size; i++) {
 
          BaryonField[HIINum][i] = CosmologySimulationInitialFractionHII *
            CoolData.HydrogenFractionByMass * BaryonField[0][i] *
            sqrt(OmegaMatterNow)/
            (CosmologySimulationOmegaBaryonNow*HubbleConstantNow);
 
          BaryonField[HeIINum][i] = CosmologySimulationInitialFractionHeII*
            BaryonField[0][i] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
          BaryonField[HeIIINum][i] = CosmologySimulationInitialFractionHeIII*
            BaryonField[0][i] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
          BaryonField[HeINum][i] =
            (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][i] -
              BaryonField[HeIINum][i] - BaryonField[HeIIINum][i];
 
          if (MultiSpecies > 1) {
            BaryonField[HMNum][i] = CosmologySimulationInitialFractionHM*
              BaryonField[HIINum][i]*
              POW(CosmologySimulationInitialTemperature,float(0.88));
            BaryonField[H2IINum][i] = CosmologySimulationInitialFractionH2II*2.0*
              BaryonField[HIINum][i]*
              POW(CosmologySimulationInitialTemperature,float(1.8));
            BaryonField[H2INum][i] = CosmologySimulationInitialFractionH2I*
              BaryonField[0][i]*CoolData.HydrogenFractionByMass*POW(301.0,5.1)*
              POW(OmegaMatterNow, float(1.5))/
              CosmologySimulationOmegaBaryonNow/HubbleConstantNow*2.0;
          }
 
          BaryonField[HINum][i] = CoolData.HydrogenFractionByMass*BaryonField[0][i]
                                  - BaryonField[HIINum][i];

          if (MultiSpecies > 1)
            BaryonField[HINum][i] -= BaryonField[HMNum][i]
                                     + BaryonField[H2IINum][i]
                                     + BaryonField[H2INum][i];
 
          BaryonField[DeNum][i] = BaryonField[HIINum][i] +
                                  0.25*BaryonField[HeIINum][i] + 0.5*BaryonField[HeIIINum][i];

          if (MultiSpecies > 1)
            BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] -
                                     BaryonField[HMNum][i];
 
          // Set Deuterium species (assumed to be negligible)
 
          if (MultiSpecies > 2) {
            BaryonField[DINum][i] = CoolData.DeuteriumToHydrogenRatio*
                                    BaryonField[HINum][i];
            BaryonField[DIINum][i] = CoolData.DeuteriumToHydrogenRatio*
                                     BaryonField[HIINum][i];
            BaryonField[HDINum][i] = CoolData.DeuteriumToHydrogenRatio*
                                     BaryonField[H2INum][i];
          }
 
        }

      }
 
      // If using metallicity, set the field
 
      if (UseMetallicityField && ReadData) {
        for (i = 0; i < size; i++) {
          BaryonField[MetalNum][i] = 1.0e-10 * BaryonField[0][i];
          BaryonField[ExtraField[0]][i] = 1.0e-10 * BaryonField[0][i];
          BaryonField[ExtraField[1]][i] = 1.0e-10 * BaryonField[0][i];
        }
      }
 
      // If they were not read in above, set the total & gas energy fields now
 
      if (CosmologySimulationDensityName != NULL && ReadData) {

        if (CosmologySimulationTotalEnergyName == NULL)
          for (i = 0; i < size; i++)
            BaryonField[1][i] = CosmologySimulationInitialTemperature/
                                TemperatureUnits/DEFAULT_MU/(Gamma-1.0);
 
        if (CosmologySimulationGasEnergyName == NULL && DualEnergyFormalism)
          for (i = 0; i < size; i++)
            BaryonField[2][i] = BaryonField[1][i];
 
        if (CosmologySimulationTotalEnergyName == NULL && HydroMethod != Zeus_Hydro)
          for (dim = 0; dim < GridRank; dim++)
            for (i = 0; i < size; i++)
              BaryonField[1][i] += 0.5 * BaryonField[vel+dim][i] * BaryonField[vel+dim][i];

      }




    } // end: if (NumberOfBaryonFields > 0)
 




 
    // Read particle fields if the was a name specified
 
    // This routine reads data from Inits and Ring: 32- or 64-bit

/*  REPLACED ABOVE 
    int ii = sizeof(inits_type);
 
    //  fprintf(stderr, "CSIG size of inits_type is %"ISYM"\n", ii);
 
    switch(ii)
    {
 
      case 4:
        type_id = HDF5_R4;
        break;
 
      case 8:
        type_id = HDF5_R8;
        break;
 
      default:
        type_id = HDF5_R4;
 
    }
 
    int_type_id = HDF5_INT;
*/ 



 
    if (CosmologySimulationParticlePositionName != NULL && ReadData) {
 
      // Get the number of particles by reading the file attributes
 
      int TempInt, Dim[1], Start[1] = {0}, End[1], Zero[1] = {0};
 
      int TempIntArray[MAX_DIMENSION];
 
      ReadAttr(CosmologySimulationParticlePositionName,
                    &TempInt, TempIntArray, &NSeg, &LSeg, log_fptr);
 
      // Error check
 
      if (TempInt != 1) {
        ENZO_VFAIL("Rank (%"ISYM") is not one in file %s.\n", TempInt,
                CosmologySimulationParticlePositionName)
      }
 
      // If doing parallel root grid IO then read in the full list of particle
      // positions and construct a mask of those within this grid
 
      int TotalParticleCount = TempIntArray[0];
 
/*
      Eunsigned_int *BitMask = NULL;
      Eunsigned_int BitsPerInt;
      Eunsigned_int BitMaskSize;
      Eunsigned_int BitMaskTrue;
      Eunsigned_int *BitArray = NULL;
      Eunsigned_int TestBit;
      int XMask;
 
      BitsPerInt = 8 * sizeof(BitsPerInt);
      BitMaskSize = (TotalParticleCount/BitsPerInt)+1;
*/
 
      if (ParallelRootGridIO == TRUE && TotalRefinement == -1) {
 
//      for(i=0; i<GridRank;i++)
//      {
//        if (io_log) fprintf(log_fptr, "Proc %"ISYM" Offset %"ISYM" Left %10.4"FSYM" Right %10.4"FSYM"\n",
//                    MyProcessorNumber, Offset[i], GridLeftEdge[i], GridRightEdge[i]);
//      }
 
 
        int PreSortedParticles = 0;
        int NumSortedParticles = 0;
        int TotParticleCount = 0;
 
        if (ParallelParticleIO)
        {
          PreSortedParticles = 1;
        }
 
        if (PreSortedParticles == 1)
        {
 
          //  This section only used with pre-sorted particle input (ring i/o)
 
          //  Generate subgrid to MPI task map for filenames

          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          NOP = mpi_size;
          RNK = 3;

          mpi_layout[0] = 0;
          mpi_layout[1] = 0;
          mpi_layout[2] = 0;

          Enzo_Dims_create(NOP, RNK, mpi_layout);

          for (idim = 0; idim < GridRank; idim++)
          {
            enzo_layout[idim] = mpi_layout[GridRank-1-idim];
          }

          if (MyProcessorNumber == ROOT_PROCESSOR) {
            fprintf(stderr, "ENZO_layout %"ISYM" x %"ISYM" x %"ISYM"\n", enzo_layout[0], enzo_layout[1], enzo_layout[2]);
          }

          A = enzo_layout[0];
          B = enzo_layout[1];
          C = enzo_layout[2];
          I = nint( GridLeftEdge[0]/((DomainRightEdge[0]-DomainLeftEdge[0])/((float) A)));
          J = nint( GridLeftEdge[1]/((DomainRightEdge[1]-DomainLeftEdge[1])/((float) B)));
          K = nint( GridLeftEdge[2]/((DomainRightEdge[2]-DomainLeftEdge[2])/((float) C)));
//        M = ((I*B*C) + J*C) + K;
          M = ((K*B*A) + J*A) + I;
 
//        if (io_log) fprintf(log_fptr, "ABC %"ISYM" %"ISYM" %"ISYM";  IJK %"ISYM" %"ISYM" %"ISYM";  M = %"ISYM"\n", A,B,C,I,J,K,M);
 
          sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, M);

          if (debug) 
            printf("PreSortedParticles - ParallelRootGridIO\n");
 
          lex = 0;
          fln = strlen(CosmologySimulationParticlePositionName);
          ext = strcspn(CosmologySimulationParticlePositionName, ".");
 
          if ( fln-ext > 0 )
          {
            lex = strlen(strstr(CosmologySimulationParticlePositionName, "."));
            Extension = new char[lex+1];
            strcpy(Extension, strstr(CosmologySimulationParticlePositionName, "."));
            lex = strlen(Extension);
          }
 
          char *PPos = new char[MAX_NAME_LENGTH];
          strcpy(PPos, "PPos");
          strcat(PPos, pid);
          if (lex > 0)
            strcat(PPos, Extension);
 
          char *PVel = new char[MAX_NAME_LENGTH];
          strcpy(PVel, "PVel");
          strcat(PVel, pid);
          if (lex > 0)
            strcat(PVel, Extension);
 
          char *PPro = new char[MAX_NAME_LENGTH];
          strcpy(PPro, "PPro");
          strcat(PPro, pid);
          if (lex > 0)
            strcat(PPro, Extension);
 
          char *PMass = new char[MAX_NAME_LENGTH];
          strcpy(PMass, "PMass");
          strcat(PMass, pid);
          if (lex > 0)
            strcat(PMass, Extension);
 
          char *PType = new char[MAX_NAME_LENGTH];
          strcpy(PType, "PType");
          strcat(PType, pid);
          if (lex > 0)
            strcat(PType, Extension);
 
          if (io_log) fprintf(log_fptr, "Proc %"ISYM" reads particle grid %s\n", MyProcessorNumber, PPos);
          if (io_log) fprintf(log_fptr, "Proc %"ISYM" reads particle grid %s\n", MyProcessorNumber, PVel);
 
          if ( fln-ext > 0 )
            delete [] Extension;

 


          // Read attributes for parallel particle I/O
 
          file_id = H5Fopen(PPos, H5F_ACC_RDONLY, H5P_DEFAULT);
          fprintf(stderr, "H5Fopen %s on proc %"ISYM" status %"ISYM"\n", PPos, MyProcessorNumber, file_id);
            if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
            if( file_id == h5_error){my_exit(EXIT_FAILURE);}
 
          dset_id = H5Dopen(file_id, PPos);
            if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
            if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
          attr_id = H5Aopen_name(dset_id, "NumberOfParticles");
            if (io_log) fprintf(log_fptr, "H5Aopen (NumberOfParticles) id: %"ISYM"\n", attr_id);
            if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
          h5_status = H5Aread(attr_id, HDF5_INT, &NumSortedParticles);
            if (io_log) fprintf(log_fptr, "H5Aread NumSortedParticles: %"ISYM"\n", h5_status);
            if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
          h5_status = H5Aclose(attr_id);
            if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
            if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
          attr_id = H5Aopen_name(dset_id, "TotalParticleCount");
            if (io_log) fprintf(log_fptr, "H5Aopen (TotalParticleCount) id: %"ISYM"\n", attr_id);
            if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
          h5_status = H5Aread(attr_id, HDF5_INT, &TotParticleCount);
            if (io_log) fprintf(log_fptr, "H5Aread TotParticleCount: %"ISYM"\n", h5_status);
            if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
          h5_status = H5Aclose(attr_id);
            if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
            if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
          h5_status = H5Dclose(dset_id);
            if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
            if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
          h5_status = H5Fclose(file_id);
            if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
            if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
          if ( TotParticleCount != TotalParticleCount )
          {
            ENZO_FAIL("DISASTER! Inconsistent particle count\n");
          }
 
          NumberOfParticles = NumSortedParticles;
 
          // Allocate space for the particles
 
          printf("Calling AllocateNewParticles with %"ISYM"\n", NumberOfParticles);
 
          this->AllocateNewParticles(NumberOfParticles);
 
          // Allocate 32- or 64-bit buffer for ring input
 
          if (io_log) fprintf(log_fptr, "Allocate %"ISYM" inits_type for RingTemp\n", NumberOfParticles);
 
          printf("Allocating %"ISYM" (NumberOfParticles)\n",NumberOfParticles);
 
 
 
 
          // Allocate buffer for floating point particle data
 
          inits_type *RingTemp = new inits_type[NumberOfParticles];
 
          // Allocate buffer for integer particle data
 
          int *IntRingTemp = new int[NumberOfParticles];
 
 
 
 
          // Read pre-sorted particle positions

          if (CosmologySimulationParticlePositionName != NULL) {
 
            for (dim = 0; dim < GridRank; dim++)
            {
 
              Slab_Rank = 2;
              Slab_Dims[0] = 3;
              Slab_Dims[1] = NumSortedParticles;
 
              if ( NumSortedParticles == 0 )
              {
                Slab_Dims[1] = 1;
              }
 
              // Data in memory is considered 1D, stride 1, with zero offset
 
              mem_stride = 1;                   // contiguous elements
              mem_count = NumSortedParticles;   // number of elements in field
              mem_offset = 0;                   // zero offset in buffer
 
              if ( NumSortedParticles == 0 )
              {
                mem_count = 1;
              }
 
              // 1D memory model
 
              mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
                if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
                if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              // Data in the file is (1+Rank)D with Npart components per grid point.
              // Offset[0] is the component Part of Npart components.  Data for each
              // Part are contiguous in the file, so stride = 1.
 
              slab_stride[0] = 1;      // contiguous elements
              slab_count[0] = 1;       // one component per call
              slab_offset[0] = dim;    // component Part of Npart
 
              slab_stride[1] = 1;                   // contiguous elements
              slab_count[1] = NumSortedParticles;   // field dimensions
              slab_offset[1] = 0;                   // complete field, no offset
 
              if ( NumSortedParticles == 0 )
              {
                slab_count[1] = 1;
              }
 
              file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
                if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
                if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              file_id = H5Fopen(PPos, H5F_ACC_RDWR, H5P_DEFAULT);
              fprintf(stderr, "H5Fopen %s on proc %"ISYM" status %"ISYM"\n", PPos, MyProcessorNumber, file_id);
                if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
                if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              dset_id =  H5Dopen(file_id, PPos);
                if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
                if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              if ( NumSortedParticles > 0 )
              {
                h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) RingTemp);
                  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
              }
 
              h5_status = H5Dclose(dset_id);
                if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(mem_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(file_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Fclose(file_id);
                if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              // Convert to FLOAT
 
              for (j = 0; j < NumberOfParticles; j++)
              {
                ParticlePosition[dim][j] = RingTemp[j];
              }
 
            } // end of loop over dims
 
          }

 
 
 
          // Read pre-sorted particle velocities

          if (CosmologySimulationParticleVelocityName != NULL) {
 
            for (dim = 0; dim < GridRank; dim++)
            {
 
              Slab_Rank = 2;
              Slab_Dims[0] = 3;
              Slab_Dims[1] = NumSortedParticles;
 
              if ( NumSortedParticles == 0 )
              {
                Slab_Dims[1] = 1;
              }
 
              // Data in memory is considered 1D, stride 1, with zero offset
 
              mem_stride = 1;                   // contiguous elements
              mem_count = NumSortedParticles;   // number of elements in field
              mem_offset = 0;                   // zero offset in buffer
 
              if ( NumSortedParticles == 0 )
              {
                mem_count = 1;
              }
 
              // 1D memory model
 
              mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
                if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
                if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              // Data in the file is (1+Rank)D with Npart components per grid point.
              // Offset[0] is the component Part of Npart components.  Data for each
              // Part are contiguous in the file, so stride = 1.
 
              slab_stride[0] = 1;      // contiguous elements
              slab_count[0] = 1;       // one component per call
              slab_offset[0] = dim;    // component Part of Npart
 
              slab_stride[1] = 1;                   // contiguous elements
              slab_count[1] = NumSortedParticles;    // field dimensions
              slab_offset[1] = 0;                   // complete field, no offset
 
              if ( NumSortedParticles == 0 )
              {
                slab_count[1] = 1;
              }
 
              file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
                if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
                if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              file_id = H5Fopen(PVel, H5F_ACC_RDWR, H5P_DEFAULT);
              fprintf(stderr, "H5Fopen %s on proc %"ISYM" status %"ISYM"\n", PVel, MyProcessorNumber, file_id);
                if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
                if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              dset_id =  H5Dopen(file_id, PVel);
                if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
                if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              if ( NumSortedParticles > 0 )
              {
                h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) RingTemp);
                  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
              }
 
              h5_status = H5Dclose(dset_id);
                if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(mem_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(file_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Fclose(file_id);
                if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              // Convert to float
 
              for (j = 0; j < NumberOfParticles; j++)
              {
                ParticleVelocity[dim][j] = RingTemp[j];
              }
 
            } // end of loop over dims
 
          }

 
 
 
          // Read pre-sorted particle mass

          if (CosmologySimulationParticleMassName != NULL) {
 
            dim = 0;
 
            Slab_Rank = 2;
            Slab_Dims[0] = 1;
            Slab_Dims[1] = NumSortedParticles;
 
            if ( NumSortedParticles == 0 )
            {
              Slab_Dims[1] = 1;
            }
 
            // Data in memory is considered 1D, stride 1, with zero offset
 
            mem_stride = 1;                   // contiguous elements
            mem_count = NumSortedParticles;   // number of elements in field
            mem_offset = 0;                   // zero offset in buffer
 
            if ( NumSortedParticles == 0 )
            {
              mem_count = 1;
            }
 
            // 1D memory model
 
            mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
              if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            // Data in the file is (1+Rank)D with Npart components per grid point.
            // Offset[0] is the component Part of Npart components.  Data for each
            // Part are contiguous in the file, so stride = 1.
 
            slab_stride[0] = 1;      // contiguous elements
            slab_count[0] = 1;       // one component per call
            slab_offset[0] = dim;    // component Part of Npart
 
            slab_stride[1] = 1;                   // contiguous elements
            slab_count[1] = NumSortedParticles;    // field dimensions
            slab_offset[1] = 0;                   // complete field, no offset
 
            if ( NumSortedParticles == 0 )
            {
              slab_count[1] = 1;
            }
 
            file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
              if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
              if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            file_id = H5Fopen(PMass, H5F_ACC_RDWR, H5P_DEFAULT);
              if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
              if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            dset_id =  H5Dopen(file_id, PMass);
              if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
              if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            if ( NumSortedParticles > 0 )
            {
              h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) RingTemp);
                if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
            }
 
            h5_status = H5Dclose(dset_id);
              if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sclose(mem_dsp_id);
              if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sclose(file_dsp_id);
              if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Fclose(file_id);
              if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            // Convert to float
 
            for (j = 0; j < NumberOfParticles; j++)
            {
              ParticleMass[j] = RingTemp[j];
            }
 
          }

 
 
 
          // Read pre-sorted particle type

          if (CosmologySimulationParticleTypeName != NULL) {
 
            dim = 0;
 
            Slab_Rank = 2;
            Slab_Dims[0] = 1;
            Slab_Dims[1] = NumSortedParticles;
 
            if ( NumSortedParticles == 0 )
            {
              Slab_Dims[1] = 1;
            }
 
            // Data in memory is considered 1D, stride 1, with zero offset
 
            mem_stride = 1;                   // contiguous elements
            mem_count = NumSortedParticles;   // number of elements in field
            mem_offset = 0;                   // zero offset in buffer
 
            if ( NumSortedParticles == 0 )
            {
              mem_count = 1;
            }
 
            // 1D memory model
 
            mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
              if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            // Data in the file is (1+Rank)D with Npart components per grid point.
            // Offset[0] is the component Part of Npart components.  Data for each
            // Part are contiguous in the file, so stride = 1.
 
            slab_stride[0] = 1;      // contiguous elements
            slab_count[0] = 1;       // one component per call
            slab_offset[0] = dim;    // component Part of Npart
 
            slab_stride[1] = 1;                   // contiguous elements
            slab_count[1] = NumSortedParticles;    // field dimensions
            slab_offset[1] = 0;                   // complete field, no offset
 
            if ( NumSortedParticles == 0 )
            {
              slab_count[1] = 1;
            }
 
            file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
              if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
              if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
              if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            file_id = H5Fopen(PType, H5F_ACC_RDWR, H5P_DEFAULT);
              if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
              if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            dset_id =  H5Dopen(file_id, PType);
              if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
              if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            if ( NumSortedParticles > 0 )
            {
              h5_status = H5Dread(dset_id, int_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) IntRingTemp);
                if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
            }
 
            h5_status = H5Dclose(dset_id);
              if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sclose(mem_dsp_id);
              if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sclose(file_dsp_id);
              if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Fclose(file_id);
              if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            // Convert to float
 
            for (j = 0; j < NumberOfParticles; j++)
            {
              ParticleType[j] = IntRingTemp[j];
            }
 
          } 
 
 

          if (debug) 
            printf("De-allocate (TotalParticleCount)\n");
 
          delete [] RingTemp;
          delete [] IntRingTemp;
 
        } // End of || particle i/o for pre-sorted particles
 
 
 
        // Normal ||rgio - obsolete - use pre-sort
 
        if (PreSortedParticles == 0)
        {

          if (debug) 
            printf("UnsortedParticles - ParallelRootGridIO\n");
 
          if (io_log) fprintf(log_fptr, "CSIG ParallelRootGridIO\n");
          if (io_log) fprintf(log_fptr, "TotalParticleCount %"ISYM"\n", TotalParticleCount);
 
          Eunsigned_int *BitMask = NULL;
          Eunsigned_int BitsPerInt;
          Eunsigned_int BitMaskSize;
          Eunsigned_int BitMaskTrue;
          Eunsigned_int *BitArray = NULL;
          Eunsigned_int TestBit;
          int XMask;
 
          BitsPerInt = 8 * sizeof(BitsPerInt);
          BitMaskSize = (TotalParticleCount/BitsPerInt)+1;
 
          // Generate a mask indicating if the particle is inside or outside
 
          if (io_log) fprintf(log_fptr, "Allocate %"ISYM" ints for Mask\n", TotalParticleCount);
 
          BitMask = new Eunsigned_int[BitMaskSize];
 
          BitArray = new Eunsigned_int[BitsPerInt];
 
          for (i = 0; i < BitMaskSize; i++)
            BitMask[i] = 0;
 
          BitArray[BitsPerInt-1] = 1;
          BitMaskTrue = 1;

          for (i=BitsPerInt-1; i>0; i--)
          {
            BitArray[i-1] = 2 * BitArray[i];
            BitMaskTrue = (BitMaskTrue | BitArray[i-1]);
          }
 
          for (i = 0; i < BitMaskSize; i++)
            BitMask[i] = BitMaskTrue;
 
          Eunsigned_int MaskAddr;
          Eunsigned_int WordAddr;
          Eunsigned_int BitAddr;
 
//        if (io_log) fprintf(log_fptr, "Allocate %"ISYM" inits_type for TempField\n", TotalParticleCount);
 
//        inits_type *TempField = new inits_type[TotalParticleCount];
 
          End[0] = TotalParticleCount - 1;
          Dim[0] = TotalParticleCount;
 
          int ppbuffsize, pp1, pp2, ppoffset, ppcount;
 
          ppbuffsize = min(1024,TotalParticleCount);
 
//        printf("ppbuffsize = %"ISYM"\n", ppbuffsize);
 
          if (io_log) fprintf(log_fptr, "Allocate %"ISYM" inits_type for TempField\n", ppbuffsize);
 
          inits_type *TempField = new inits_type[ppbuffsize];
          int *IntTempField = new int[ppbuffsize];
 
          for (pp1 = 0; pp1 < TotalParticleCount; pp1 = pp1 + ppbuffsize)
          {
            pp2 = min(pp1+ppbuffsize-1,TotalParticleCount-1);
            ppoffset = pp1;
            ppcount = pp2-pp1+1;
//          if (io_log) fprintf(log_fptr, "PP1 = %"ISYM", PP2 = %"ISYM", PPoffset = %"ISYM", PPCount = %"ISYM"\n", pp1, pp2, ppoffset, ppcount);
          }
 
 
 
 
          int n;
 
          for (dim = 0; dim < GridRank; dim++) {
 
            n = 0; // counter for particles in grid
 
            for (pp1 = 0; pp1 < TotalParticleCount; pp1 = pp1 + ppbuffsize)
            {
              pp2 = min(pp1+ppbuffsize-1,TotalParticleCount-1);
              ppoffset = pp1;
              ppcount = pp2-pp1+1;
 
              Part = dim;
 
              if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", CosmologySimulationParticlePositionName);
 
              file_id = H5Fopen(CosmologySimulationParticlePositionName, H5F_ACC_RDONLY, H5P_DEFAULT);
                if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
                if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", CosmologySimulationParticlePositionName);
 
              dset_id = H5Dopen(file_id, CosmologySimulationParticlePositionName);
                if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
                if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
 
              if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Rank\n");
 
              attr_id = H5Aopen_name(dset_id, "Component_Rank");
                if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aread(attr_id, HDF5_INT, &component_rank_attr);
                if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aclose(attr_id);
                if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              if (io_log) fprintf(log_fptr, "COMPONENT_RANK %"ISYM"\n", component_rank_attr);
 
 
              if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Size\n");
 
              attr_id = H5Aopen_name(dset_id, "Component_Size");
                if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aread(attr_id, HDF5_INT, &component_size_attr);
                if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aclose(attr_id);
                if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              if (io_log) fprintf(log_fptr, "COMPONENT_SIZE %"ISYM"\n", component_size_attr);
 
 
              if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Rank\n");
 
              attr_id = H5Aopen_name(dset_id, "Rank");
                if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aread(attr_id, HDF5_INT, &field_rank_attr);
                if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aclose(attr_id);
                if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              if (io_log) fprintf(log_fptr, "RANK %"ISYM"\n", field_rank_attr);
 
 
              if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Dimensions\n");
 
              attr_id = H5Aopen_name(dset_id, "Dimensions");
                if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              attr_count = field_rank_attr;
 
              attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
                if (io_log) fprintf(log_fptr, "Attr_dsp_id: %"ISYM"\n", attr_dsp_id);
                if( attr_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aread(attr_id, HDF5_INT, field_dims_attr);
                if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aclose(attr_id);
                if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(attr_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              for (idim = 0; idim < field_rank_attr; idim++)
              {
                if (io_log) fprintf(log_fptr, "DIMS %"ISYM":  %"ISYM"\n", idim, (int) field_dims_attr[idim]);
              }
 
              xfer_size = 1;
 
              for ( idim = 0; idim < field_rank_attr; idim++ )
              {
                xfer_size = xfer_size * field_dims_attr[idim];
              }
 
              if (io_log) fprintf(log_fptr, "  Grid Elements %"ISYM"\n", (int) xfer_size);
 
 
 
 
              Slab_Rank = field_rank_attr+1;
 
              if ( Slab_Rank != 2 )
              {
                ENZO_FAIL(" MAJOR ERROR!! Particle Slab_Rank != 2\n");
              }
 
              Slab_Dims[0] = component_rank_attr;
              Slab_Dims[1] = field_dims_attr[0];
 
              if (io_log) fprintf(log_fptr, "  Extended Rank %"ISYM"\n", (int) Slab_Rank);
 
              for ( idim = 0; idim < Slab_Rank; idim++ )
              {
                if (io_log) fprintf(log_fptr, "    %"ISYM":  %"ISYM"\n", idim, (int) Slab_Dims[idim]);
              }
 
              slab_offset[0] = Part;
              slab_stride[0] = 1;
              slab_count[0] = 1;
              slab_block[0] = 1;
 
              xfer_size = ppcount;
 
              slab_offset[1] = ppoffset;
              slab_stride[1] = 1;
              slab_count[1] = ppcount;
              slab_block[1] = 1;
 
              for ( idim=0; idim < Slab_Rank; idim++)
              {
                if (io_log) fprintf(log_fptr, "Slab [%"ISYM"] Offset %"ISYM", Count %"ISYM", Stride %"ISYM", Block %"ISYM"\n",
                        idim, (int) slab_offset[idim], (int) slab_count[idim],
                             (int) slab_stride[idim], (int) slab_block[idim]);
              }

              if (io_log) fprintf(log_fptr, "Xfer_size %"ISYM"\n", (int) xfer_size);
              if (io_log) fprintf(log_fptr, "Comp_size %"ISYM"\n", (int) component_size_attr);
 
              // Data in memory is considered 1D, stride 1, with zero offset
 
              mem_stride = 1;
              mem_count = xfer_size;
              mem_offset = 0;
              mem_block = 1;
 
              // 1D memory model
 
              mem_dsp_id = H5Screate_simple((Eint32) 1, &xfer_size, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
                if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, &mem_block);
                if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
                if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
                if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) TempField);
                if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", (int) h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(mem_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(file_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Dclose(dset_id);
                if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Fclose(file_id);
                if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              for (i = pp1, j = 0; i < pp1 + ppcount; i++, j++)
              {
 
//              if (TempField[i] < GridLeftEdge[dim] || TempField[i] > GridRightEdge[dim] )
//              {
//                Mask[i] = FALSE;
//              }
 
                if (TempField[j] < GridLeftEdge[dim] || TempField[j] > GridRightEdge[dim] )
                {
                  MaskAddr = i;
                  WordAddr = MaskAddr/BitsPerInt;
                  BitAddr  = MaskAddr%BitsPerInt;
                  TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
                  if ( TestBit != 0 )
                    BitMask[WordAddr] = (BitMask[WordAddr] ^ BitArray[BitAddr]);
//                  fprintf(log_fptr, "%4"ISYM"  %4u  %4u  %16o %16o\n", i, WordAddr, BitAddr, BitMask[WordAddr], BitArray[BitAddr]);
                }

              }
 
            }  // End of loop over particle position slices
 
          }  // End of loop over particle position dimensions
 
 
//        for (i=0; i<BitMaskSize; i++)
//          fprintf(log_fptr, "%4"ISYM"  %16o\n", i, BitMask[i]);
//
//        icol(Mask, TotalParticleCount, 128, log_fptr);
 
          NumberOfParticles = 0;
 
          for (i = 0; i < TotalParticleCount; i++) {

            XMask = 1;
            MaskAddr = i;
            WordAddr = MaskAddr/BitsPerInt;
            BitAddr  = MaskAddr%BitsPerInt;
            TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
            if ( TestBit == 0 )
              XMask = 0;
            NumberOfParticles += XMask;
 
//          fprintf(log_fptr, "%4"ISYM"    %2"ISYM" %2"ISYM"    %16o  %16o  %16o\n",
//                  i, Mask[i], XMask[i], BitMask[WordAddr], BitArray[BitAddr], TestBit);
 
          }
 
          // Compute the number of particles left
 
//        NumberOfParticles = 0;
//        for (i = 0; i < TotalParticleCount; i++)
//	    NumberOfParticles += Mask[i];
 
          if (io_log) fprintf(log_fptr, "NumberOfParticles %"ISYM"\n", NumberOfParticles);
 
 
 
 
          // Allocate space for the particles
 
          this->AllocateNewParticles(NumberOfParticles);
 
          // Read particle positions using mask and delete temp positions
 
          for (dim = 0; dim < GridRank; dim++) {
 
            n = 0; // counter for particles in grid
 
            for (pp1 = 0; pp1 < TotalParticleCount; pp1 = pp1 + ppbuffsize)
            {
              pp2 = min(pp1+ppbuffsize-1,TotalParticleCount-1);
              ppoffset = pp1;
              ppcount = pp2-pp1+1;
 
              Part = dim;
 
              if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", CosmologySimulationParticlePositionName);
 
              file_id = H5Fopen(CosmologySimulationParticlePositionName, H5F_ACC_RDONLY, H5P_DEFAULT);
                if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
                if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", CosmologySimulationParticlePositionName);
 
              dset_id = H5Dopen(file_id, CosmologySimulationParticlePositionName);
                if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
                if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              slab_offset[0] = Part;
              slab_stride[0] = 1;
              slab_count[0] = 1;
              slab_block[0] = 1;
 
              xfer_size = ppcount;
 
              slab_offset[1] = ppoffset;
              slab_stride[1] = 1;
              slab_count[1] = ppcount;
              slab_block[1] = 1;
 
              for ( idim=0; idim < Slab_Rank; idim++)
              {
                if (io_log) fprintf(log_fptr, "Slab [%"ISYM"] Offset %"ISYM", Count %"ISYM", Stride %"ISYM", Block %"ISYM"\n",
                        idim, (int) slab_offset[idim], (int) slab_count[idim],
                         (int) slab_stride[idim], (int) slab_block[idim]);
              }

              if (io_log) fprintf(log_fptr, "Xfer_size %"ISYM"\n", (int) xfer_size);
              if (io_log) fprintf(log_fptr, "Comp_size %"ISYM"\n", (int) component_size_attr);
 
              // Data in memory is considered 1D, stride 1, with zero offset
 
              mem_stride = 1;
              mem_count = xfer_size;
              mem_offset = 0;
              mem_block = 1;
 
              // 1D memory model
 
              mem_dsp_id = H5Screate_simple((Eint32) 1, &xfer_size, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
                if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, &mem_block);
                if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
                if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
                if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) TempField);
                if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", (int) h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(mem_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(file_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Dclose(dset_id);
                if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Fclose(file_id);
                if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              for (i = pp1, j = 0; i < pp1 + ppcount; i++, j++)
              {
                XMask = TRUE;
                MaskAddr = i;
                WordAddr = MaskAddr/BitsPerInt;
                BitAddr  = MaskAddr%BitsPerInt;
                TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
 
                if ( TestBit == 0 )
                  XMask = FALSE;
 
                if (XMask)
	          ParticlePosition[dim][n++] = TempField[j];
              }
 
            }  // End of loop over particle position slices
 
          }  // End of loop over particle position dimensions
 
 
 
 
          // Read particle velocities using mask
 
          if (CosmologySimulationParticleVelocityName != NULL)
          {
 
            for (dim = 0; dim < GridRank; dim++) {
 
              n = 0; // counter for particles in grid
 
              for (pp1 = 0; pp1 < TotalParticleCount; pp1 = pp1 + ppbuffsize)
              {
                pp2 = min(pp1+ppbuffsize-1,TotalParticleCount-1);
                ppoffset = pp1;
                ppcount = pp2-pp1+1;
 
                Part = dim;
 
                if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", CosmologySimulationParticleVelocityName);
 
                file_id = H5Fopen(CosmologySimulationParticleVelocityName, H5F_ACC_RDONLY, H5P_DEFAULT);
                  if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
                  if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
                if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", CosmologySimulationParticleVelocityName);
 
                dset_id = H5Dopen(file_id, CosmologySimulationParticleVelocityName);
                  if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
                  if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
 
                if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Rank\n");
 
                attr_id = H5Aopen_name(dset_id, "Component_Rank");
                  if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                  if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Aread(attr_id, HDF5_INT, &component_rank_attr);
                  if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Aclose(attr_id);
                  if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                if (io_log) fprintf(log_fptr, "COMPONENT_RANK %"ISYM"\n", component_rank_attr);
 
 
                if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Size\n");
 
                attr_id = H5Aopen_name(dset_id, "Component_Size");
                  if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                  if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Aread(attr_id, HDF5_INT, &component_size_attr);
                  if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Aclose(attr_id);
                  if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                if (io_log) fprintf(log_fptr, "COMPONENT_SIZE %"ISYM"\n", component_size_attr);
 
 
                if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Rank\n");
 
                attr_id = H5Aopen_name(dset_id, "Rank");
                  if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                  if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Aread(attr_id, HDF5_INT, &field_rank_attr);
                  if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Aclose(attr_id);
                  if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                if (io_log) fprintf(log_fptr, "RANK %"ISYM"\n", field_rank_attr);
 
 
                if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Dimensions\n");
 
                attr_id = H5Aopen_name(dset_id, "Dimensions");
                  if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                  if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
                attr_count = field_rank_attr;
 
                attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
                  if (io_log) fprintf(log_fptr, "Attr_dsp_id: %"ISYM"\n", attr_dsp_id);
                  if( attr_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Aread(attr_id, HDF5_INT, field_dims_attr);
                  if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Aclose(attr_id);
                  if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Sclose(attr_dsp_id);
                  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
 
                for (idim = 0; idim < field_rank_attr; idim++)
                {
                  if (io_log) fprintf(log_fptr, "DIMS %"ISYM":  %"ISYM"\n", idim, (int) field_dims_attr[idim]);
                }
 
                xfer_size = 1;
 
                for ( idim = 0; idim < field_rank_attr; idim++ )
                {
                  xfer_size = xfer_size * field_dims_attr[idim];
                }
 
                if (io_log) fprintf(log_fptr, "  Grid Elements %"ISYM"\n", (int) xfer_size);
 
                Slab_Rank = field_rank_attr+1;
 
                Slab_Dims[0] = component_rank_attr;
                Slab_Dims[1] = field_dims_attr[0];

                if (io_log) fprintf(log_fptr, "  Extended Rank %"ISYM"\n", (int) Slab_Rank);
                for ( idim = 0; idim < Slab_Rank; idim++ )
                {
                  if (io_log) fprintf(log_fptr, "    %"ISYM":  %"ISYM"\n", idim, (int) Slab_Dims[idim]);
                }
 
                slab_offset[0] = Part;
                slab_stride[0] = 1;
                slab_count[0] = 1;
                slab_block[0] = 1;
 
                xfer_size = ppcount;
 
                slab_offset[1] = ppoffset;
                slab_stride[1] = 1;
                slab_count[1] = ppcount;
                slab_block[1] = 1;
 
                for ( idim=0; idim < Slab_Rank; idim++)
                {
                  if (io_log) fprintf(log_fptr, "Slab [%"ISYM"] Offset %"ISYM", Count %"ISYM", Stride %"ISYM", Block %"ISYM"\n",
                      idim, (int) slab_offset[idim], (int) slab_count[idim],
                       (int) slab_stride[idim], (int) slab_block[idim]);
                }

                if (io_log) fprintf(log_fptr, "Xfer_size %"ISYM"\n", (int) xfer_size);
                if (io_log) fprintf(log_fptr, "Comp_size %"ISYM"\n", (int) component_size_attr);
 
                // Data in memory is considered 1D, stride 1, with zero offset
 
                mem_stride = 1;
                mem_count = xfer_size;
                mem_offset = 0;
                mem_block = 1;
 
                // 1D memory model
 
                mem_dsp_id = H5Screate_simple((Eint32) 1, &xfer_size, NULL);
                  if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
                  if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, &mem_block);
                  if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
                  if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
                  if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
                  if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) TempField);
                  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", (int) h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Sclose(mem_dsp_id);
                  if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Sclose(file_dsp_id);
                  if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Dclose(dset_id);
                  if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                h5_status = H5Fclose(file_id);
                  if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
                  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
                for (i = pp1, j = 0; i < pp1 + ppcount; i++, j++)
                {
                  XMask = TRUE;
                  MaskAddr = i;
                  WordAddr = MaskAddr/BitsPerInt;
                  BitAddr  = MaskAddr%BitsPerInt;
                  TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
 
                  if ( TestBit == 0 )
                    XMask = FALSE;
 
                  if (XMask)
                    ParticleVelocity[dim][n++] = TempField[j];
                }
 
              }  // End of loop over particle velocity slices
 
            }  // End of loop over particle velocity dimensions
 
          }
 
 
 
 
          // Read particle mass
 
          if (CosmologySimulationParticleMassName != NULL) {
 
            n = 0; // counter for particles in grid
 
            for (pp1 = 0; pp1 < TotalParticleCount; pp1 = pp1 + ppbuffsize)
            {
              pp2 = min(pp1+ppbuffsize-1,TotalParticleCount-1);
              ppoffset = pp1;
              ppcount = pp2-pp1+1;
 
              Part = 0;
 
              if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", CosmologySimulationParticleMassName);
 
              file_id = H5Fopen(CosmologySimulationParticleMassName, H5F_ACC_RDONLY, H5P_DEFAULT);
                if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
                if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", CosmologySimulationParticleMassName);
 
              dset_id = H5Dopen(file_id, CosmologySimulationParticleMassName);
                if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
                if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
 
              if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Rank\n");
 
              attr_id = H5Aopen_name(dset_id, "Component_Rank");
                if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aread(attr_id, HDF5_INT, &component_rank_attr);
                if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aclose(attr_id);
                if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              if (io_log) fprintf(log_fptr, "COMPONENT_RANK %"ISYM"\n", component_rank_attr);

 
              if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Size\n");
 
              attr_id = H5Aopen_name(dset_id, "Component_Size");
                if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aread(attr_id, HDF5_INT, &component_size_attr);
                if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aclose(attr_id);
                if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              if (io_log) fprintf(log_fptr, "COMPONENT_SIZE %"ISYM"\n", component_size_attr);

 
              if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Rank\n");
 
              attr_id = H5Aopen_name(dset_id, "Rank");
                if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aread(attr_id, HDF5_INT, &field_rank_attr);
                if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aclose(attr_id);
                if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              if (io_log) fprintf(log_fptr, "RANK %"ISYM"\n", field_rank_attr);
 
 
              if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Dimensions\n");
 
              attr_id = H5Aopen_name(dset_id, "Dimensions");
                if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
                if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              attr_count = field_rank_attr;
 
              attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
                if (io_log) fprintf(log_fptr, "Attr_dsp_id: %"ISYM"\n", attr_dsp_id);
                if( attr_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aread(attr_id, HDF5_INT, field_dims_attr);
                if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Aclose(attr_id);
                if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              for (idim = 0; idim < field_rank_attr; idim++)
              {
                if (io_log) fprintf(log_fptr, "DIMS %"ISYM":  %"ISYM"\n", idim, (int) field_dims_attr[idim]);
              }
 
              xfer_size = 1;
 
              for ( idim = 0; idim < field_rank_attr; idim++ )
              {
                xfer_size = xfer_size * field_dims_attr[idim];
              }
 
              if (io_log) fprintf(log_fptr, "  Grid Elements %"ISYM"\n", (int) xfer_size);
 
              Slab_Rank = field_rank_attr+1;
 
              Slab_Dims[0] = component_rank_attr;
              Slab_Dims[1] = field_dims_attr[0];
 
              if (io_log) fprintf(log_fptr, "  Extended Rank %"ISYM"\n", (int) Slab_Rank);
              for ( idim = 0; idim < Slab_Rank; idim++ )
              {
                if (io_log) fprintf(log_fptr, "    %"ISYM":  %"ISYM"\n", idim, (int) Slab_Dims[idim]);
              }
 
              slab_offset[0] = Part;
              slab_stride[0] = 1;
              slab_count[0] = 1;
              slab_block[0] = 1;
 
              xfer_size = ppcount;
 
              slab_offset[1] = ppoffset;
              slab_stride[1] = 1;
              slab_count[1] = ppcount;
              slab_block[1] = 1;
 
              for ( idim=0; idim < Slab_Rank; idim++)
              {
                if (io_log) fprintf(log_fptr, "Slab [%"ISYM"] Offset %"ISYM", Count %"ISYM", Stride %"ISYM", Block %"ISYM"\n",
                        idim, (int) slab_offset[idim], (int) slab_count[idim],
                         (int) slab_stride[idim], (int) slab_block[idim]);
              }

              if (io_log) fprintf(log_fptr, "Xfer_size %"ISYM"\n", (int) xfer_size);
              if (io_log) fprintf(log_fptr, "Comp_size %"ISYM"\n", (int) component_size_attr);
 
              // Data in memory is considered 1D, stride 1, with zero offset
 
              mem_stride = 1;
              mem_count = xfer_size;
              mem_offset = 0;
              mem_block = 1;
 
              // 1D memory model
 
              mem_dsp_id = H5Screate_simple((Eint32) 1, &xfer_size, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
                if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, &mem_block);
                if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
                if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
                if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
                if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) TempField);
                if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", (int) h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(mem_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Sclose(file_dsp_id);
                if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Dclose(dset_id);
                if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              h5_status = H5Fclose(file_id);
                if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
                if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
              for (i = pp1, j = 0; i < pp1 + ppcount; i++, j++)
              {
                XMask = TRUE;
                MaskAddr = i;
                WordAddr = MaskAddr/BitsPerInt;
                BitAddr  = MaskAddr%BitsPerInt;
                TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
 
                if ( TestBit == 0 )
                  XMask = FALSE;
 
                if (XMask)
                  ParticleMass[n++] = TempField[j];
              }
 
            }  // End of loop over particle mass slices
 
          }  // If there are particle masses
 
 
 
 
          // Read particle type
 
          if (CosmologySimulationParticleTypeName != NULL) {
 
          n = 0; // counter for particles in grid
 
          for (pp1 = 0; pp1 < TotalParticleCount; pp1 = pp1 + ppbuffsize)
          {
            pp2 = min(pp1+ppbuffsize-1,TotalParticleCount-1);
            ppoffset = pp1;
            ppcount = pp2-pp1+1;
 
            Part = 0;
 
            if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", CosmologySimulationParticleTypeName);
 
            file_id = H5Fopen(CosmologySimulationParticleTypeName, H5F_ACC_RDONLY, H5P_DEFAULT);
              if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
              if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", CosmologySimulationParticleTypeName);
 
            dset_id = H5Dopen(file_id, CosmologySimulationParticleTypeName);
              if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
              if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
 
            if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Rank\n");
 
            attr_id = H5Aopen_name(dset_id, "Component_Rank");
              if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
              if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Aread(attr_id, HDF5_INT, &component_rank_attr);
              if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Aclose(attr_id);
              if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            if (io_log) fprintf(log_fptr, "COMPONENT_RANK %"ISYM"\n", component_rank_attr);
 
 
            if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Size\n");
 
            attr_id = H5Aopen_name(dset_id, "Component_Size");
              if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
              if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Aread(attr_id, HDF5_INT, &component_size_attr);
              if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Aclose(attr_id);
              if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            if (io_log) fprintf(log_fptr, "COMPONENT_SIZE %"ISYM"\n", component_size_attr);
 
 
            if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Rank\n");
 
            attr_id = H5Aopen_name(dset_id, "Rank");
              if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
              if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Aread(attr_id, HDF5_INT, &field_rank_attr);
              if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Aclose(attr_id);
              if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            if (io_log) fprintf(log_fptr, "RANK %"ISYM"\n", field_rank_attr);
 
 
            if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Dimensions\n");
 
            attr_id = H5Aopen_name(dset_id, "Dimensions");
              if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
              if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            attr_count = field_rank_attr;
 
            attr_dsp_id = H5Screate_simple((Eint32) 1, &attr_count, NULL);
              if (io_log) fprintf(log_fptr, "Attr_dsp_id: %"ISYM"\n", attr_dsp_id);
              if( attr_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Aread(attr_id, HDF5_INT, field_dims_attr);
              if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Aclose(attr_id);
              if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            for (idim = 0; idim < field_rank_attr; idim++)
            {
              if (io_log) fprintf(log_fptr, "DIMS %"ISYM":  %"ISYM"\n", idim, (int) field_dims_attr[idim]);
            }
 
            xfer_size = 1;
 
            for ( idim = 0; idim < field_rank_attr; idim++ )
            {
              xfer_size = xfer_size * field_dims_attr[idim];
            }
 
            if (io_log) fprintf(log_fptr, "  Grid Elements %"ISYM"\n", (int) xfer_size);
 
            Slab_Rank = field_rank_attr+1;
 
            Slab_Dims[0] = component_rank_attr;
            Slab_Dims[1] = field_dims_attr[0];
 
            if (io_log) fprintf(log_fptr, "  Extended Rank %"ISYM"\n", (int) Slab_Rank);
            for ( idim = 0; idim < Slab_Rank; idim++ )
            {
              if (io_log) fprintf(log_fptr, "    %"ISYM":  %"ISYM"\n", idim, (int) Slab_Dims[idim]);
            }
 
            slab_offset[0] = Part;
            slab_stride[0] = 1;
            slab_count[0] = 1;
            slab_block[0] = 1;
 
            xfer_size = ppcount;
 
            slab_offset[1] = ppoffset;
            slab_stride[1] = 1;
            slab_count[1] = ppcount;
            slab_block[1] = 1;
 
            for ( idim=0; idim < Slab_Rank; idim++)
            {
              if (io_log) fprintf(log_fptr, "Slab [%"ISYM"] Offset %"ISYM", Count %"ISYM", Stride %"ISYM", Block %"ISYM"\n",
                  idim, (int) slab_offset[idim], (int) slab_count[idim],
                           (int) slab_stride[idim], (int) slab_block[idim]);
            }

            if (io_log) fprintf(log_fptr, "Xfer_size %"ISYM"\n", (int) xfer_size);
            if (io_log) fprintf(log_fptr, "Comp_size %"ISYM"\n", (int) component_size_attr);
 
            // Data in memory is considered 1D, stride 1, with zero offset
 
            mem_stride = 1;
            mem_count = xfer_size;
            mem_offset = 0;
            mem_block = 1;
 
            // 1D memory model
 
            mem_dsp_id = H5Screate_simple((Eint32) 1, &xfer_size, NULL);
              if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
              if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, &mem_block);
              if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
              if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
              if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, slab_block);
              if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Dread(dset_id, int_type_id, mem_dsp_id, file_dsp_id, H5P_DEFAULT, (VOIDP) IntTempField);
              if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", (int) h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sclose(mem_dsp_id);
              if (io_log) fprintf(log_fptr, "H5Sclose mem_dsp: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Sclose(file_dsp_id);
              if (io_log) fprintf(log_fptr, "H5Sclose file_dsp: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Dclose(dset_id);
              if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            h5_status = H5Fclose(file_id);
              if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
              if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
            for (i = pp1, j = 0; i < pp1 + ppcount; i++, j++)
            {
              XMask = TRUE;
              MaskAddr = i;
              WordAddr = MaskAddr/BitsPerInt;
              BitAddr  = MaskAddr%BitsPerInt;
              TestBit = (BitMask[WordAddr] & BitArray[BitAddr]);
 
              if ( TestBit == 0 )
                XMask = FALSE;
 
              if (XMask)
                ParticleType[n++] = IntTempField[j];
            }
 
          }  // End of loop over particle type slices
 
        }  // If there is particle type

 
        // Delete mask and temp field
 
        if (io_log) fprintf(log_fptr, "De-Allocate Mask and TempField\n");
 
        delete [] TempField;
        delete [] IntTempField;
 
        delete [] BitMask;
        delete [] BitArray;
 
      }  // End of normal ||rgio (Particles NOT pre-sorted)
 
 
 
 
    } else {
 
 
 
 
      // For regular IO (non-parallel)

      if (debug) 
        printf("Standard Serial I/O\n");
 
      if (io_log) fprintf(log_fptr, "CSIG Serial I/O\n");
 
      // Set Number Of particles
 
      NumberOfParticles = TempIntArray[0];
 
      End[0] = NumberOfParticles - 1;
      Dim[0] = NumberOfParticles;
 
      // Allocate space for the particles
 
      this->AllocateNewParticles(NumberOfParticles);
 
      // Read particle positions

      if (CosmologySimulationParticlePositionName != NULL) { 
        for (dim = 0; dim < GridRank; dim++) {
          if (ReadFile(CosmologySimulationParticlePositionName, 1, Dim,
                       Start, End, Zero, NULL, &tempbuffer, dim, 3) == FAIL) {
            ENZO_VFAIL("Error reading particle position %"ISYM".\n", dim)
          }
          for (i = Start[0]; i <= End[0]; i++) {
            ParticlePosition[dim][i] = FLOAT(tempbuffer[i]);
          }
 
          if (io_log) fprintf(log_fptr, "De-Allocate [] tempbuffer for dim = %"ISYM"\n", dim);
 
          delete [] tempbuffer;  // WHY???

        }
      }
 
      // Read particle velocities
 
      if (CosmologySimulationParticleVelocityName != NULL) {
	for (dim = 0; dim < GridRank; dim++) {
	  if (ReadFile(CosmologySimulationParticleVelocityName, 1, Dim,
	     Start, End, Zero, ParticleVelocity[dim], &tempbuffer, dim, 3) == FAIL) {
	    ENZO_VFAIL("Error reading particle velocity %"ISYM".\n", dim)
	  }
	}
      }
 
      // Read particle mass
 
      if (CosmologySimulationParticleMassName != NULL) {
	if (ReadFile(CosmologySimulationParticleMassName, 1, Dim, Start, End,
		     Zero, ParticleMass, &tempbuffer, 0, 1) == FAIL) {
	  ENZO_FAIL("Error reading particle masses.\n");
	}
      }
 
      // Read particle type
 
      if (CosmologySimulationParticleTypeName != NULL) {
	if (ReadIntFile(CosmologySimulationParticleTypeName, 1, Dim, Start, End,
		     Zero, ParticleType, &int_tempbuffer, 0, 1) == FAIL) {
	  ENZO_FAIL("Error reading particle types.\n");
	}
      }
 
 
    } // end: if ParallelRootGridIO == TRUE
 
 
 
 
    // If there are particles, but no velocity file, then set the velocities to zero
 
    if (NumberOfParticles > 0 && CosmologySimulationParticleVelocityName == NULL) {
      ENZO_FAIL("Error -- no velocity field specified.\n");
      //  printf("CosmologySimulation warning: setting velocities to zero.\n");
      //      for (dim = 0; dim < GridRank; dim++)
      //	for (i = 0; i < NumberOfParticles; i++)
      //	  ParticleVelocity[dim][i] = 0.0;
    }
 
    // If there are particles, but a mass file name wasn't specified, make
    // up the particles mass
 
    if (NumberOfParticles > 0 && CosmologySimulationParticleMassName == NULL) {
 
      // Compute mass of each particle
 
      float UniformParticleMass = CosmologySimulationOmegaCDMNow/
	                          OmegaMatterNow;
 
      // If there are exactly 1/8 as many particles as cells, then set the
      // particle mass to 8 times the usual
 
      int NumberOfActiveCells = (GridEndIndex[0]-GridStartIndex[0]+1)*
	                        (GridEndIndex[1]-GridStartIndex[1]+1)*
	                        (GridEndIndex[2]-GridStartIndex[2]+1);
      if (NumberOfParticles*8 == NumberOfActiveCells)
	UniformParticleMass *= 8;
      if (NumberOfParticles == NumberOfActiveCells*8)
        UniformParticleMass /= 8;
 
      // UniformParticleMass *= float(POW(TotalRefinement, GridRank));
	/*   for (dim = 0; dim < GridRank; dim++)
               UniformParticleMass *= float(GridEndIndex[dim]-GridStartIndex[dim]+1);
               UniformParticleMass /= float(NumberOfParticles); */
 
      // Set all particle masses to this value
 
      if (debug)
        printf("CosmologySimulationInitializeGrid: particle mass = %"GSYM"\n", UniformParticleMass);

      for (i = 0; i < NumberOfParticles; i++)
	ParticleMass[i] = UniformParticleMass;
 
    }
 
    // Set Particle attributes to FLOAT_UNDEFINED
 
    for (j = 0; j < NumberOfParticleAttributes; j++)
      for (i = 0; i < NumberOfParticles; i++)
	ParticleAttribute[j][i] = FLOAT_UNDEFINED;
 
    // If no type file was read, then set all particles to dm particles
 
    if (CosmologySimulationParticleTypeName == NULL)
      for (i = 0; i < NumberOfParticles; i++)
        ParticleType[i] = PARTICLE_TYPE_DARK_MATTER;
 
    // Assign particles unique identifiers
 
    if (ParallelRootGridIO == TRUE)
      CurrentParticleNumber = 0;  // Unique ID's will be added later
 
    for (i = 0; i < NumberOfParticles; i++)
      ParticleNumber[i] = CurrentParticleNumber + i;
 
    CurrentParticleNumber += NumberOfParticles;



 
    } // end: if (CosmologySimulationParticleName != NULL && ReadData)

 
  } // end: if (ProcessorNumber == MyProcessorNumber)




  // All processors particpate from here on
 
  // Share number of particles amoung processors
 
  //  The following causes I/O serialization and should be bypassed
  //  for ParallelRootGridIO = 1
 
  if ( ParallelRootGridIO != TRUE ) {
    if (ReadData == TRUE) {
      CommunicationBroadcastValue(&NumberOfParticles, ProcessorNumber);
    }
  }
 
  OldTime = Time;
 
  if (io_log) fclose(log_fptr);

 
  return SUCCESS;
}
