/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COSMOLOGY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Robert Harkness, July 2002
/              Changed default initialization of Total_Energy
/              The divisors DEFAULT_MU and Gamma-1.0 were missing,
/              leading to an effectively incorrect initial T.
/  modifiedN:  Robert Harkness, December 2003
/              Nested grids, 64-bit Inits
/  modified2:  Robert Harkness, March 2006
/              Use MPICH V1.x Dims_create code
/  modified3:  Robert Harkness, April 2008
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif
 
#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

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
#define DENSITY_FLOOR 0.01
 
#define READFILE ReadFile
 
// Function prototypes
 
int ReadFile(char *name, int Rank, int Dims[], int StartIndex[],
	     int EndIndex[], int BufferOffset[], float *buffer,
	     inits_type **tempbuffer, int Part, int Npart);
 
int ReadIntFile(char *name, int Rank, int Dims[], int StartIndex[],
		int EndIndex[], int BufferOffset[], int *buffer,
		int **tempbuffer, int Part, int Npart);

void ReadAttribute(hid_t dset_id, int *Attribute, char *AttributeName, FILE *log_fptr, int io_log);
int ReadAttr(char *Fname, int *Rank, int Dims[], int *NSeg, int *LSeg, FILE *log_fptr);
 
void fcol(float *x, int n, int m, FILE *log_fptr);
void pcol(FLOAT *x, int n, int m, FILE *log_fptr);
void icol(int *x, int n, int m, FILE *log_fptr);
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
 
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
 
int Enzo_Dims_create(int nnodes, int ndims, int *dims); 


int grid::NestedCosmologySimulationInitializeGrid(
                          int   InitialGridNumber,
                          float CosmologySimulationOmegaBaryonNow,
                          float CosmologySimulationOmegaCDMNow,
                          float CosmologySimulationInitialTemperature,
                          char *CosmologySimulationDensityName,
                          char *CosmologySimulationTotalEnergyName,
                          char *CosmologySimulationGasEnergyName,
                          char *CosmologySimulationVelocityNames[],
                          char *CosmologySimulationParticlePositionName,
                          char *CosmologySimulationParticleVelocityName,
 			  char *CosmologySimulationParticleDisplacementName,
                          char *CosmologySimulationParticleMassName,
                          char *CosmologySimulationParticleTypeName,
                          char *CosmologySimulationParticleVelocityNames[],
 			  char *CosmologySimulationParticleDisplacementNames[],
                          int   CosmologySimulationSubgridsAreStatic,
                          int   TotalRefinement,
                          float CosmologySimulationInitialFractionHII,
                          float CosmologySimulationInitialFractionHeII,
                          float CosmologySimulationInitialFractionHeIII,
                          float CosmologySimulationInitialFractionHM,
                          float CosmologySimulationInitialFractionH2I,
                          float CosmologySimulationInitialFractionH2II,
			  float CosmologySimulationInitialFractionMetal,
			  float CosmologySimulationInitialFractionMetalIa,
                          int   UseMetallicityField,
                          PINT &CurrentParticleNumber,
                          int CosmologySimulationManuallySetParticleMassRatio,
                          float CosmologySimulationManualParticleMassRatio,
                          int CosmologySimulationCalculatePositions,
			  FLOAT SubDomainLeftEdge[],
			  FLOAT SubDomainRightEdge[],
			  float CosmologySimulationInitialUniformBField[])

{ 
 
  // Declarations
 
  int idim, ndim, dim, i, j, vel, OneComponentPerFile, level;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum, MetalIaNum;
 
  int iTE = ietot;
  int ExtraField[2];
  int ForbidNum;
  int MachNum, PSTempNum, PSDenNum;
 
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
 
  hsize_t    mem_offset;
  hsize_t    slab_offset[4];
 
  int NSeg, LSeg, Part;
 
  int component_rank_attr;
  int component_size_attr;
  int field_rank_attr;
  int field_dims_attr[3];
 
  int TopGridStart[3];
  int TopGridEnd[3];
  int TopGridDims[3];
 
  //  float SubDomainLeftEdge[3];
  //  float SubDomainRightEdge[3];
  float SubCellWidth[3];
 
  int enzo_layout[3];
  int A,B,CC,I,J,K,M;

  //  MPI_Arg NOP,RNK;
  //  MPI_Arg mpi_layout[3];

  FILE *log_fptr;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

  // Warn if number of processors is too large for serial IO
 
  if ( NumberOfProcessors > 64 )
    {
      if ( (ParallelRootGridIO != TRUE) && (ParallelParticleIO != TRUE) )
	{
	  if (MyProcessorNumber == ROOT_PROCESSOR)
        fprintf(stderr, "WARNING: ParallelRootGridIO and ParallelParticleIO are recommended for > 64 cpus to shrink data read-in time.");
	}
    }
 
  char pid[MAX_TASK_TAG_SIZE];
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  char *logname = new char[MAX_NAME_LENGTH];
  strcpy(logname, "CSlog.");
  strcat(logname, pid);
 
  if (io_log) log_fptr = fopen(logname, "a");
 
  delete [] logname;
 
  if (io_log) fprintf(log_fptr, "\n");
  if (io_log) fprintf(log_fptr, "Processor %"ISYM", Target processor %"ISYM"\n", MyProcessorNumber, ProcessorNumber);
 
  // Determine if the data should be loaded in or not
  // TotalRefinement == 1 indicated the top grid
  // In the current version, nested grids below top are read on the
  // first pass using serial I/O.  Only top uses ||rgio and needs a
  // subsequent call to CosmologySimulationReInitialize which sets
  // TotalRefinement = -1 to signal a "real" read
 
  int ReadData = TRUE, Offset[] = {0,0,0};
 
  //  if (ParallelRootGridIO == TRUE && TotalRefinement == 1)
  //    ReadData = FALSE;
 
  if (ParallelRootGridIO == TRUE && TotalRefinement > 0)
    ReadData = FALSE;
 
  if (io_log) fprintf(log_fptr, "CSIG ParallelRootGridIO = %"ISYM"\n", ParallelRootGridIO);
  if (io_log) fprintf(log_fptr, "TotalRefinement = %"ISYM"\n", TotalRefinement);
  if (io_log) fprintf(log_fptr, "ReadData = %"ISYM"\n", ReadData);
 
  /* Calculate buffer Offset (same as Grid unless doing ParallelRootGridIO
     (TotalRefinement = -1 if used as a signal that we should really load
     in the data regardless of the value of ParallelRootGridIO). */
 
  // The offsets need to be calculated at each level of initial grid so the default
  // expression for position within the topgrid won't do...
 
  //  if (ParallelRootGridIO == TRUE && TotalRefinement == -1)
  //    for (dim = 0; dim < GridRank; dim++)
  //      Offset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/CellWidth[dim][0]);
 
  if (ParallelRootGridIO == TRUE && TotalRefinement < 0) {

    if (CosmologySimulationDensityName != NULL) {
      file_id = H5Fopen(CosmologySimulationDensityName, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
      if( file_id == h5_error){ENZO_FAIL("IO Problem");}
 
      dset_id = H5Dopen(file_id, CosmologySimulationDensityName);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error){ENZO_FAIL("IO Problem");}
    } else if (CosmologySimulationParticleVelocityName != NULL) {
      file_id = H5Fopen(CosmologySimulationParticleVelocityName, 
			H5F_ACC_RDONLY, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
      if (file_id == h5_error){ENZO_FAIL("IO Problem");}
      dset_id = H5Dopen(file_id, CosmologySimulationParticleVelocityName);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if (dset_id == h5_error){ENZO_FAIL("IO Problem");}
    } else if (CosmologySimulationParticleVelocityNames[0] != NULL) {
      file_id = H5Fopen(CosmologySimulationParticleVelocityNames[0], 
			H5F_ACC_RDONLY, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
      if (file_id == h5_error){ENZO_FAIL("IO Problem");}
      dset_id = H5Dopen(file_id, CosmologySimulationParticleVelocityNames[0]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if (dset_id == h5_error){ENZO_FAIL("IO Problem");}
    } else
      ENZO_FAIL("Cannot find initial density or particle velocity datafile"
		" for nested cosmology run.");

    

    /* Read Rank, Dimensions, TopGridStart, TopGridEnd and TopGridDims attributes */
 
    ReadAttribute(dset_id, &field_rank_attr, "Rank", log_fptr, io_log);
    if (io_log) fprintf(log_fptr, "RANK %"ISYM"\n", field_rank_attr);
 
    ReadAttribute(dset_id, field_dims_attr, "Dimensions", log_fptr, io_log);
 
    ReadAttribute(dset_id, TopGridStart, "TopGridStart", log_fptr, io_log);
    if (io_log) fprintf(log_fptr, "TopGridStart: %"ISYM" %"ISYM" %"ISYM"\n", TopGridStart[0], TopGridStart[1], TopGridStart[2]);

    ReadAttribute(dset_id, TopGridEnd, "TopGridEnd", log_fptr, io_log);
    if (io_log) fprintf(log_fptr, "TopGridEnd: %"ISYM" %"ISYM" %"ISYM"\n", TopGridEnd[0], TopGridEnd[1], TopGridEnd[2]);
 
    ReadAttribute(dset_id, TopGridDims, "TopGridDims", log_fptr, io_log);
    if (io_log) fprintf(log_fptr, "TopGridDims: %"ISYM" %"ISYM" %"ISYM"\n", TopGridDims[0], TopGridDims[1], TopGridDims[2]);

    /* Close dataset and file. */
 
    h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
    h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
    for (dim = 0; dim < GridRank; dim++) {
      if (TopGridStart[dim] == INT_UNDEFINED)
	TopGridStart[dim] = 0;
      if (TopGridEnd[dim] == INT_UNDEFINED)
	TopGridEnd[dim] = TopGridDims[dim] - 1;
      // Correct for a "feature" in mpgrafic
      if (CosmologySimulationCalculatePositions)
	TopGridEnd[dim] -= 1;
    }
 
    //        fprintf(stderr, "TopGridDims = %"ISYM" %"ISYM" %"ISYM"\n", TopGridDims[0], TopGridDims[1], TopGridDims[2]);
    //        fprintf(stderr, "TopGridStart: %"ISYM" %"ISYM" %"ISYM"\n", TopGridStart[0], TopGridStart[1], TopGridStart[2]);
    //        fprintf(stderr, "TopGridEnd: %"ISYM" %"ISYM" %"ISYM"\n", TopGridEnd[0], TopGridEnd[1], TopGridEnd[2]);

    /* Calculate the start of the entire subgrid region in problem units. */
    /* This is now passed in from above rather than calculated here because this
       assumes that initial refined regions on any level start and end of top grid
       cell boundaries which is not necessarily true. */
 
    for (dim = 0; dim < GridRank; dim++) {
      //      SubDomainLeftEdge[dim] = TopGridStart[dim] * (DomainRightEdge[dim]-DomainLeftEdge[dim])/((float) TopGridDims[dim]);
      //      SubDomainRightEdge[dim] = (TopGridEnd[dim]+1) * (DomainRightEdge[dim]-DomainLeftEdge[dim])/((float) TopGridDims[dim]);


      SubCellWidth[dim] = (SubDomainRightEdge[dim]-SubDomainLeftEdge[dim])/((float) field_dims_attr[dim]);
    }
 
    for (dim = 0; dim < GridRank; dim++)
      if (io_log) fprintf(log_fptr, "DDim %"ISYM" : LHS %10.4"FSYM"  RHS %10.4"FSYM"\n", dim, DomainLeftEdge[dim], DomainRightEdge[dim]);
 
    for (dim = 0; dim < GridRank; dim++)
      if (io_log) fprintf(log_fptr, "SDim %"ISYM" : LHS %10.4"FSYM"  RHS %10.4"FSYM"\n", dim, SubDomainLeftEdge[dim], SubDomainRightEdge[dim]);
 
    for (dim = 0; dim < GridRank; dim++)
      if (io_log) fprintf(log_fptr, "GDim %"ISYM" : LHS %10.4"FSYM"  RHS %10.4"FSYM"\n", dim, GridLeftEdge[dim], GridRightEdge[dim]);
 
  } // end: if (ParallelRootGridIO == TRUE && TotalRefinement < 0)
 
  if (ParallelRootGridIO == TRUE && TotalRefinement < 0)
    for (dim = 0; dim < GridRank; dim++){
      Offset[dim] = nint((GridLeftEdge[dim] - SubDomainLeftEdge[dim])/SubCellWidth[dim]);
    }
  /*
    if (ParallelRootGridIO == TRUE && TotalRefinement < 0)
    for (dim = 0; dim < GridRank; dim++)
    Offset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/CellWidth[dim][0]);
  */
  if (io_log) fprintf(log_fptr, "CSIG Offsets %"ISYM" %"ISYM" %"ISYM"\n", Offset[0], Offset[1], Offset[2]);
 
 
  level = (InitialGridNumber > 0) ? 
    StaticRefineRegionLevel[InitialGridNumber-1] + 1 : 0;
 
  // Create baryon fields (unless they are not needed)
 
  NumberOfBaryonFields = 0;
 
  if (CosmologySimulationDensityName != NULL) {
 
    //  fprintf(stderr, "Create baryon fields for %s on CPU %"ISYM"\n", CosmologySimulationDensityName, MyProcessorNumber);
 
    FieldType[NumberOfBaryonFields++] = Density;
  vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1 || (HydroMethod == MHD_RK) || (HydroMethod == HD_RK))
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2 || (HydroMethod == MHD_RK) || (HydroMethod == HD_RK))
    FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[iTE = NumberOfBaryonFields++] = TotalEnergy;

  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;

  if (HydroMethod == MHD_RK) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }
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
      if (StarMakerTypeIaSNe)
	FieldType[MetalIaNum = NumberOfBaryonFields++] = MetalSNIaDensity;
      if(MultiMetals){
	FieldType[ExtraField[0] = NumberOfBaryonFields++] = ExtraType0;
	FieldType[ExtraField[1] = NumberOfBaryonFields++] = ExtraType1;
      }
    }
    if (WritePotential)
      FieldType[NumberOfBaryonFields++] = GravPotential;
    if(STARMAKE_METHOD(COLORED_POP3_STAR)){
      fprintf(stderr, "Initializing Forbidden Refinement color field\n");
      FieldType[ForbidNum = NumberOfBaryonFields++] = ForbiddenRefinement;
    }
    if(ShockMethod){
      FieldType[MachNum   = NumberOfBaryonFields++] = Mach;
      if(StorePreShockFields){
	FieldType[PSTempNum = NumberOfBaryonFields++] = PreShockTemperature;
	FieldType[PSDenNum = NumberOfBaryonFields++] = PreShockDensity;
      }
    }    
  }


  //  fprintf(stderr, "Total Baryon Fields in VVV: %"ISYM" on CPU %"ISYM"\n", NumberOfBaryonFields, MyProcessorNumber);
 
  // Set the subgrid static flag
 
  SubgridsAreStatic = CosmologySimulationSubgridsAreStatic;
 
  // Return if this doesn't concern us
  // The grids are initilized on the owning processor only
  // When this is called by InitializeNew, it is only done for
  // the nested grids because ReadData != -1 on the top grid
  
  if (ProcessorNumber == MyProcessorNumber) {
 
    if (NumberOfBaryonFields > 0) {
 
      // Get the units so we can convert temperature later
 
      float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
	VelocityUnits=1;
      double MassUnits=1;
 
      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		   &TimeUnits, &VelocityUnits, &MassUnits, 
		   InitialTimeInCodeUnits) == FAIL) {
	ENZO_FAIL("Error in GetUnits.\n");
      }
 
      // Determine the size of the fields
 
      int size = 1;
 
      for (dim = 0; dim < GridRank; dim++)
	size *= GridDimension[dim];
 
      // Allocate space for the fields ONLY if ReadData is TRUE!
     
      if (ReadData == TRUE) {
	if (io_log) fprintf(log_fptr, "Allocate %"ISYM" fields, %"ISYM" floats per field\n", NumberOfBaryonFields, size);
	for (dim=0; dim < GridRank; dim++) {
	  if (io_log) fprintf(log_fptr, "  Field dim %"ISYM" size %"ISYM"\n", dim, GridDimension[dim]);
	}
      }
 
      if (ReadData == TRUE && debug)
	printf("Allocating %"ISYM" baryon fields of size %"ISYM"\n", NumberOfBaryonFields, size);
 
      if (ReadData == TRUE)
	for (int field = 0; field < NumberOfBaryonFields; field++)
	  BaryonField[field] = new float[size];
 
      // Read the density field
 
      if (CosmologySimulationDensityName != NULL && ReadData) {
	if (READFILE(CosmologySimulationDensityName, GridRank, GridDimension,
		     GridStartIndex, GridEndIndex, Offset, BaryonField[0],
		     &tempbuffer, 0, 1) == FAIL) {
	  ENZO_FAIL("Error reading density field.\n");
	}
	//  fcol(BaryonField[0], size, 10, log_fptr);
	for (i = 0; i < size; i++)
	  BaryonField[0][i] = max(BaryonField[0][i], DENSITY_FLOOR);
      }
 
      // Read the total energy field
 
      if (CosmologySimulationTotalEnergyName != NULL && ReadData)
	if (READFILE(CosmologySimulationTotalEnergyName, GridRank,
		     GridDimension, GridStartIndex, GridEndIndex, Offset,
		     BaryonField[iTE], &tempbuffer, 0, 1) == FAIL) {
	  ENZO_FAIL("Error reading total energy field.\n");
	}
 
      // Read the gas energy field
 
      if (CosmologySimulationGasEnergyName != NULL && DualEnergyFormalism && ReadData)
	if (READFILE(CosmologySimulationGasEnergyName, GridRank, GridDimension,
		     GridStartIndex, GridEndIndex, Offset, BaryonField[2],
		     &tempbuffer, 0, 1) == FAIL) {
	  ENZO_FAIL("Error reading gas energy field.\n");
	}
 
      // Read the velocity fields
 
      if (CosmologySimulationVelocityNames[0] != NULL && ReadData) {
	// Determine if we're reading different files for each component
	if (GridRank > 1)
	  if (strstr(CosmologySimulationVelocityNames[0], 
		     CosmologySimulationVelocityNames[1]) == NULL)
	    OneComponentPerFile = TRUE;
	  else
	    OneComponentPerFile = FALSE;

	for (dim = 0; dim < GridRank; dim++) {
	  if (OneComponentPerFile) {
	    ndim = 1;
	    idim = 0;
	  } else {
	    ndim = 3;
	    idim = dim;
	  }
	  if (READFILE(CosmologySimulationVelocityNames[dim], GridRank,
		       GridDimension, GridStartIndex, GridEndIndex, Offset,
		       BaryonField[vel+dim], &tempbuffer, idim, ndim) == FAIL) {
	    ENZO_VFAIL("Error reading velocity field %"ISYM".\n", dim)
	  }
	} // ENDFOR dim
      } // ENDIF grid velocities
      
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
 
	} // end: loop over i
      } // end: if (MultiSpecies && ReadData)
 
      // If using metallicity, set the field
      
      if (UseMetallicityField && ReadData) {
	for (i = 0; i < size; i++)
	  BaryonField[MetalNum][i] = CosmologySimulationInitialFractionMetal
	    * BaryonField[0][i];

	if (StarMakerTypeIaSNe)
	  for (i = 0; i < size; i++)
	    BaryonField[MetalIaNum][i] = CosmologySimulationInitialFractionMetalIa
	      * BaryonField[0][i];
	
	if (MultiMetals) {
	  for (i = 0; i < size; i++) {
	    BaryonField[ExtraField[0]][i] = CosmologySimulationInitialFractionMetal
	      * BaryonField[0][i];
	    BaryonField[ExtraField[1]][i] = CosmologySimulationInitialFractionMetal
	      * BaryonField[0][i];
	  }
	}

	if (STARMAKE_METHOD(COLORED_POP3_STAR) && ReadData) {
	  for (i = 0; i < size; i++)
	    BaryonField[ForbidNum][i] = 0.0;
	}
      } // ENDIF UseMetallicityField
 
      // If they were not read in above, set the total & gas energy fields now
 
      if (CosmologySimulationDensityName != NULL && ReadData) {

	if (StringKick > 0.) { // gives only baryons a uniform kick velocity in x direction
	  // models http://adsabs.harvard.edu/abs/2010PhRvD..82h3520T
	  printf("adding string kick %"FSYM" %"FSYM"\n", StringKick, 
		 StringKick/VelocityUnits*1e5);
	  int dim0 = vel + StringKickDimension;
	  int dim1 = vel + (StringKickDimension+1) % GridRank;
	  int dim2 = vel + (StringKickDimension+2) % GridRank;
	  for (i = 0; i < size; i++) {
	    BaryonField[0][i]   = 	    
	      (CosmologySimulationOmegaBaryonNow)/(OmegaMatterNow);
	    BaryonField[dim0][i] = StringKick/VelocityUnits*1e5; // input in km/s
	    BaryonField[dim1][i] = 0.; // do not neglect initial perturbations. (below jeans length)
	    BaryonField[dim2][i] = 0.;
	  }
	}

	if (CosmologySimulationTotalEnergyName == NULL)
	  for (i = 0; i < size; i++)
	    BaryonField[iTE][i] = CosmologySimulationInitialTemperature/
	      TemperatureUnits/DEFAULT_MU/(Gamma-1.0);
 
	/*          * POW(BaryonField[0][i]/CosmologySimulationOmegaBaryonNow,Gamma-1)
		    / (Gamma-1); */
 
	if (CosmologySimulationGasEnergyName == NULL && DualEnergyFormalism)
	  for (i = 0; i < size; i++)
	    BaryonField[iTE+1][i] = BaryonField[iTE][i];
 
	if (CosmologySimulationTotalEnergyName == NULL &&
	    HydroMethod != Zeus_Hydro) {
	  for (dim = 0; dim < GridRank; dim++)
	    for (i = 0; i < size; i++) {
	      BaryonField[iTE][i] +=
		0.5 * BaryonField[vel+dim][i] * BaryonField[vel+dim][i];
	      
	      if (HydroMethod == MHD_RK) {
		BaryonField[iBx  ][i] = CosmologySimulationInitialUniformBField[0];
		BaryonField[iBy  ][i] = CosmologySimulationInitialUniformBField[1];
		BaryonField[iBz  ][i] = CosmologySimulationInitialUniformBField[2];
		BaryonField[iPhi ][i] = 0.0;
		BaryonField[iTE][i] += 0.5*(BaryonField[iBx][i] * BaryonField[iBx][i]+
					    BaryonField[iBy][i] * BaryonField[iBy][i]+
					    BaryonField[iBz][i] * BaryonField[iBz][i])/
		  BaryonField[iden][i];
	      }
	    }
	}
      } // end: if (CosmologySimulationDensityName != NULL)

      // Shock/Cosmic Ray Model
      if (ShockMethod && ReadData) {
	for (i = 0; i < size; i++) {
	  BaryonField[MachNum][i] = tiny_number;
	  if (StorePreShockFields) {
	    BaryonField[PSTempNum][i] = tiny_number;
	    BaryonField[PSDenNum][i] = tiny_number;
	    }
	}
      } // end: if (ShockMethod && ReadData)

    } // end: if (NumberOfBaryonFields > 0)

    // ------------------------------------------------------------ 
    // Read particle fields if the was a name specified
 
    // This routine reads data from Inits and Ring: 32- or 64-bit
 
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
  
    if ((CosmologySimulationParticlePositionName != NULL ||
	 CosmologySimulationCalculatePositions) && ReadData) {


      // Get the total number of particles from this file by reading the file attributes
 
      int TempInt, Dim[1], Start[1] = {0}, End[1], Zero[1] = {0};
      int TempIntArray[MAX_DIMENSION], TotalParticleCount;
      if (!CosmologySimulationCalculatePositions) {

	ReadAttr(CosmologySimulationParticlePositionName,
		 &TempInt, TempIntArray, &NSeg, &LSeg, log_fptr);
 
	if (TempInt != 1) {
	  ENZO_VFAIL("Rank (%"ISYM") is not one in file %s.\n", TempInt,
		  CosmologySimulationParticlePositionName)
	}
 
	/* If doing parallel root grid IO then read in the full list of particle
	   positions and construct a mask of those within this grid. */
 
	TotalParticleCount = TempIntArray[0];
 
      } // end: if (!CosmologySimulationCalculatePositions)
 
      int PreSortedParticles = 0;
      if (ParallelParticleIO && !CosmologySimulationCalculatePositions)
	PreSortedParticles = 1;
      

      // ---------------------------------------------------
      //  This section only used with
      //  Pre-sorted particle input (ring i/o)

      if (ParallelRootGridIO == TRUE && PreSortedParticles == 1) {
 

	int NumSortedParticles = 0;
	int TotParticleCount = 0;
 
	//  Generate subgrid to MPI task map for filenames

	MPI_Arg mpi_size;
 
#ifdef USE_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else
	mpi_size = 1;
#endif

	int NOP, RNK;
	int mpi_layout[3];

	NOP = mpi_size;
	RNK = 3;

	mpi_layout[0] = 0;
	mpi_layout[1] = 0;
	mpi_layout[2] = 0;

	Enzo_Dims_create(NOP, RNK, mpi_layout);
 
	for (idim = 0; idim < GridRank; idim++)
	  enzo_layout[idim] = mpi_layout[GridRank-1-idim];

	if (MyProcessorNumber == ROOT_PROCESSOR) {
	  fprintf(stderr, "ENZO_layout %"ISYM" x %"ISYM" x %"ISYM"\n", enzo_layout[0], enzo_layout[1], enzo_layout[2]);
	}

	/* This calculates the I,J,K index of the particle region to read in
	   (must match calculation in ring). */
 
	A = enzo_layout[0];
	B = enzo_layout[1];
	CC = enzo_layout[2];
	I = nint( (GridLeftEdge[0]-SubDomainLeftEdge[0])/((SubDomainRightEdge[0]-SubDomainLeftEdge[0])/((float) A)));
	J = nint( (GridLeftEdge[1]-SubDomainLeftEdge[1])/((SubDomainRightEdge[1]-SubDomainLeftEdge[1])/((float) B)));
	K = nint( (GridLeftEdge[2]-SubDomainLeftEdge[2])/((SubDomainRightEdge[2]-SubDomainLeftEdge[2])/((float) CC)));
	//    M = ((I*B*CC) + J*CC) + K;
	M = ((K*B*A) + J*A) + I;
 
	if (io_log) fprintf(log_fptr, "ABC %"ISYM" %"ISYM" %"ISYM";  IJK %"ISYM" %"ISYM" %"ISYM";  M = %"ISYM"\n", A,B,CC,I,J,K,M);
 
	sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, M);
 
	printf("PreSortedParticles - ParallelRootGridIO\n");

	/* Generate names for ring particle files. */
 
	int fln, ext, lex;
	char *Extension;
 
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
	if( file_id == h5_error){ENZO_FAIL("IO Problem");}
 
	dset_id = H5Dopen(file_id, PPos);
	if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	if( dset_id == h5_error ){ENZO_FAIL("IO Problem");}

	ReadAttribute(dset_id, &NumSortedParticles, "NumberOfParticles", log_fptr, io_log);
 
	ReadAttribute(dset_id, &TotParticleCount, "TotalParticleCount", log_fptr, io_log);
 
	h5_status = H5Dclose(dset_id);
	if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	h5_status = H5Fclose(file_id);
	if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
	if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
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

	if (debug)
	  printf("Allocating %"ISYM" (NumberOfParticles)\n",NumberOfParticles);
 
	// Allocate buffer for floating point particle data
 
	inits_type *RingTemp = new inits_type[NumberOfParticles];
 
	// Allocate buffer for integer particle data
 
	int *IntRingTemp = new int[NumberOfParticles];
 
	// Read pre-sorted particle positions

	if (CosmologySimulationParticlePositionName != NULL) {
 
	  for (dim = 0; dim < GridRank; dim++) {
 
	    Slab_Rank = 2;
	    Slab_Dims[0] = 3;
	    Slab_Dims[1] = NumSortedParticles;
 
	    if ( NumSortedParticles == 0 )
	      Slab_Dims[1] = 1;
 
	    // Data in memory is considered 1D, stride 1, with zero offset
 
	    mem_stride = 1;                   // contiguous elements
	    mem_count = NumSortedParticles;   // number of elements in field
	    mem_offset = 0;                   // zero offset in buffer
 
	    if ( NumSortedParticles == 0 )
	      mem_count = 1;
 
	    // 1D memory model
 
	    mem_dsp_id = H5Screate_simple((Eint32) 1, &mem_count, NULL);
	    if (io_log) fprintf(log_fptr, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
	    if( mem_dsp_id == h5_error ){ENZO_FAIL("IO Problem");}
	      
	    h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
	    if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
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
	      slab_count[1] = 1;
 
	    file_dsp_id = H5Screate_simple((Eint32) Slab_Rank, Slab_Dims, NULL);
	    if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
	    if( file_dsp_id == h5_error ){ENZO_FAIL("IO Problem");}
	    
	    h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
	    if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	    file_id = H5Fopen(PPos, H5F_ACC_RDWR, H5P_DEFAULT);
	    fprintf(stderr, "H5Fopen %s on proc %"ISYM" status %"ISYM"\n", PPos, MyProcessorNumber, file_id);
	    if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
	    if( file_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	    dset_id =  H5Dopen(file_id, PPos);
	    if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	    if( dset_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	    if ( NumSortedParticles > 0 ) {
	      h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) RingTemp);
	      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
	    }
 
	    h5_status = H5Dclose(dset_id);
	    if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	    h5_status = H5Sclose(mem_dsp_id);
	    if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	    h5_status = H5Sclose(file_dsp_id);
	    if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	    h5_status = H5Fclose(file_id);
	    if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	    // Convert to FLOAT
 
	    for (j = 0; j < NumberOfParticles; j++)
	      ParticlePosition[dim][j] = RingTemp[j];
 
	  } // end of loop over dims
 
	} // end: if (CosmologySimulationParticlePositionName != NULL)

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
	      if( mem_dsp_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	      h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
	      if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
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
	      if( file_dsp_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	      h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
	      if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	      file_id = H5Fopen(PVel, H5F_ACC_RDWR, H5P_DEFAULT);
	      fprintf(stderr, "H5Fopen %s on proc %"ISYM" status %"ISYM"\n", PVel, MyProcessorNumber, file_id);
	      if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
	      if( file_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	      dset_id =  H5Dopen(file_id, PVel);
	      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	      if( file_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	      if ( NumSortedParticles > 0 )
		{
		  h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) RingTemp);
		  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
		  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
		}
 
	      h5_status = H5Dclose(dset_id);
	      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	      h5_status = H5Sclose(mem_dsp_id);
	      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	      h5_status = H5Sclose(file_dsp_id);
	      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	      h5_status = H5Fclose(file_id);
	      if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	      // Convert to float
 
	      for (j = 0; j < NumberOfParticles; j++)
		{
		  ParticleVelocity[dim][j] = RingTemp[j];
		}
 
	    } // end of loop over dims
 
	} // end: if (CosmologySimulationParticleVelocityName != NULL)

 
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
	  if( mem_dsp_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
	  if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
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
	  if( file_dsp_id == h5_error ){ENZO_FAIL("IO Problem");}

	  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
	  if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	  file_id = H5Fopen(PMass, H5F_ACC_RDWR, H5P_DEFAULT);
	  if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
	  if( file_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	  dset_id =  H5Dopen(file_id, PMass);
	  if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	  if( file_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	  if ( NumSortedParticles > 0 )
	    {
	      h5_status = H5Dread(dset_id, type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) RingTemp);
	      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
	    }

	  h5_status = H5Dclose(dset_id);
	  if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}

	  h5_status = H5Sclose(mem_dsp_id);
	  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}

	  h5_status = H5Sclose(file_dsp_id);
	  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}

	  h5_status = H5Fclose(file_id);
	  if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}

	  // Convert to float
 
	  for (j = 0; j < NumberOfParticles; j++)
	    {
	      ParticleMass[j] = RingTemp[j];
	    }
 
	} // end: if (CosmologySimulationParticleMassName != NULL)

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
	  if( mem_dsp_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
	  if (io_log) fprintf(log_fptr, "H5Sselect mem slab: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
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
	  if( file_dsp_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, slab_offset, slab_stride, slab_count, NULL);
	  if (io_log) fprintf(log_fptr, "H5Sselect file slab: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	  file_id = H5Fopen(PType, H5F_ACC_RDWR, H5P_DEFAULT);
	  if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
	  if( file_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	  dset_id =  H5Dopen(file_id, PType);
	  if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	  if( file_id == h5_error ){ENZO_FAIL("IO Problem");}
 
	  if ( NumSortedParticles > 0 )
	    {
	      h5_status = H5Dread(dset_id, int_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, (VOIDP) IntRingTemp);
	      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	      if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
	    }
 
	  h5_status = H5Dclose(dset_id);
	  if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	  h5_status = H5Sclose(mem_dsp_id);
	  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	  h5_status = H5Sclose(file_dsp_id);
	  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	  h5_status = H5Fclose(file_id);
	  if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){ENZO_FAIL("IO Problem");}
 
	  // Convert to float
 
	  for (j = 0; j < NumberOfParticles; j++)
	    {
	      ParticleType[j] = IntRingTemp[j];
	    }
 
	}

	printf("De-allocate (TotalParticleCount)\n");
 
	delete [] RingTemp;
	delete [] IntRingTemp;

      } // end: if (ParallelRootGridIO && PresortedParticles == 1)


      if (CosmologySimulationCalculatePositions) {
	if (CosmologyInitializeParticles(CosmologySimulationParticleVelocityName, 
					 CosmologySimulationParticleDisplacementName,
					 CosmologySimulationParticleMassName,
					 CosmologySimulationParticleTypeName,
					 CosmologySimulationParticleVelocityNames,
					 CosmologySimulationParticleDisplacementNames,
					 CosmologySimulationOmegaBaryonNow,
					 Offset, level) == FAIL) {
	  ENZO_FAIL("Error in grid::CosmologyInitializePositions.\n");
	}
      }


      if (PreSortedParticles == 0 && !CosmologySimulationCalculatePositions) {

	// For regular IO (parallel or non-parallel)
      
	printf("Standard Parallel or Serial I/O\n");
 
	// This code currently used in reading nested grids below top
 
	if (io_log) fprintf(log_fptr, "CSIG Parallel/Serial I/O\n");
 
	if (ParallelRootGridIO == TRUE && ReadData) {

	  /* If doing parallel root grid IO then divide up the list of particles.  These
	     will not be the correct particles for this processors in general, but we
	     will rearrange them during the rebuild.  (For ring I/O see above). */

	  Start[0] = nint(double(MyProcessorNumber)*double(TotalParticleCount)/double(NumberOfProcessors));
	  End[0] = nint(double(MyProcessorNumber+1)*double(TotalParticleCount)/double(NumberOfProcessors)) - 1;
	  printf("P(%d): ParticleCount = %d\n", MyProcessorNumber, TotalParticleCount);
	  printf("P(%d): reading particles %d to %d\n", MyProcessorNumber, Start[0], End[0]);

	  if (MyProcessorNumber+1 == NumberOfProcessors)
	    assert (End[0]+1 == TotalParticleCount);

	} else {

	  /* For regular IO (non-parallel) */

	  End[0] = TotalParticleCount - 1;
	  Dim[0] = TotalParticleCount;

	} // end: if (ParallelRootGridIO)

	// Set Number Of particles
 
	Dim[0] = End[0] - Start[0] + 1;
	NumberOfParticles = Dim[0];
 
	/* Set Offset and reset start and end so that start is 0. */

	Offset[0] = Start[0];
	Start[0] -= Offset[0];
	End[0]   -= Offset[0];

	// Allocate space for the particles
 
	this->AllocateNewParticles(NumberOfParticles);
 
	// Read particle positions
 
	for (dim = 0; dim < GridRank; dim++) {
	  if (READFILE(CosmologySimulationParticlePositionName, 1, Dim,
		       Start, End, Offset, NULL, &tempbuffer, dim, 3) == FAIL) {
	    ENZO_VFAIL("Error reading particle position %"ISYM".\n", dim)
	  }
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticlePosition[dim][i] = FLOAT(tempbuffer[i]);
 
	  if (io_log) fprintf(log_fptr, "De-Allocate [] tempbuffer for dim = %"ISYM"\n", dim);
 
	  delete [] tempbuffer;
	}
 
	// Read particle velocities
 
	if (CosmologySimulationParticleVelocityName != NULL)
	  for (dim = 0; dim < GridRank; dim++) {
	    if (READFILE(CosmologySimulationParticleVelocityName, 1, Dim,
			 Start, End, Offset, ParticleVelocity[dim], &tempbuffer, dim, 3) == FAIL) {
	      ENZO_VFAIL("Error reading particle velocity %"ISYM".\n", dim)
	    }
	    //        fcol(ParticleVelocity[dim], NumberOfParticles, 10, log_fptr);
	  }
 
	// Read particle mass
 
	if (CosmologySimulationParticleMassName != NULL) {
	  if (READFILE(CosmologySimulationParticleMassName, 1, Dim, Start, End,
		       Offset, ParticleMass, &tempbuffer, 0, 1) == FAIL) {
	    ENZO_FAIL("Error reading particle masses.\n");
	  }
	  //      fcol(ParticleMass, NumberOfParticles, 10, log_fptr);
	}
 
	// Read particle type
 
	if (CosmologySimulationParticleTypeName != NULL) {
	  if (ReadIntFile(CosmologySimulationParticleTypeName, 1, Dim, Start, End,
			  Offset, ParticleType, &int_tempbuffer, 0, 1) == FAIL) {
	    ENZO_FAIL("Error reading particle types.\n");
	  }
	  //      icol(ParticleType, NumberOfParticles, 10, log_fptr);
	}

      } // end: if (PreSortedParticles == 0 && !CosmologySimCalculatePositions)
 

      //    } // end: if ParallelRootGridIO == TRUE
 
      // If there are particles, but no velocity file then abort
      // or set the velocities to zero
 
      if (NumberOfParticles > 0 && 
	  (CosmologySimulationParticleVelocityName == NULL &&
	   CosmologySimulationParticleVelocityNames[0] == NULL)) {
	ENZO_FAIL("Error -- no velocity field specified.\n");
	//  printf("CosmologySimulation warning: setting velocities to zero.\n");
	//      for (dim = 0; dim < GridRank; dim++)
	//	for (i = 0; i < NumberOfParticles; i++)
	//	  ParticleVelocity[dim][i] = 0.0;
      }
 
      // If there are particles, but a mass file name wasn't specified,
      // make up the particles mass
 
      if (NumberOfParticles > 0 && CosmologySimulationParticleMassName == NULL) {
 
	// Compute mass of each particle
 
	float UniformParticleMass = CosmologySimulationOmegaCDMNow/
	  OmegaMatterNow;
 
	// if user wants to manually specify the mass ratio, great.  Otherwise,
	// go ahead and try to figure out what's going on.

	if (CosmologySimulationManuallySetParticleMassRatio == FALSE) {

	  // If there are exactly 1/8 as many particles as cells,
	  // then set the particle mass to 8 times the usual
    
	  int NumberOfActiveCells = (GridEndIndex[0]-GridStartIndex[0]+1)*
	    (GridEndIndex[1]-GridStartIndex[1]+1)*
	    (GridEndIndex[2]-GridStartIndex[2]+1);
	  if (NumberOfParticles*8 == NumberOfActiveCells)
	    UniformParticleMass *= 8;
	  if (NumberOfParticles == NumberOfActiveCells*8)
	    UniformParticleMass /= 8;
 
	  //      UniformParticleMass *= float(POW(TotalRefinement, GridRank));
 
	  // Issue a warning if PPIO or PRGIO are on (possibility of errors
	  // being caused)
	  if( ((ParallelParticleIO == TRUE) || (ParallelRootGridIO == TRUE)) &&
	      MyProcessorNumber == ROOT_PROCESSOR) {
	    fprintf(stderr,"\n\n\n*********************************************\n");
	    fprintf(stderr,"NestedCosmologySimulationInitializeGrid: WARNING!\n");
	    fprintf(stderr,"ParallelParticleIO or ParallelRootGridIO are ON and\n");
	    fprintf(stderr,"you are not manually setting particle mass.  You want\n");
	    fprintf(stderr,"to double-check that particle masses are being \n");
	    fprintf(stderr,"calculated correctly!\n");
	    fprintf(stderr,"Your particle mass is being calculated as:  %"GSYM"\n",
		    UniformParticleMass);
	    fprintf(stderr,"\n*********************************************\n\n\n");
	  }

	} else {

	  UniformParticleMass *= CosmologySimulationManualParticleMassRatio;

	  if (debug) {
	    printf("Particle mass ratio is: %"GSYM" so particle mass is %"GSYM"\n",
		   CosmologySimulationManualParticleMassRatio,
		   UniformParticleMass);
	    fflush(stdout);
	  }
	}


	// Set all particle masses to this value
 
	if (debug) 
	  printf("NestedCosmologySimulationInitializeGrid: particle mass = %"GSYM"\n",
		 UniformParticleMass);
 
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleMass[i] = UniformParticleMass;
 
      } // end: if (NumberOfParticles > 0 && CosmologySimulationParticleMassName == NULL)
 



      // Set Particle attributes to FLOAT_UNDEFINED
 
      for (j = 0; j < NumberOfParticleAttributes; j++)
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleAttribute[j][i] = FLOAT_UNDEFINED;


 
      // If no type file was read, then set all particles to dm particles
     

      if (CosmologySimulationParticleTypeName == NULL){
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleType[i] = PARTICLE_TYPE_DARK_MATTER;
      }

 
      // Assign particles unique identifiers
 
      if (ParallelRootGridIO == TRUE) {
	if (ParallelParticleIO || CosmologySimulationCalculatePositions) {
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticleNumber[i] = i; // Unique ID's will be added later
	  CurrentParticleNumber += NumberOfParticles; // ignored
	} else {
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticleNumber[i] = CurrentParticleNumber + i + Offset[0]; // Unique ID's calculated here!
          printf("P(%d): CurrentParticleNumber = %d  Offset = %d\n", MyProcessorNumber, CurrentParticleNumber, Offset[0]);
	  CurrentParticleNumber += TotalParticleCount; // set this so the calling routine knows the total number of particles on this level
	}
      } else {
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleNumber[i] = CurrentParticleNumber + i;
	CurrentParticleNumber += NumberOfParticles; // ignored
      }
 



  
    } // end: if (CosmologySimulationParticleName != NULL)

  } // end: if (ProcessorNumber == MyProcessorNumber)
 
  // Share number of particles amoung processors
 
  //  The following causes I/O serialization and should be bypassed
  //  for ParallelRootGridIO = 1

#ifdef UNUSED 
  if ( ParallelRootGridIO != TRUE )
    {
      if (ReadData == TRUE)
	CommunicationBroadcastValue(&NumberOfParticles, ProcessorNumber);
    }
#endif
 
  OldTime = Time;

  if (io_log) fclose(log_fptr);

 
  return SUCCESS;
 
}
