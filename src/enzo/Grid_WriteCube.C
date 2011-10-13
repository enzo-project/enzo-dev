/***********************************************************************
/
/  GRID CLASS (WRITE OUT UNIGRID DATA CUBES)
/
/  written by: Robert Harkness
/  date:       August, 2004
/  modified1:  Robert Harkness
/              May, 2008
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include <mpi.h>
#endif
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
#include "CommunicationUtilities.h"

void my_exit(int status);
 
// HDF5 function prototypes
 

 
// function prototypes
 
int FindCube(char *cube_name);
int FindField(int f, int farray[], int n);
int WriteStringAttr(hid_t dset_id, char *Alabel, char *String, FILE *log_fptr);

 
 
 
int grid::WriteCube(char *base_name, int grid_id, int TGdims[])
{

  int i, j, k, dim, field, size, ActiveDim[MAX_DIMENSION];
  int file_status;
  int output_cube;
 
  int nop[MAX_NUMBER_OF_TASKS], nopout[MAX_NUMBER_OF_TASKS];
  int pof[MAX_NUMBER_OF_TASKS];
  int TCount;
  int npe, ipe;
 
  float32 *temp;
  int *tempint;
  PINT *tempPINT;
 
  FILE *log_fptr;
 
  hid_t       file_id, dset_id;
  hid_t       float_type_id, FLOAT_type_id;
  hid_t       file_type_id, FILE_type_id;
  hid_t       mem_type_id;
  hid_t       mem_dsp_id, file_dsp_id;
  hid_t       attr_id, attr_type, attr_dsp_id;
 
  hid_t       file_access_template;
  hid_t       xfer_prop_list;
 
  hsize_t     OutDim[MAX_DIMENSION];
  hsize_t     InDim[MAX_DIMENSION];
  hsize_t     TempIntArray[1];
 
  hsize_t     dbuff_size;
  hsize_t     gbuff_size;
  hsize_t     m_size, l_size;
 
  hsize_t     mem_stride, mem_count;
  hsize_t     m_file_stride, m_file_count;
  hsize_t     file_stride[3], file_count[3], file_block[3];
  hsize_t     cube_stride[3], cube_count[3], cube_block[3];
 
  hssize_t    mem_offset;
  hssize_t    m_file_offset;
  hssize_t    file_offset[3];
  hssize_t    cube_offset[3];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int gridsize[3];
  int ndims;
  int StartIndex[3], EndIndex[3];
 
  float Left[3];
  float Right[3];
  float eps = 1.0e-06;
 
  char id[10];
  char FieldName[80+1];
  char GlueFile[80+1];
  char PartName[80+1];
  char LogName[80+1];
 
  char *ParticlePositionLabel[] =
     {"particle_position_x", "particle_position_y", "particle_position_z"};
 
  char *ParticleVelocityLabel[] =
     {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
 
  char *ParticleMassLabel = "particle_mass";
  char *ParticleTypeLabel = "particle_type";
  char *ParticleIndexLabel = "particle_index";
#ifdef WINDS
    char *ParticleAttributeLabel[] = 
      {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x", 
       "particle_jet_y", "particle_jet_z", "typeia_fraction"};
#else
    char *ParticleAttributeLabel[] = 
      {"creation_time", "dynamical_time", "metallicity_fraction", "typeia_fraction"};
#endif
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
 
  io_log = 1;
 
// Begin
 
#ifdef PARALLEL_HDF5
 
  int ii = sizeof(float32);
 
  switch(ii)
  {
 
    case 4:
      mem_type_id = HDF5_R4;
      float_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
      break;
 
    case 8:
      mem_type_id = HDF5_R8;
      float_type_id = HDF5_R8;
      file_type_id = HDF5_FILE_R8;
      break;
 
    default:
      mem_type_id = HDF5_R4;
      float_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
 
  }
 
  int jj = sizeof(FLOAT);
 
  switch(jj)
  {
 
    case 4:
      FLOAT_type_id = HDF5_R4;
      FILE_type_id = HDF5_FILE_R4;
      break;
 
    case 8:
      FLOAT_type_id = HDF5_R8;
      FILE_type_id = HDF5_FILE_R8;
      break;
 
    case 16:
      FLOAT_type_id = HDF5_R16;
      FILE_type_id = H5Tcopy(HDF5_FILE_B8);
                     H5Tset_size(FILE_type_id,16);
      break;
 
    default:
      fprintf(stderr, "INCORRECT FLOAT DEFINITION\n");
 
  }
 
  // Make sure quantities defined at least for 3d
 
  for (dim = GridRank; dim < 3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }
 
  for (dim = 0; dim < 3; dim++)
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;
 
  ndims = GridRank;
 
  if (MyProcessorNumber == ProcessorNumber)
  {
 
    sprintf(id, "%"GRID_TAG_FORMAT""ISYM, grid_id);
 
    strcpy(LogName, base_name);
    strcat(LogName, ".cubelog");
    strcat(LogName, id);
 
    if (io_log) {
      log_fptr = fopen(LogName, "a");
      fprintf(log_fptr, "Grid_WriteCube\n");
      fprintf(log_fptr, "  ID %"ISYM"  %s\n", grid_id, id);
      fprintf(log_fptr, "  BASE %s\n", base_name);
      fprintf(log_fptr, "  Left  %12.4"FSYM"  %12.4"FSYM"  %12.4"FSYM"\n", GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
      fprintf(log_fptr, "  Right %12.4"FSYM"  %12.4"FSYM"  %12.4"FSYM"\n", GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
    }
 
  }
 
 
  if (MyProcessorNumber == ProcessorNumber)
  {
#ifdef USE_MPI
    MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
#endif
    MPI_Arg mpi_size;
    MPI_Arg Count;
 
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else
    mpi_size = 1;
#endif
    npe = mpi_size;
    Count = npe;
 
    for( ipe=0; ipe < npe; ipe++)
    {
      nop[ipe] = 0;
      pof[ipe] = 0;
    }
 
    nop[MyProcessorNumber] = NumberOfParticles;
 
#ifdef USE_MPI
    MPI_Allreduce(nop, nopout, Count, DataTypeInt, MPI_SUM, MPI_COMM_WORLD);
    for( ipe=0; ipe < npe; ipe++)
    {
      nop[ipe] = nopout[ipe];
    }
#endif
 
    TCount = 0;
 
    for( ipe=0; ipe < npe; ipe++)
    {
      pof[ipe] = TCount;
      TCount = TCount + nop[ipe];
    }
 
//    fprintf(stderr, "nop %"ISYM": %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n", MyProcessorNumber,
//      nop[0], nop[1], nop[2], nop[3], nop[4], nop[5], nop[6], nop[7]);
//    fprintf(stderr, "pof %"ISYM": %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n", MyProcessorNumber,
//      pof[0], pof[1], pof[2], pof[3], pof[4], pof[5], pof[6], pof[7]);
//    fprintf(stderr, "total particle count: %"ISYM"\n", TCount);
 
    //  fprintf(stderr, "Call Barrier 1 on task %"ISYM", grid %"ISYM"\n", MyProcessorNumber, grid_id);
    CommunicationBarrier();
    //  fprintf(stderr, "Call Barrier 1 on task %"ISYM", grid %"ISYM" complete\n", MyProcessorNumber, grid_id);
  }
 
 
 
  if (MyProcessorNumber == ProcessorNumber)
  {
 
    // Every task makes it through the barrier once with its local grid boundary.
    // It should be safe to do a collective open now
 
    // Set output_dims[]
 
    gridsize[0] = TGdims[0];
    gridsize[1] = TGdims[1];
    gridsize[2] = TGdims[2];
 
    OutDim[0] = TGdims[0];  //TopGridDims[0];
    OutDim[1] = TGdims[1];  //TopGridDims[1];
    OutDim[2] = TGdims[2];  //TopGridDims[2];
 
 
    Left[0] = GridLeftEdge[0];
    Left[1] = GridLeftEdge[1];
    Left[2] = GridLeftEdge[2];
 
    Right[0] = GridRightEdge[0];
    Right[1] = GridRightEdge[1];
    Right[2] = GridRightEdge[2];
 
    StartIndex[0] = (int)(gridsize[0] * Left[0] + eps);
    EndIndex[0] = (int)(gridsize[0] * Right[0] + eps) - 1;
    StartIndex[1] = (int)(gridsize[1] * Left[1] + eps);
    EndIndex[1] = (int)(gridsize[1] * Right[1] + eps) - 1;
    StartIndex[2] = (int)(gridsize[2] * Left[2] + eps);
    EndIndex[2] = (int)(gridsize[2] * Right[2] + eps) - 1;
 
    InDim[0] = (GridEndIndex[0] - GridStartIndex[0]) + 1;
    InDim[1] = (GridEndIndex[1] - GridStartIndex[1]) + 1;
    InDim[2] = (GridEndIndex[2] - GridStartIndex[2]) + 1;
 
    dbuff_size = 1;
 
    for (i = 0; i < ndims; i++)
    {
      dbuff_size = dbuff_size * InDim[i];
    }
 
  }
 
 
  if (MyProcessorNumber == ProcessorNumber)
  {
 
    // For each Baryon field required
 
    for (field = 0; field < NumberOfBaryonFields; field++) {
 
    output_cube = FindCube(DataLabel[field]);
 
    if ( output_cube > -1 ) {
 
    // Set new GlueFile with field name
 
    // field = FindField(Density, FieldType, NumberOfBaryonFields);
 
    // DataLabel[field];
 
    // strcpy(FieldName, "Density");
    // strcpy(GlueFile, base_name);
    // strcat(GlueFile, ".Density");
 
    strcpy(FieldName, DataLabel[field]);
    strcpy(GlueFile, base_name);
    strcat(GlueFile, ".");
    strcat(GlueFile, DataLabel[field]);
 
    if (io_log) fprintf(log_fptr, "Field name = %s\n", FieldName);
    if (io_log) fprintf(log_fptr, "GlueFile name = %s\n", GlueFile);
 
    // Set file access template and switch to MPI I/O mode
 
    file_access_template = H5Pcreate (H5P_FILE_ACCESS);
      if( file_access_template == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, MPI_INFO_NULL);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    mem_dsp_id = H5Screate_simple((Eint32) 3, InDim, NULL);
        if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    file_dsp_id = H5Screate_simple((Eint32) 3, OutDim, NULL);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_template);
      if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Pclose(file_access_template);
      fprintf(log_fptr, "H5Pclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    dset_id = H5Dcreate(file_id, FieldName, file_type_id, file_dsp_id, H5P_DEFAULT);
      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    // Write out co-ordinate values.  Use the centre of each cell.
 
    size = dbuff_size;
    temp = new float32[size];
 
    // Copy active part of field into grid
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           +
	         (j-GridStartIndex[1])*ActiveDim[0]              +
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		       float32(
	      BaryonField[field][i + j*GridDimension[0] +
		                     k*GridDimension[0]*GridDimension[1]]
                              );
 
    // write dataset
 
    file_offset[0] = StartIndex[2];
    file_stride[0] = 1;
    file_count[0] = InDim[2];
    file_block[0] = 1;
 
    file_offset[1] = StartIndex[1];
    file_stride[1] = 1;
    file_count[1] = InDim[1];
    file_block[2] = 1;
 
    file_offset[2] = StartIndex[0];
    file_stride[2] = 1;
    file_count[2] = InDim[0];
    file_block[2] = 1;
 
    mem_offset = 0;
    mem_stride = 1;
    mem_count = dbuff_size;
 
    cube_offset[0] = 0;
    cube_stride[0] = 1;
    cube_count[0] = InDim[0];
 
    cube_offset[1] = 0;
    cube_stride[1] = 1;
    cube_count[1] = InDim[1];
 
    cube_offset[2] = 0;
    cube_stride[2] = 1;
    cube_count[2] = InDim[2];
 
    h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, cube_offset, cube_stride, cube_count, NULL);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
    h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
      if( xfer_prop_list == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, (VOIDP) temp);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Dclose(dset_id);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Sclose(mem_dsp_id);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Sclose(file_dsp_id);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Pclose(xfer_prop_list);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Fclose(file_id);
      fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    delete temp;
 
    } // if cube is in active output list
 
    } // end of loop over fields
 
  } // if MyProcessorNumber
 
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    //  fprintf(stderr, "Call Barrier 2 on task %"ISYM", grid %"ISYM"\n", MyProcessorNumber, grid_id);
    CommunicationBarrier();
    //  fprintf(stderr, "Call Barrier 2 on task %"ISYM", grid %"ISYM" complete\n", MyProcessorNumber, grid_id);
  }
 
 
  // repeat for temperature
 
  if (MyProcessorNumber == ProcessorNumber)
  {
 
    output_cube = FindCube("Temperature");
 
    if ( output_cube > -1 ) {
 
    if (ComovingCoordinates) {
 
      strcpy(FieldName, "Temperature");
      strcpy(GlueFile, base_name);
      strcat(GlueFile, ".Temperature");
 
      if (io_log) fprintf(log_fptr, "Field name = %s\n", FieldName);
      if (io_log) fprintf(log_fptr, "GlueFile name = %s\n", GlueFile);
 
      // Set file access template and switch to MPI I/O mode
 
      file_access_template = H5Pcreate (H5P_FILE_ACCESS);
        if( file_access_template == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, MPI_INFO_NULL);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      mem_dsp_id = H5Screate_simple((Eint32) 3, InDim, NULL);
        if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      file_dsp_id = H5Screate_simple((Eint32) 3, OutDim, NULL);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_template);
        if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Pclose(file_access_template);
        fprintf(log_fptr, "H5Pclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      dset_id = H5Dcreate(file_id, FieldName, file_type_id, file_dsp_id, H5P_DEFAULT);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      // Write out co-ordinate values.  Use the centre of each cell.
 
      size = dbuff_size;
      temp = new float32[size];
 
      int tsize;
      tsize = 1;
 
      for (dim = 0; dim < GridRank; dim++)
        tsize *= GridDimension[dim];
 
      float *temperature = new float[tsize];
 
      if (this->ComputeTemperatureField(temperature) == FAIL) {
	ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");
      }
 
      // Copy active part of field into grid
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           +
	         (j-GridStartIndex[1])*ActiveDim[0]              +
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     float32(
		   temperature[(k*GridDimension[1] + j)*GridDimension[0] + i]
			     );
 
      file_offset[0] = StartIndex[2];
      file_stride[0] = 1;
      file_count[0] = InDim[2];
      file_block[0] = 1;
 
      file_offset[1] = StartIndex[1];
      file_stride[1] = 1;
      file_count[1] = InDim[1];
      file_block[2] = 1;
 
      file_offset[2] = StartIndex[0];
      file_stride[2] = 1;
      file_count[2] = InDim[0];
      file_block[2] = 1;
 
      mem_offset = 0;
      mem_stride = 1;
      mem_count = dbuff_size;
 
      cube_offset[0] = 0;
      cube_stride[0] = 1;
      cube_count[0] = InDim[0];
 
      cube_offset[1] = 0;
      cube_stride[1] = 1;
      cube_count[1] = InDim[1];
 
      cube_offset[2] = 0;
      cube_stride[2] = 1;
      cube_count[2] = InDim[2];
 
      h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, cube_offset, cube_stride, cube_count, NULL);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
      h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
        if( xfer_prop_list == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, (VOIDP) temp);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(mem_dsp_id);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Pclose(xfer_prop_list);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Fclose(file_id);
        fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      delete temperature;
      delete temp;
 
    } // end: if (ComovingCoordinates)
 
    } // if output cube active
 
  } // end: if MyProc
 
 
 
 
  // repeat for dark matter
 
  if (MyProcessorNumber == ProcessorNumber)
  {
 
    output_cube = FindCube("Dark_Matter_Density");
 
    if ( output_cube > -1 ) {
 
    if (SelfGravity && NumberOfParticles > 0) {
      float SaveGravityResolution = GravityResolution;
      GravityResolution = 1;
      this->InitializeGravitatingMassFieldParticles(RefineBy);
      this->ClearGravitatingMassFieldParticles();
      this->DepositParticlePositions(this, Time, GRAVITATING_MASS_FIELD_PARTICLES);
      GravityResolution = SaveGravityResolution;
    }
 
    // If present, write out the GravitatingMassFieldParticles
 
    if (GravitatingMassFieldParticles != NULL) {
 
      int GravStartIndex[] = {0,0,0}, GravEndIndex[] = {0,0,0};
 
      for (dim = 0; dim < GridRank; dim++) {
 
	GravStartIndex[dim] = nint((GridLeftEdge[dim] -
			            GravitatingMassFieldParticlesLeftEdge[dim])/
			           GravitatingMassFieldParticlesCellSize);
 
	GravEndIndex[dim] = nint((GridRightEdge[dim] -
			          GravitatingMassFieldParticlesLeftEdge[dim])/
			         GravitatingMassFieldParticlesCellSize) - 1;
 
        //  fprintf(stderr, "%"ISYM" %"ISYM" %10.4"FSYM" %10.4"FSYM" %10.4"FSYM" %10.4"FSYM"\n",
        //          GravStartIndex[dim], GravEndIndex[dim],
        //          GridLeftEdge[dim], GridRightEdge[dim],
        //          GravitatingMassFieldParticlesLeftEdge[dim],
        //          GravitatingMassFieldParticlesCellSize);
 
      }
 
      strcpy(FieldName, "Dark_Matter_Density");
      strcpy(GlueFile, base_name);
      strcat(GlueFile, ".Dark_Matter_Density");
 
      if (io_log) fprintf(log_fptr, "Field name = %s\n", FieldName);
      if (io_log) fprintf(log_fptr, "GlueFile name = %s\n", GlueFile);
 
      // Set file access template and switch to MPI I/O mode
 
      file_access_template = H5Pcreate (H5P_FILE_ACCESS);
        if( file_access_template == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, MPI_INFO_NULL);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      mem_dsp_id = H5Screate_simple((Eint32) 3, InDim, NULL);
        if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      file_dsp_id = H5Screate_simple((Eint32) 3, OutDim, NULL);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_template);
        if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Pclose(file_access_template);
        fprintf(log_fptr, "H5Pclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      dset_id = H5Dcreate(file_id, FieldName, file_type_id, file_dsp_id, H5P_DEFAULT);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      // Write out co-ordinate values.  Use the centre of each cell.
 
      size = dbuff_size;
 
      temp = new float32[size];
 
      // Copy active part of field into grid
 
      for (k = GravStartIndex[2]; k <= GravEndIndex[2]; k++)
	for (j = GravStartIndex[1]; j <= GravEndIndex[1]; j++)
	  for (i = GravStartIndex[0]; i <= GravEndIndex[0]; i++)
	    temp[(i-GravStartIndex[0])                           +
	         (j-GravStartIndex[1])*ActiveDim[0]              +
	         (k-GravStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     float32(
			     GravitatingMassFieldParticles[ i +
			       j*GravitatingMassFieldParticlesDimension[0] +
			       k*GravitatingMassFieldParticlesDimension[0]*
			         GravitatingMassFieldParticlesDimension[1]]
			     );
 
 
      // write dataset
 
      file_offset[0] = StartIndex[2];
      file_stride[0] = 1;
      file_count[0] = InDim[2];
      file_block[0] = 1;
 
      file_offset[1] = StartIndex[1];
      file_stride[1] = 1;
      file_count[1] = InDim[1];
      file_block[2] = 1;
 
      file_offset[2] = StartIndex[0];
      file_stride[2] = 1;
      file_count[2] = InDim[0];
      file_block[2] = 1;
 
      mem_offset = 0;
      mem_stride = 1;
      mem_count = dbuff_size;
 
      cube_offset[0] = 0;
      cube_stride[0] = 1;
      cube_count[0] = InDim[0];
 
      cube_offset[1] = 0;
      cube_stride[1] = 1;
      cube_count[1] = InDim[1];
 
      cube_offset[2] = 0;
      cube_stride[2] = 1;
      cube_count[2] = InDim[2];
 
      h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, cube_offset, cube_stride, cube_count, NULL);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
      h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
        if( xfer_prop_list == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dwrite(dset_id, mem_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, (VOIDP) temp);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(mem_dsp_id);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Pclose(xfer_prop_list);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Fclose(file_id);
        fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      // Clean up if we modified the resolution
 
      if (SelfGravity && GravityResolution != 1)
	this->DeleteGravitatingMassFieldParticles();
 
    } // end of (if GravitatingMassFieldParticles != NULL)
 
    delete temp;
 
    } // if output cube active
 
  } // end: if MyProcessor
 
 
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    //  fprintf(stderr, "Call Barrier 3 on task %"ISYM", grid %"ISYM"\n", MyProcessorNumber, grid_id);
    CommunicationBarrier();
    //  fprintf(stderr, "Call Barrier 3 on task %"ISYM", grid %"ISYM" complete\n", MyProcessorNumber, grid_id);
  }
 
 
 
  // particles
 
  if (NumberOfParticles > 0) {
 
    // get the total particle count TCount
 
    if (MyProcessorNumber == ProcessorNumber) {
 
      // Sort particles according to their identifier
 
      this->SortParticlesByNumber();
 
      // Particle positions are not converted to 32 bit first.
      // (128 bit numbers are not supported by HDF so convert to 64)
 
      float64 *temp_pointer = NULL;
 
      if (sizeof(FLOAT) == 16)
        temp_pointer = new float64[NumberOfParticles];
 
      TempIntArray[0] = NumberOfParticles;
 
      for (dim = 0; dim < GridRank; dim++) {
 
        output_cube = FindCube(ParticlePositionLabel[dim]);
 
        if ( output_cube > -1 ) {
 
        // Convert to 64 if 128, otherwise just write out
 
        if (sizeof(FLOAT) == 16)
	  for (i = 0; i < NumberOfParticles; i++)
	    temp_pointer[i] = float64(ParticlePosition[dim][i]);
        else
	  temp_pointer = (float64*) ParticlePosition[dim];
 
        m_size = TCount;
 
        strcpy(PartName, ParticlePositionLabel[dim]);
        strcpy(GlueFile, base_name);
        strcat(GlueFile, ".");
        strcat(GlueFile, ParticlePositionLabel[dim]);
 
        if (io_log) fprintf(log_fptr, "Field name = %s\n", PartName);
        if (io_log) fprintf(log_fptr, "GlueFile name = %s\n", GlueFile);
 
        file_access_template = H5Pcreate (H5P_FILE_ACCESS);
          if( file_access_template == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, MPI_INFO_NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_dsp_id = H5Screate_simple((Eint32) 1, &m_size, NULL);
          if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_template);
          if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        dset_id = H5Dcreate(file_id, PartName, FILE_type_id, file_dsp_id, H5P_DEFAULT);
          if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
          if( xfer_prop_list == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        l_size = nop[MyProcessorNumber];
 
        mem_dsp_id = H5Screate_simple((Eint32) 1, &l_size, NULL);
          if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        m_file_offset = pof[MyProcessorNumber];
        m_file_stride = 1;
        m_file_count = nop[MyProcessorNumber];
 
        mem_offset = 0;
        mem_stride = 1;
        mem_count = nop[MyProcessorNumber];
 
        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &m_file_offset, &m_file_stride, &m_file_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dwrite(dset_id, FLOAT_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, (VOIDP) temp_pointer);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(mem_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(file_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dclose(dset_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(file_access_template);
          fprintf(log_fptr, "H5Pclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(xfer_prop_list);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
          fprintf(log_fptr, "H5Fflush: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fclose(file_id);
          fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        } // if output cube active for this dim
 
      } // end of loop over dim
 
 
    if (sizeof(FLOAT) == 16)
      delete [] temp_pointer;
 
 
    // particle velocity
 
    // Create a temporary buffer (32 bit)
 
    temp = new float32[NumberOfParticles];
 
    // Copy particle velocities to temp and write them
 
    for (dim = 0; dim < GridRank; dim++) {
 
      output_cube = FindCube(ParticleVelocityLabel[dim]);
 
      if ( output_cube > -1 ) {
 
      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleVelocity[dim][i]);
 
        strcpy(PartName, ParticleVelocityLabel[dim]);
        strcpy(GlueFile, base_name);
        strcat(GlueFile, ".");
        strcat(GlueFile, ParticleVelocityLabel[dim]);
 
        if (io_log) fprintf(log_fptr, "Field name = %s\n", PartName);
        if (io_log) fprintf(log_fptr, "GlueFile name = %s\n", GlueFile);
 
        file_access_template = H5Pcreate (H5P_FILE_ACCESS);
          if( file_access_template == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, MPI_INFO_NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_dsp_id = H5Screate_simple((Eint32) 1, &m_size, NULL);
          if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_template);
          if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        dset_id = H5Dcreate(file_id, PartName, file_type_id, file_dsp_id, H5P_DEFAULT);
          if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
          if( xfer_prop_list == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        l_size = nop[MyProcessorNumber];
 
        mem_dsp_id = H5Screate_simple((Eint32) 1, &l_size, NULL);
          if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        m_file_offset = pof[MyProcessorNumber];
        m_file_stride = 1;
        m_file_count = nop[MyProcessorNumber];
 
        mem_offset = 0;
        mem_stride = 1;
        mem_count = nop[MyProcessorNumber];
 
        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &m_file_offset, &m_file_stride, &m_file_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dwrite(dset_id, float_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, (VOIDP) temp);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(mem_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(file_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dclose(dset_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(file_access_template);
          fprintf(log_fptr, "H5Pclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(xfer_prop_list);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
          fprintf(log_fptr, "H5Fflush: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fclose(file_id);
          fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      } // if output cube active for this dim
 
    }
 
    delete temp;
 
 
    // particle mass
 
    output_cube = FindCube("particle_mass");
 
    if ( output_cube > -1 ) {
 
    temp = new float32[NumberOfParticles];
 
    for (i = 0; i < NumberOfParticles; i++)
      temp[i] = float32(ParticleMass[i]);
 
    strcpy(PartName, "particle_mass");
    strcpy(GlueFile, base_name);
    strcat(GlueFile, ".particle_mass");
 
    if (io_log) fprintf(log_fptr, "Field name = %s\n", PartName);
    if (io_log) fprintf(log_fptr, "GlueFile name = %s\n", GlueFile);
 
        file_access_template = H5Pcreate (H5P_FILE_ACCESS);
          if( file_access_template == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, MPI_INFO_NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_dsp_id = H5Screate_simple((Eint32) 1, &m_size, NULL);
          if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_template);
          if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        dset_id = H5Dcreate(file_id, PartName, file_type_id, file_dsp_id, H5P_DEFAULT);
          if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
          if( xfer_prop_list == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        l_size = nop[MyProcessorNumber];
 
        mem_dsp_id = H5Screate_simple((Eint32) 1, &l_size, NULL);
          if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        m_file_offset = pof[MyProcessorNumber];
        m_file_stride = 1;
        m_file_count = nop[MyProcessorNumber];
 
        mem_offset = 0;
        mem_stride = 1;
        mem_count = nop[MyProcessorNumber];
 
        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &m_file_offset, &m_file_stride, &m_file_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dwrite(dset_id, float_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, (VOIDP) temp);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(mem_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(file_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dclose(dset_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(file_access_template);
          fprintf(log_fptr, "H5Pclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(xfer_prop_list);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
          fprintf(log_fptr, "H5Fflush: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fclose(file_id);
          fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    delete temp;
 
    } // if output cube active
 
    // particle_index and particle_type
 
    output_cube = FindCube("particle_index");
 
    if ( output_cube > -1 ) {
 
    tempPINT = new PINT[NumberOfParticles];
 
    for (i = 0; i < NumberOfParticles; i++)
      tempint[i] = ParticleNumber[i];
 
    strcpy(PartName, "particle_index");
    strcpy(GlueFile, base_name);
    strcat(GlueFile, ".particle_index");
 
    if (io_log) fprintf(log_fptr, "Field name = %s\n", PartName);
    if (io_log) fprintf(log_fptr, "GlueFile name = %s\n", GlueFile);
 
        file_access_template = H5Pcreate (H5P_FILE_ACCESS);
          if( file_access_template == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, MPI_INFO_NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_dsp_id = H5Screate_simple((Eint32) 1, &m_size, NULL);
          if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_template);
          if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        dset_id = H5Dcreate(file_id, PartName, HDF5_FILE_PINT, file_dsp_id, H5P_DEFAULT);
          if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
          if( xfer_prop_list == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        l_size = nop[MyProcessorNumber];
 
        mem_dsp_id = H5Screate_simple((Eint32) 1, &l_size, NULL);
          if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        m_file_offset = pof[MyProcessorNumber];
        m_file_stride = 1;
        m_file_count = nop[MyProcessorNumber];
 
        mem_offset = 0;
        mem_stride = 1;
        mem_count = nop[MyProcessorNumber];
 
        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &m_file_offset, &m_file_stride, &m_file_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dwrite(dset_id, HDF5_PINT, mem_dsp_id, file_dsp_id, xfer_prop_list, (VOIDP) tempPINT);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(mem_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(file_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dclose(dset_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(file_access_template);
          fprintf(log_fptr, "H5Pclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(xfer_prop_list);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
          fprintf(log_fptr, "H5Fflush: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fclose(file_id);
          fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    delete tempint;
 
    } // if output cube active
 
    output_cube = FindCube("particle_type");
 
    if ( output_cube > -1 ) {
 
    tempint = new int[NumberOfParticles];
 
    for (i = 0; i < NumberOfParticles; i++)
      tempint[i] = ParticleType[i];
 
    strcpy(PartName, "particle_type");
    strcpy(GlueFile, base_name);
    strcat(GlueFile, ".particle_type");
 
    if (io_log) fprintf(log_fptr, "Field name = %s\n", PartName);
    if (io_log) fprintf(log_fptr, "GlueFile name = %s\n", GlueFile);
 
        file_access_template = H5Pcreate (H5P_FILE_ACCESS);
          if( file_access_template == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, MPI_INFO_NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_dsp_id = H5Screate_simple((Eint32) 1, &m_size, NULL);
          if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_template);
          if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        dset_id = H5Dcreate(file_id, PartName, HDF5_FILE_INT, file_dsp_id, H5P_DEFAULT);
          if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
          if( xfer_prop_list == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        l_size = nop[MyProcessorNumber];
 
        mem_dsp_id = H5Screate_simple((Eint32) 1, &l_size, NULL);
          if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        m_file_offset = pof[MyProcessorNumber];
        m_file_stride = 1;
        m_file_count = nop[MyProcessorNumber];
 
        mem_offset = 0;
        mem_stride = 1;
        mem_count = nop[MyProcessorNumber];
 
        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &m_file_offset, &m_file_stride, &m_file_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dwrite(dset_id, HDF5_INT, mem_dsp_id, file_dsp_id, xfer_prop_list, (VOIDP) tempint);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(mem_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(file_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dclose(dset_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(file_access_template);
          fprintf(log_fptr, "H5Pclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(xfer_prop_list);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
          fprintf(log_fptr, "H5Fflush: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fclose(file_id);
          fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    delete [] tempint;
    delete [] tempPINT;
 
    } // if output cube active
 
    // particle attributes
 
    temp = new float32[NumberOfParticles];
 
    for (j = 0; j < NumberOfParticleAttributes; j++) {
 
      output_cube = FindCube(ParticleAttributeLabel[j]);
 
      if ( output_cube > -1 ) {
 
      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleAttribute[j][i]);
 
      strcpy(PartName, ParticleAttributeLabel[j]);
      strcpy(GlueFile, base_name);
      strcat(GlueFile, ".");
      strcat(GlueFile, ParticleAttributeLabel[j]);
 
      if (io_log) fprintf(log_fptr, "Field name = %s\n", PartName);
      if (io_log) fprintf(log_fptr, "GlueFile name = %s\n", GlueFile);
 
        file_access_template = H5Pcreate (H5P_FILE_ACCESS);
          if( file_access_template == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_fapl_mpio(file_access_template, MPI_COMM_WORLD, MPI_INFO_NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_dsp_id = H5Screate_simple((Eint32) 1, &m_size, NULL);
          if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        file_id = H5Fcreate(GlueFile, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_template);
          if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        dset_id = H5Dcreate(file_id, PartName, file_type_id, file_dsp_id, H5P_DEFAULT);
          if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        xfer_prop_list = H5Pcreate (H5P_DATASET_XFER);
          if( xfer_prop_list == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pset_dxpl_mpio(xfer_prop_list, H5FD_MPIO_COLLECTIVE);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        l_size = nop[MyProcessorNumber];
 
        mem_dsp_id = H5Screate_simple((Eint32) 1, &l_size, NULL);
          if( mem_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
        m_file_offset = pof[MyProcessorNumber];
        m_file_stride = 1;
        m_file_count = nop[MyProcessorNumber];
 
        mem_offset = 0;
        mem_stride = 1;
        mem_count = nop[MyProcessorNumber];
 
        h5_status = H5Sselect_hyperslab(mem_dsp_id, H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, &m_file_offset, &m_file_stride, &m_file_count, NULL);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dwrite(dset_id, float_type_id, mem_dsp_id, file_dsp_id, xfer_prop_list, (VOIDP) temp);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(mem_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Sclose(file_dsp_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Dclose(dset_id);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(file_access_template);
          fprintf(log_fptr, "H5Pclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Pclose(xfer_prop_list);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fflush(file_id, H5F_SCOPE_GLOBAL);
          fprintf(log_fptr, "H5Fflush: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
        h5_status = H5Fclose(file_id);
          fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
          if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      } // if output cube active for this attribute
 
    }
 
    delete temp;
 
 
    } // end: if (MyProcessorNumber)
 
  } // end: if (NumberOfParticles > 0)
 
 
 
 
 
 
 
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    if (io_log) fclose(log_fptr);

    //  fprintf(stderr, "Call Barrier 4 on task %"ISYM", grid %"ISYM"\n", MyProcessorNumber, grid_id);
    CommunicationBarrier();
    //  fprintf(stderr, "Call Barrier 4 on task %"ISYM", grid %"ISYM" complete\n", MyProcessorNumber, grid_id);
  }
 
  //  fprintf(stderr, "Exit Grid_WriteCube on cpu %"ISYM"\n", MyProcessorNumber);
 
#endif
 
  return SUCCESS;
 
}
