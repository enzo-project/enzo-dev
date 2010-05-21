/***********************************************************************
/
/  GRID CLASS (WRITE OUT GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  Robert Harkness, April 2008
/  modified3:  Robert Harkness, May 2008
/
/  PURPOSE:
/
************************************************************************/
 
//  Write grid to file pointer fptr
//     (we assume that the grid is at an appropriate stopping point,
//      where the Old values aren't required)
 
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
 
// HDF5 function prototypes
 

 
// function prototypes
 
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int WriteStringAttr(hid_t dset_id, char *Alabel, char *String, FILE *log_fptr);
 
 
 
 
int grid::WriteGridX(FILE *fptr, char *base_name, int grid_id)
{
 
  int i, j, k, dim, field, size, ActiveDim[MAX_DIMENSION];
  int file_status;
 
  float32 *temp;
 
  FILE *log_fptr;
 
  hid_t       file_id, dset_id;
  hid_t       float_type_id, FLOAT_type_id;
  hid_t       file_type_id, FILE_type_id;
  hid_t       file_dsp_id;
 
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     TempIntArray[1];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  char *ParticlePositionLabel[] =
     {"particle_position_x", "particle_position_y", "particle_position_z"};
  char *ParticleVelocityLabel[] =
     {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
  char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
				    "metallicity_fraction", "particle_jet_x", "particle_jet_y", "particle_jet_z", "alpha_fraction"};
  /*  char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
      "metallicity_fraction", "alpha_fraction"};*/
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
  printf("THIS IS WGX \n");
 
  /* initialize */
 
  char id[MAX_GRID_TAG_SIZE];
  sprintf(id, "%"GRID_TAG_FORMAT""ISYM, grid_id);
 
  /* make sure quantities defined at least for 3d */
 
  for (dim = GridRank; dim < 3; dim++) {
    GridDimension[dim] = 1;
    GridStartIndex[dim] = 0;
    GridEndIndex[dim] = 0;
  }
 
  for (dim = 0; dim < 3; dim++)
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;
 
  /* ------------------------------------------------------------------- */
  /* 1) Save general grid class data */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
 
    fprintf(fptr, "GridRank          = %"ISYM"\n", GridRank);
 
    fprintf(fptr, "GridDimension     = ");
    WriteListOfInts(fptr, GridRank, GridDimension);
 
    fprintf(fptr, "GridStartIndex    = ");
    WriteListOfInts(fptr, GridRank, GridStartIndex);
 
    fprintf(fptr, "GridEndIndex      = ");
    WriteListOfInts(fptr, GridRank, GridEndIndex);
 
    fprintf(fptr, "GridLeftEdge      = ");
    WriteListOfFloats(fptr, GridRank, GridLeftEdge);
 
    fprintf(fptr, "GridRightEdge     = ");
    WriteListOfFloats(fptr, GridRank, GridRightEdge);
 
    fprintf(fptr, "Time              = %"GOUTSYM"\n", Time);
 
    fprintf(fptr, "SubgridsAreStatic = %"ISYM"\n", SubgridsAreStatic);
 
    fprintf(fptr, "NumberOfBaryonFields = %"ISYM"\n", NumberOfBaryonFields);
 
  }
 
  char *name = new char[MAX_NAME_LENGTH];
  strcpy(name, base_name);
  strcat(name, id);
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    char *logname = new char[MAX_NAME_LENGTH];
    strcpy(logname, base_name);
    strcat(logname, id);
    strcat(logname, ".log");
    if (io_log) log_fptr = fopen(logname, "a");
    delete logname;
 
    if (io_log) fprintf(log_fptr, "Grid_WriteGrid\n");
    if (io_log) fprintf(log_fptr, "  ID %"ISYM"  %s\n", grid_id, id);
    if (io_log) fprintf(log_fptr, "  BASE %s\n", base_name);
    if (io_log) fprintf(log_fptr, "  HDF %s\n", name);
  }
 
  /* Open HDF file for writing. */
 
 
//  if (MyProcessorNumber == ProcessorNumber)
//  {
//    if (io_log) fprintf(log_fptr, "H5Fcreate with Name = %s\n", name);
//
//    file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//
//    fprintf(stderr,"GPFS create: file_id = %"ISYM"  %s\n", file_id, name);
//
//    if (io_log) fprintf(log_fptr,"H5Fcreate id: %"ISYM"\n",file_id);
//    if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
//  }
 
 
  /* ------------------------------------------------------------------- */
  /* 2) save baryon field quantities (including fields). */
 
  int ii = sizeof(float32);
 
  switch(ii)
  {
 
    case 4:
      float_type_id = HDF5_R4;
      file_type_id = HDF5_FILE_R4;
      break;
 
    case 8:
      float_type_id = HDF5_R8;
      file_type_id = HDF5_FILE_R8;
      break;
 
    default:
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
      fprintf(stderr, "INCORRECT FLOAT DEFINITION in Grid_WriteGridX\n");
 
  }
 
  if (NumberOfBaryonFields > 0) {
 
    if (MyProcessorNumber == ROOT_PROCESSOR) {
 
      fprintf(fptr, "FieldType = ");
      WriteListOfInts(fptr, NumberOfBaryonFields, FieldType);
 
      fprintf(fptr, "BaryonFileName = %s\n", name);
 
      fprintf(fptr, "CourantSafetyNumber    = %"FSYM"\n", CourantSafetyNumber);
      fprintf(fptr, "PPMFlatteningParameter = %"ISYM"\n", PPMFlatteningParameter);
      fprintf(fptr, "PPMDiffusionParameter  = %"ISYM"\n", PPMDiffusionParameter);
      fprintf(fptr, "PPMSteepeningParameter = %"ISYM"\n", PPMSteepeningParameter);
 
    }
 
    if (MyProcessorNumber == ProcessorNumber) {
 
    /* 2a) Set HDF file dimensions (use FORTRAN ordering). */
 
    for (dim = 0; dim < GridRank; dim++)
      OutDims[GridRank-dim-1] = ActiveDim[dim];
 
    /* 2b) Write out co-ordinate values.  Use the centre of each cell. */
 
    size = 1;
    float32 *tempdim[MAX_DIMENSION];
 
    for (dim = GridRank-1; dim >= 0; dim--) {
 
      /* Compute cell centers and put them in temp. */
 
      tempdim[dim] = new float32[GridDimension[dim]];

      for (i = 0; i <= GridEndIndex[dim] - GridStartIndex[dim]; i++)
	tempdim[dim][i] = CellLeftEdge[dim][GridStartIndex[dim] + i] +
	          0.5 * CellWidth[dim][GridStartIndex[dim] + i];

      size *= GridDimension[dim];
    }
 
    /* create temporary buffer */
 
    temp = new float32[size];
 
    /* 2c) Loop over fields, writing each one. */
 
    for (field = 0; field < NumberOfBaryonFields; field++) {
 
      /* copy active part of field into grid */
 
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
 
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Fcreate with Name = %s\n", DataLabel[field]);
 
      file_id = H5Fcreate(DataLabel[field], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
        if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", DataLabel[field]);
 
      dset_id =  H5Dcreate(file_id, DataLabel[field], file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      /* set datafield name and units, etc. */
 
      if ( DataUnits[field] == NULL )
      {
        DataUnits[field] = "none";
      }
 
      WriteStringAttr(dset_id, "Label", DataLabel[field], log_fptr);
      WriteStringAttr(dset_id, "Units", DataUnits[field], log_fptr);
      WriteStringAttr(dset_id, "Format", "e10.4", log_fptr);
      WriteStringAttr(dset_id, "Geometry", "Cartesian", log_fptr);
 
      h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Fclose(file_id);
        if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    }   // end of loop over fields
 
    /* If this is cosmology, compute the temperature field as well since
       its such a pain to compute after the fact. */
 
    if (OutputTemperature) {
 
      /* Allocate field and compute temperature. */
 
      float *temperature = new float[size];
      if (this->ComputeTemperatureField(temperature) == FAIL) {
	ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");
      }
 
      /* Copy active part of field into grid */
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           +
	         (j-GridStartIndex[1])*ActiveDim[0]              +
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     float32(
		   temperature[(k*GridDimension[1] + j)*GridDimension[0] + i]
			     );
 
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Fcreate with Name = Temperature\n");
 
      file_id = H5Fcreate("Temperature", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
        if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Dcreate with Name = Temperature\n");
 
      dset_id = H5Dcreate(file_id, "Temperature", file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if ( DataUnits[field] == NULL )
      {
        DataUnits[field] = "none";
      }
 
      WriteStringAttr(dset_id, "Label", "Temperature", log_fptr);
      WriteStringAttr(dset_id, "Units", "K", log_fptr);
      WriteStringAttr(dset_id, "Format", "e10.4", log_fptr);
      WriteStringAttr(dset_id, "Geometry", "Cartesian", log_fptr);
 
      h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Fclose(file_id);
        if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      delete temperature;
 
    } // end: if (ComovingCoordinates)
 
    /* Make sure that there is a copy of dark matter field to save
       (and at the right resolution). */
 
    if (SelfGravity && NumberOfParticles > 0) {
      float SaveGravityResolution = GravityResolution;
      GravityResolution = 1;
      this->InitializeGravitatingMassFieldParticles(RefineBy);
      this->ClearGravitatingMassFieldParticles();
      this->DepositParticlePositions(this, Time,
				     GRAVITATING_MASS_FIELD_PARTICLES);
      GravityResolution = SaveGravityResolution;
    }
 
    /* If present, write out the GravitatingMassFieldParticles. */
 
    if (GravitatingMassFieldParticles != NULL) {
 
      /* Set dimensions. */
 
      int StartIndex[] = {0,0,0}, EndIndex[] = {0,0,0};
      for (dim = 0; dim < GridRank; dim++) {
	StartIndex[dim] = nint((GridLeftEdge[dim] -
				GravitatingMassFieldParticlesLeftEdge[dim])/
			       GravitatingMassFieldParticlesCellSize);
	EndIndex[dim] = nint((GridRightEdge[dim] -
			      GravitatingMassFieldParticlesLeftEdge[dim])/
			     GravitatingMassFieldParticlesCellSize) - 1;
      }
 
      /* Copy active part of field into grid */
 
      for (k = StartIndex[2]; k <= EndIndex[2]; k++)
	for (j = StartIndex[1]; j <= EndIndex[1]; j++)
	  for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	    temp[(i-StartIndex[0])                           +
	         (j-StartIndex[1])*ActiveDim[0]              +
	         (k-StartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     float32(
			     GravitatingMassFieldParticles[ i +
			       j*GravitatingMassFieldParticlesDimension[0] +
			       k*GravitatingMassFieldParticlesDimension[0]*
			         GravitatingMassFieldParticlesDimension[1]]
			     );
 
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Fcreate with Name = Dark_Matter_Density\n");
 
      file_id = H5Fcreate("Dark_Matter_Density", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
        if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Dcreate with Name = Dark_Matter_Density\n");
 
      dset_id = H5Dcreate(file_id, "Dark_Matter_Density", file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if ( DataUnits[field] == NULL )
      {
        DataUnits[field] = "none";
      }
 
      WriteStringAttr(dset_id, "Label", "Dark_Matter_Density", log_fptr);
      WriteStringAttr(dset_id, "Units", "", log_fptr);
      WriteStringAttr(dset_id, "Format", "e10.4", log_fptr);
      WriteStringAttr(dset_id, "Geometry", "Cartesian", log_fptr);
 
      h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Fclose(file_id);
        if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      /* Clean up if we modified the resolution. */
 
      if (SelfGravity && GravityResolution != 1)
	this->DeleteGravitatingMassFieldParticles();
 
    } // end of (if GravitatingMassFieldParticles != NULL)
 
    delete temp;
    for (dim = 0; dim < GridRank; dim++)
      delete tempdim[dim];
 
    /* Write BoundaryFluxes info (why? it's just recreated when the grid
                                  is read in) */
 
   }  // end: if (ProcessorNumber == MyProcessorNumber)
  } // end: if (NumberOfBaryonFields > 0)
 
 
//  DISABLE PARTICLE SORT & OUTPUT FOR NOW
 
if (MyProcessorNumber == ROOT_PROCESSOR)
  printf("Particle sort & output DISABLED\n");
 
if ( 0 == 1 )
{
 
  /* ------------------------------------------------------------------- */
  /* 3) Save particle quantities. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fprintf(fptr, "NumberOfParticles   = %"ISYM"\n", NumberOfParticles);
 
  if (NumberOfParticles > 0) {
 
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(fptr, "ParticleFileName = %s\n", name); // must be same as above
 
    if (MyProcessorNumber == ProcessorNumber) {
 
    /* Sort particles according to their identifier. */
 
    this->SortParticlesByNumber();
 
    /* Create a temporary buffer (32 bit). */
 
    float32 *temp = new float32[NumberOfParticles];
 
    /* Particle positions are not converted to 32 bit first.
       (128 bit numbers are not supported by HDF so convert to 64). */
 
    float64 *temp_pointer = NULL;
 
    if (sizeof(FLOAT) == 16)
      temp_pointer = new float64[NumberOfParticles];
 
    TempIntArray[0] = NumberOfParticles;
 
    for (dim = 0; dim < GridRank; dim++) {
 
      /* Convert to 64 if 128, either just write out. */
 
      if (sizeof(FLOAT) == 16)
	for (i = 0; i < NumberOfParticles; i++)
	  temp_pointer[i] = float64(ParticlePosition[dim][i]);
      else
	temp_pointer = (float64*) ParticlePosition[dim];
 
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Fcreate with Name = %s\n", ParticlePositionLabel[dim]);
 
      file_id = H5Fcreate(ParticlePositionLabel[dim], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
        if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", ParticlePositionLabel[dim]);
 
      dset_id =  H5Dcreate(file_id, ParticlePositionLabel[dim],  FILE_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dwrite(dset_id, FLOAT_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp_pointer);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Fclose(file_id);
        if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    }
 
    if (sizeof(FLOAT) == 16)
      delete [] temp_pointer;  /* clean up if allocated. */
 
    /* Copy particle velocities to temp and write them. */
 
    for (dim = 0; dim < GridRank; dim++) {
 
      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleVelocity[dim][i]);
 
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n",file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Fcreate with Name = %s\n", ParticleVelocityLabel[dim]);
 
      file_id = H5Fcreate(ParticleVelocityLabel[dim], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
        if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", ParticleVelocityLabel[dim]);
 
      dset_id = H5Dcreate(file_id, ParticleVelocityLabel[dim], file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Fclose(file_id);
        if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    }
 
    /* Copy mass to temp and write it. */
 
    for (i = 0; i < NumberOfParticles; i++)
      temp[i] = float32(ParticleMass[i]);
 
 
    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    if (io_log) fprintf(log_fptr, "H5Fcreate with Name = particle_mass\n");
 
    file_id = H5Fcreate("particle_mass", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
      if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    if (io_log) fprintf(log_fptr, "H5Dcreate with Name = particle_mass\n");
 
    dset_id = H5Dcreate(file_id, "particle_mass", file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Fclose(file_id);
      if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    /* Copy number (ID) to temp and write it. */
 
    PINT *tempPINT = new PINT[NumberOfParticles];
 
    for (i = 0; i < NumberOfParticles; i++)
      tempPINT[i] = ParticleNumber[i];
 
 
    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
    if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    if (io_log) fprintf(log_fptr, "H5Fcreate with Name = particle_index\n");
 
    file_id = H5Fcreate("particle_index", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
      if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    if (io_log) fprintf(log_fptr, "H5Dcreate with Name = particle_index\n");
 
    dset_id = H5Dcreate(file_id, "particle_index", HDF5_FILE_PINT, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Dwrite(dset_id, HDF5_PINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempPINT);
      if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Fclose(file_id);
      if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    /* Copy particle attributes to temp and write them. */
 
    for (j = 0; j < NumberOfParticleAttributes; j++) {
 
      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = float32(ParticleAttribute[j][i]);
 
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Fcreate with Name = %s\n", ParticleAttributeLabel[j]);
 
      file_id = H5Fcreate(ParticleAttributeLabel[j], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
        if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Dcreate with Name = %s\n", ParticleAttributeLabel[j]);
 
      dset_id = H5Dcreate(file_id, ParticleAttributeLabel[j], file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Fclose(file_id);
        if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    }
 
    /* clean up */
 
    delete temp;
    delete tempPINT;
 
  } // end: if (MyProcessorNumber...)
  } // end: if (NumberOfParticles > 0)
 
}
 
// PARTICLE OUTPUT DISABLED
 
  /* Close HDF file. */
 
 
//  if (MyProcessorNumber == ProcessorNumber)
//  {
//     h5_status = H5Fclose(file_id);
//       if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
//       if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
//  }
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    if (io_log) fclose(log_fptr);
  }
 
  /* 4) Save Gravity info. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    if (SelfGravity)

      fprintf(fptr, "GravityBoundaryType = %"ISYM"\n", GravityBoundaryType);
 
  /* Clean up. */
 
  delete name;
 
  return SUCCESS;
 
}
