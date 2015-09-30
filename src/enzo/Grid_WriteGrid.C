/***********************************************************************
/
/  GRID CLASS (WRITE OUT GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  Michael Kuhlen, October 2010, HDF5 hierarchy
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
#include <math.h>


 
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
 
int FindField(int field, int farray[], int numfields);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
#ifdef IO_64
#define io_type float64
#else
#define io_type float32
#endif

int WriteDataset(hid_t WriteLoc, float * data_buffer, io_type * tmp_buffer,
		 int * DataDims, int GridRank,
		 int *WriteStartIndex, int *WriteEndIndex, int * WriteDims,
		 char * Label, char * Units,hid_t file_type_id,hid_t float_type_id,FILE *log_fptr ) ;

int grid::WriteGrid(FILE *fptr, char *base_name, int grid_id)
{
  printf("In WriteGrid**");
 
  int i, j, k, dim, field, size, active_size, ActiveDim[MAX_DIMENSION];
  int file_status;
  float *temperature, *dust_temperature,
    *cooling_time;

 
  io_type *temp;
 
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
#ifdef WINDS
    char *ParticleAttributeLabel[] = 
      {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x", 
       "particle_jet_y", "particle_jet_z", "typeia_fraction"};
#else
    char *ParticleAttributeLabel[] = 
      {"creation_time", "dynamical_time", "metallicity_fraction", "typeia_fraction"};
#endif
  char *SmoothedDMLabel[] = {"Dark_Matter_Density", "Velocity_Dispersion",
			     "Particle_x-velocity", "Particle_y-velocity",
			     "Particle_z-velocity"};
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
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
 
  if (MyProcessorNumber == ROOT_PROCESSOR && HierarchyFileOutputFormat > 0) {

    fprintf(fptr, "Task              = %"ISYM"\n", ProcessorNumber);
 
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
 
  char *name = new char[MAX_LINE_LENGTH];
  strcpy(name, base_name);
  strcat(name, id);
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    char *logname = new char[MAX_LINE_LENGTH];
    strcpy(logname, base_name);
    strcat(logname, id);
    strcat(logname, ".log");
    if (io_log) log_fptr = fopen(logname, "a");
    delete [] logname;
 
    if (io_log) fprintf(log_fptr, "Grid_WriteGrid\n");
    if (io_log) fprintf(log_fptr, "  ID %"ISYM"  %s\n", grid_id, id);
    if (io_log) fprintf(log_fptr, "  BASE %s\n", base_name);
    if (io_log) fprintf(log_fptr, "  HDF %s\n", name);
  }
 
  /* Open HDF file for writing. */
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    if (io_log) fprintf(log_fptr,"H5Fcreate with Name = %s\n",name);
 
    file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fcreate id: %"ISYM"\n", file_id);
      if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    fprintf(stderr,"GPFS create: %s on CPU %"ISYM"\n", name, MyProcessorNumber);
  }
 

  /* ------------------------------------------------------------------- */
  /* 2) save baryon field quantities (including fields). */
 
  int ii = sizeof(io_type);
 
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
      printf("INCORRECT FLOAT DEFINITION\n");
 
  }
 
  if (NumberOfBaryonFields > 0) {
 
    if (MyProcessorNumber == ROOT_PROCESSOR && HierarchyFileOutputFormat > 0) {
 
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
 
    temp = new io_type[size];
 
    /* 2c) Loop over fields, writing each one. */
 
    for (field = 0; field < NumberOfBaryonFields; field++) {
 
      /* copy active part of field into grid */
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           +
	         (j-GridStartIndex[1])*ActiveDim[0]              +
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		       io_type(
	      BaryonField[field][i + j*GridDimension[0] +
		                     k*GridDimension[0]*GridDimension[1]]
                              );
 
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n",DataLabel[field]);
 
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
 
    }   // end of loop over fields

    if( UseMHDCT ){

      hsize_t MHDOutDims[3];
      int MHDActive[3], MHDWriteStartIndex[3], MHDWriteEndIndex[3];
      int BiggieSize = (GridDimension[0]+1)*(GridDimension[1]+1)*(GridDimension[2]+1);
      int index1, index2;
      io_type *MHDtmp = new io_type[BiggieSize];
      int WriteBoundary = FALSE; 

      for(field=0;field<nBfields;field++){
      if( WriteBoundary == TRUE){
	for(i=0;i<3;i++){
	  MHDWriteStartIndex[i] = 0;
	  MHDWriteEndIndex[i] = MagneticDims[field][i]-1;
	}
      }else{
	for(i=0;i<3;i++){
	  MHDWriteStartIndex[i] = MHDStartIndex[field][i];
	  MHDWriteEndIndex[i] = MHDEndIndex[field][i];
	}
      }
      /*
	int WriteDataset(hid_t WriteLoc, float * data_buffer, io_type * tmp_buffer,
	int * DataDims, int GridRank,
	int *WriteStartIndex, int *WriteEndIndex, int * WriteDims,
	char * Label, char * Units,hid_t file_type_id,hid_t float_type_id,FILE *log_fptr ) 
      */
      for (dim = 0; dim < 3; dim++)
	MHDActive[dim] = MHDWriteEndIndex[dim] - MHDWriteStartIndex[dim] +1;
      WriteDataset(file_id,MagneticField[field],MHDtmp,
		     MagneticDims[field],GridRank,
		     MHDWriteStartIndex,MHDWriteEndIndex,MHDActive,
		     MHDLabel[field],MHDUnits[0], file_type_id, float_type_id,log_fptr);
      }
      if( MHD_WriteElectric ){
	for(field=0;field<nBfields;field++){
	  if( WriteBoundary == TRUE ){
	    for( i=0;i<3;i++){
	      MHDWriteStartIndex[i] = 0;
	      MHDWriteEndIndex[i] = ElectricDims[field][i] - 1;
	    }
	  }else{
	    for(i=0;i<3;i++){
	      MHDWriteStartIndex[i] = MHDeStartIndex[field][i];
	      MHDWriteEndIndex[i] = MHDeEndIndex[field][i];
	    }
	  }
	  for(dim = 0; dim<3; dim++)
	    MHDActive[dim] = MHDWriteEndIndex[dim] - MHDWriteStartIndex[dim] +1;
	  
	  WriteDataset(file_id,ElectricField[field],MHDtmp,
		       ElectricDims[field],GridRank,
		       MHDWriteStartIndex,MHDWriteEndIndex,MHDActive,
		       MHDeLabel[field],MHDeUnits[0], file_type_id, float_type_id,log_fptr);

	}
	
      }//WriteElectric
      delete [] MHDtmp;
    }//UseMHDCT
    
    /* If requested compute and output the temperature field
       as well since its such a pain to compute after the fact. */
  
    if (OutputTemperature) {
 
      /* Allocate field and compute temperature. */
 
      temperature = new float[size];
 
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
		     io_type(
		   temperature[(k*GridDimension[1] + j)*GridDimension[0] + i]
			     );
 
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = Temperature\n");
 
      dset_id = H5Dcreate(file_id, "Temperature", file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
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
 
      // If outputing dust temperature, keep temperature field for the calculation.
      if (!OutputDustTemperature) {
	delete [] temperature;
      }
 
    } // end: if (OutputTemperature)


    /* If requested compute and output the dust temperature field
       as well since its such a pain to compute after the fact. */
  
    if (OutputDustTemperature != FALSE) {
 
      /* Get temperature field if we do not already have it. */

      if (!OutputTemperature) {
	temperature = new float[size];

	if (this->ComputeTemperatureField(temperature) == FAIL) {
	  ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");
	}
      }

      /* Allocate field and compute temperature. */
 
      dust_temperature = new float[size];
 
      if (this->ComputeDustTemperatureField(temperature,
					    dust_temperature) == FAIL) {
	ENZO_FAIL("Error in grid->ComputeDustTemperatureField.\n");
      }
 
      /* Copy active part of field into grid */
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           +
	         (j-GridStartIndex[1])*ActiveDim[0]              +
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     io_type(
		   temperature[(k*GridDimension[1] + j)*GridDimension[0] + i]
			     );
 
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = Dust_Temperature\n");
 
      dset_id = H5Dcreate(file_id, "Dust_Temperature", file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      WriteStringAttr(dset_id, "Label", "Dust_Temperature", log_fptr);
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
 
      if (!OutputTemperature) {
	delete [] temperature;
      }
      delete [] dust_temperature;
 
    } // end: if (OutputDustTemperature)
 
    if (VelAnyl) {

     
      int Vel1Num, Vel2Num, Vel3Num;
      
      for (int i=0; i< NumberOfBaryonFields; i++){


	switch (FieldType[i]){
	case Velocity1:
	  Vel1Num=i; break;
	case Velocity2:
	  Vel2Num=i; break;
	case Velocity3:
	  Vel3Num=i; break;
	}
      }


      io_type *curlx; io_type *curly;

      
      if(GridRank==3){
      curlx = new io_type [size];
      curly = new io_type [size];}


      io_type *curlz = new io_type [size];
      io_type *div   = new io_type [size];

      FLOAT dx = CellWidth[0][0],
	dy = CellWidth[1][0], dz;
      if (GridRank>2)
	dz = CellWidth[2][0];
    
      /* Copy active part of field into grid */
      int igrid, igridyp1, igridym1, igridzp1, igridzm1;
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {

	    igrid = (k*GridDimension[1] + j)*GridDimension[0] + i;
	    igridyp1 = (k*GridDimension[1] + j + 1)*GridDimension[0] + i;
	    igridym1 = (k*GridDimension[1] + j - 1)*GridDimension[0] + i;
	    igridzp1 = ((k+1)*GridDimension[1]+j)*GridDimension[0] + i;
	    igridzm1 = ((k-1)*GridDimension[1]+j)*GridDimension[0] + i;



	    if (GridRank==3){
	      
	    div[(i-GridStartIndex[0])+
		(j-GridStartIndex[1])*ActiveDim[0]+
		(k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
              io_type(
		      (0.5*(BaryonField[Vel1Num][igrid+1]-BaryonField[Vel1Num][igrid-1])/dx +
		      0.5*(BaryonField[Vel2Num][igridyp1]-BaryonField[Vel2Num][igridym1])/dy +
			  0.5*(BaryonField[Vel3Num][igridzp1]-BaryonField[Vel3Num][igridzm1])/dz)
                     );
	    curlx[(i-GridStartIndex[0])+
	       (j-GridStartIndex[1])*ActiveDim[0]+ 
	       (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
	       io_type(
		      (0.5*(BaryonField[Vel3Num][igridyp1]-BaryonField[Vel3Num][igridym1])/dy -		      
			  0.5*(BaryonField[Vel2Num][igridzp1]-BaryonField[Vel2Num][igridzm1])/dz)
		      );
	    curly[(i-GridStartIndex[0])+
	       (j-GridStartIndex[1])*ActiveDim[0]+ 
	       (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
	      io_type(
		       (0.5*(BaryonField[Vel1Num][igridzp1]-BaryonField[Vel1Num][igridzm1])/dz -		      
			   0.5*(BaryonField[Vel3Num][igrid+1]-BaryonField[Vel3Num][igrid-1])/dx)
		       );
	     curlz[(i-GridStartIndex[0])+
		(j-GridStartIndex[1])*ActiveDim[0]+ 
		(k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
	      io_type(
		      (0.5*(BaryonField[Vel2Num][igrid+1]-BaryonField[Vel2Num][igrid-1])/dx -		      
			  0.5*(BaryonField[Vel1Num][igridyp1]-BaryonField[Vel1Num][igridym1])/dy)
		      );
	    }

	    if (GridRank==2){
	      
	    div[(i-GridStartIndex[0])+
		(j-GridStartIndex[1])*ActiveDim[0]+
		(k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
              io_type(
		      (0.5*(BaryonField[Vel1Num][igrid+1]-BaryonField[Vel1Num][igrid-1])/dx +
		       0.5*(BaryonField[Vel2Num][igridyp1]-BaryonField[Vel2Num][igridym1])/dy 
		       ));
	     curlz[(i-GridStartIndex[0])+
		(j-GridStartIndex[1])*ActiveDim[0]+ 
		(k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
	      io_type(
		      (0.5*(BaryonField[Vel2Num][igrid+1]-BaryonField[Vel2Num][igrid-1])/dx -		      
			  0.5*(BaryonField[Vel1Num][igridyp1]-BaryonField[Vel1Num][igridym1])/dy)
		      );
	    }

	    

	  }
	}
      }

      /* output */
      
      char *DataLabelN[4];
      if (GridRank==2) {
	DataLabelN[0]="Velocity_Div";
	DataLabelN[1]="Velocity_Vorticity";
      }
      if (GridRank==3) {
	DataLabelN[0]="Velocity_Div";
	DataLabelN[1]="Velocity_Vorticity1";
	DataLabelN[2]="Velocity_Vorticity2";
	DataLabelN[3]="Velocity_Vorticity3";
      }

      

      int tFields=3;
      if (GridRank==2) tFields=1;


      for (int field=0; field<=tFields; field++){

      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n",DataLabelN[field]);
 
      dset_id =  H5Dcreate(file_id, DataLabelN[field], file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      /* set datafield name and units, etc. */
 
      // fprintf(stderr,  DataLabelN[field]);

      WriteStringAttr(dset_id, "Label", DataLabelN[field], log_fptr);
    
      WriteStringAttr(dset_id, "Units", "VortUnits", log_fptr);
    
      WriteStringAttr(dset_id, "Format", "e10.4", log_fptr);
    
      WriteStringAttr(dset_id, "Geometry", "Cartesian", log_fptr);
    

      switch (field) {
      case 0:
	temp=div;
      case 1:
	temp=curlz;
      case 2:
	temp=curly;
      case 3:
	temp=curlx;

      }
 
      h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}


      }

      fprintf(stderr,"Inside sorta \n");


      if(GridRank==3){
	delete curlx;
	delete curly;}
      delete curlz;
      delete div;
    }



    if (OutputCoolingTime != FALSE) {
 
      /* Allocate field and compute cooling time. */

      cooling_time = new float[size];
 
      float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
	VelocityUnits = 1, TimeUnits = 1, aUnits = 1;

      GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time);

      if (this->ComputeCoolingTime(cooling_time) == FAIL) {
	ENZO_FAIL("Error in grid->ComputeCoolingTime.\n");
      }

      // Make all cooling time values positive and convert to seconds.
      for (i = 0;i < size;i++) {
	cooling_time[i] = fabs(cooling_time[i]) * TimeUnits;
      }
 
      /* Copy active part of field into grid */
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           +
	         (j-GridStartIndex[1])*ActiveDim[0]              +
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		     io_type(
		   cooling_time[(k*GridDimension[1] + j)*GridDimension[0] + i]
			     );
 
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = Cooling_Time\n");
 
      dset_id = H5Dcreate(file_id, "Cooling_Time", file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      WriteStringAttr(dset_id, "Label", "Cooling_Time", log_fptr);
      WriteStringAttr(dset_id, "Units", "s", log_fptr);
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
 
      delete [] cooling_time;
 
    } // if (OutputCoolingTime)
   
 
    /* Make sure that there is a copy of dark matter field to save
       (and at the right resolution). */

    if (OutputSmoothedDarkMatter == FALSE) {
    
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
//      fprintf(stderr, "%"ISYM" %"ISYM" %10.4"FSYM" %10.4"FSYM" %10.4"FSYM" %10.4"FSYM"\n", StartIndex[dim], EndIndex[dim], GridLeftEdge[dim], GridRightEdge[dim], GravitatingMassFieldParticlesLeftEdge[dim], GravitatingMassFieldParticlesCellSize);
	}
 
	/* Copy active part of field into grid */
 
	for (k = StartIndex[2]; k <= EndIndex[2]; k++)
	  for (j = StartIndex[1]; j <= EndIndex[1]; j++)
	    for (i = StartIndex[0]; i <= EndIndex[0]; i++)
	      temp[(i-StartIndex[0])                           +
		   (j-StartIndex[1])*ActiveDim[0]              +
		   (k-StartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
		io_type(
			GravitatingMassFieldParticles[ i +
			j*GravitatingMassFieldParticlesDimension[0] +
			k*GravitatingMassFieldParticlesDimension[0]*
			GravitatingMassFieldParticlesDimension[1]]
			);
 
	file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	if (io_log) fprintf(log_fptr,"H5Dcreate with Name = Cooling_Time\n");
 
	dset_id = H5Dcreate(file_id, "Dark_Matter_Density", file_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	WriteStringAttr(dset_id, "Label", "Cooling_Time", log_fptr);
	WriteStringAttr(dset_id, "Units", "s", log_fptr);
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
  
	/* Clean up if we modified the resolution. */
 
	if (SelfGravity && GravityResolution != 1)
	  this->DeleteGravitatingMassFieldParticles();
 
      } // end of (if GravitatingMassFieldParticles != NULL)

    } // ENDIF !OutputSmoothedDarkMatter

 
    delete [] temp;
 
    for (dim = 0; dim < GridRank; dim++)
      delete [] tempdim[dim];
 
    /* Write BoundaryFluxes info (why? it's just recreated when the grid
                                  is read in) */
 
   }  // end: if (ProcessorNumber == MyProcessorNumber)
  } // end: if (NumberOfBaryonFields > 0)

  /* ------------------------------------------------------------------- */
  /* 2b) Save particle quantities smoothed to the grid. */
 
  if (OutputSmoothedDarkMatter > 0 && MyProcessorNumber == ProcessorNumber) {

    size = active_size = 1;
    for (dim = 0; dim < GridRank; dim++) {
      OutDims[GridRank-dim-1] = ActiveDim[dim];
      size *= GridDimension[dim];
      active_size *= ActiveDim[dim];
    }
 
    temp = new io_type[size];

    int NumberOfDMFields;
    switch (OutputSmoothedDarkMatter) {
    case 1: NumberOfDMFields = 1; break;  // density
    case 2: NumberOfDMFields = 5; break;  // + rms velocity + 3-velocity
    } // ENDSWITCH
      
    for (field = 0; field < NumberOfDMFields; field++) {

      // Only the active part was calculated, so just copy over.
      for (i = 0; i < active_size; i++)
	temp[i] = io_type(InterpolatedField[field][i]);
	
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n", SmoothedDMLabel[field]);
 
      dset_id = H5Dcreate(file_id, SmoothedDMLabel[field], file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      WriteStringAttr(dset_id, "Label", SmoothedDMLabel[field], log_fptr);
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
 
      delete [] InterpolatedField[field];
      InterpolatedField[field] = NULL;

    } // ENDFOR field

    delete [] temp;
      
  } // ENDIF OutputSmoothedDarkMatter
 
  /* ------------------------------------------------------------------- */
  /* 3) Save particle quantities. */
 
  if (MyProcessorNumber == ROOT_PROCESSOR && HierarchyFileOutputFormat > 0)
    fprintf(fptr, "NumberOfParticles   = %"ISYM"\n", NumberOfParticles);
 
  if (NumberOfParticles > 0) {
 
    if (MyProcessorNumber == ROOT_PROCESSOR && HierarchyFileOutputFormat > 0)
      fprintf(fptr, "ParticleFileName = %s\n", name); // must be same as above
 
    if (MyProcessorNumber == ProcessorNumber) {
 
    /* Sort particles according to their identifier. */
 
    this->SortParticlesByNumber();
 
    /* Create a temporary buffer (64 bit). */
 
    io_type *temp = new io_type[NumberOfParticles];
 
    /* Particle positions are not converted to 32 bit first.
       (128 bit numbers are not supported by HDF so convert to 64). */

    io_type *temp_pointer = NULL;
    float128 *long_temp_pointer = NULL;
 
    TempIntArray[0] = NumberOfParticles;
 
    for (dim = 0; dim < GridRank; dim++) {
 
      /* Convert to 64 if 128, either just write out. */

      if (sizeof(FLOAT) == 16) {
        long_temp_pointer = (float128*) ParticlePosition[dim];
      }
      else
      {
        temp_pointer = (io_type*) ParticlePosition[dim];
      }
 
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n", ParticlePositionLabel[dim]);
 
      dset_id =  H5Dcreate(file_id, ParticlePositionLabel[dim],  FILE_type_id, file_dsp_id, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (sizeof(FLOAT) == 16) {
//                                  NOTE: for 128bits this must be FILE_type_id and NOT FLOAT_type_id!
      h5_status = H5Dwrite(dset_id, FILE_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) long_temp_pointer);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

      }
      else
      {

      h5_status = H5Dwrite(dset_id, FLOAT_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp_pointer);
        if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

      }
 
      h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    }
 
    /* Copy particle velocities to temp and write them. */
 
    for (dim = 0; dim < GridRank; dim++) {
 
      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = io_type(ParticleVelocity[dim][i]);
 
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n",ParticleVelocityLabel[dim]);
 
      dset_id =  H5Dcreate(file_id, ParticleVelocityLabel[dim], file_type_id, file_dsp_id, H5P_DEFAULT);
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
 
    }
 
    /* Copy mass to temp and write it. */
 
    for (i = 0; i < NumberOfParticles; i++)
      temp[i] = io_type(ParticleMass[i]);
 
 
    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    if (io_log) fprintf(log_fptr,"H5Dcreate with Name = particle_mass\n");
 
    dset_id =  H5Dcreate(file_id, "particle_mass", file_type_id, file_dsp_id, H5P_DEFAULT);
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
 
    /* Copy number (ID) to temp and write it. */
 
    PINT *tempPINT = new PINT[NumberOfParticles];
 
    for (i = 0; i < NumberOfParticles; i++)
      tempPINT[i] = ParticleNumber[i];
 
 
    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    if (io_log) fprintf(log_fptr,"H5Dcreate with Name = particle_index\n");
 
    dset_id =  H5Dcreate(file_id, "particle_index", HDF5_FILE_PINT, file_dsp_id, H5P_DEFAULT);
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
 
 
    /* Copy type to temp and write it. */

    if (ParticleTypeInFile == TRUE) {

    int *tempint = new int[NumberOfParticles];
    
    if (ParticleType == NULL)
      ENZO_FAIL("Undefined ParticleType!\n");
 
    for (i = 0; i < NumberOfParticles; i++)
      tempint[i] = ParticleType[i];
 
    file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    if (io_log) fprintf(log_fptr,"H5Dcreate with Name = particle_type\n");
 
    dset_id =  H5Dcreate(file_id, "particle_type", HDF5_FILE_INT, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Dcreate id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Dwrite(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempint);
      if (io_log) fprintf(log_fptr, "H5Dwrite: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
    h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

    delete [] tempint;
 
    }

 
    /* Copy particle attributes to temp and write them. */
 
    for (j = 0; j < NumberOfParticleAttributes; j++) {
 
      for (i = 0; i < NumberOfParticles; i++)
	temp[i] = io_type(ParticleAttribute[j][i]);
 
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dcreate with Name = %s\n",ParticleAttributeLabel[j]);
 
      dset_id =  H5Dcreate(file_id, ParticleAttributeLabel[j], file_type_id, file_dsp_id, H5P_DEFAULT);
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
 
    }
 
 
    /* clean up */
 
    delete [] temp;
    delete [] tempPINT;
 
  } // end: if (MyProcessorNumber...)
  } // end: if (NumberOfParticles > 0)

 
 
  /* Close HDF file. */
 
  if (MyProcessorNumber == ProcessorNumber)
  {
 
     h5_status = H5Fclose(file_id);
   
     if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
     if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
   
  }
 
  if (MyProcessorNumber == ProcessorNumber)
  {  
    if (io_log) fclose(log_fptr);
   
  }
 
  /* 4) Save Gravity info. */


 
 
  if (MyProcessorNumber == ROOT_PROCESSOR && HierarchyFileOutputFormat > 0)
    if (SelfGravity)

      fprintf(fptr, "GravityBoundaryType = %"ISYM"\n", GravityBoundaryType);
 
  /* Clean up. */
 
  
  delete [] name;

  
 
  return SUCCESS;
 
}


