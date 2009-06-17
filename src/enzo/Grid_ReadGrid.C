/***********************************************************************
/
/  GRID CLASS (READ GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  Alexei Kritsuk, Jan 2004   a trick for RandomForcing //AK
/
/  PURPOSE:
/
************************************************************************/
 
//  Input a grid from file pointer fpt

#ifdef USE_HDF4
#include "mfhdf.h"
#endif /* USE_HDF4 */
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
 
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif
 
// function prototypes
 
int ReadListOfFloats(FILE *fptr, int N, FLOAT floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
 
// extern int ParticleTypeInFile; // declared and set in ReadParameterFile
 
#ifdef USE_HDF4
int ReadField(float *temp, int Dims[], int Rank, char *name,
	      char *field_name);
static int32 sd_id, sds_index; // HDF4 (SD) handlers                                               
#endif 
 
int grid::ReadGrid(FILE *fptr, int GridID, 
		   int ReadText, int ReadData)
{
 
  int i, j, k, dim, field, size, active_size;
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
  char logname[MAX_LINE_LENGTH];
 
  int ActiveDim[MAX_DIMENSION];
 
  FILE *log_fptr;
 
  hid_t       file_id, dset_id;
  hid_t       float_type_id, FLOAT_type_id;
  hid_t       file_type_id, FILE_type_id;
  hid_t       file_dsp_id;
  hid_t       num_type;
 
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     TempIntArray[MAX_DIMENSION];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int         num_size;
 
  char *ParticlePositionLabel[] =
    {"particle_position_x", "particle_position_y", "particle_position_z"};
  char *ParticleVelocityLabel[] =
    {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
  char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
                                    "metallicity_fraction", "alpha_fraction"};

#ifdef USE_HDF4
    int32 TempIntArray2[MAX_DIMENSION];
    int32 sds_id, num_type2, attributes, TempInt;  
    sds_index = 0;  // start at first SDS                                                                            
#endif 
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

#ifdef IO_64
#define io_type float64
#else
#define io_type float32
#endif
 
  if(ReadText){

    /* Read general grid class data */

    /* make sure quantities defined at least for 3d */
 
    for (dim = GridRank; dim < 3; dim++) {
      GridDimension[dim] = 1;
      GridStartIndex[dim] = 0;
      GridEndIndex[dim] = 0;
    }

    if (fscanf(fptr, "GridRank = %"ISYM"\n", &GridRank) != 1) {
      fprintf(stderr, "Error reading GridRank.\n");
      ENZO_FAIL("");
    }
 
    if (fscanf(fptr, "GridDimension = ") != 0) {
      fprintf(stderr, "Error reading GridDimension(0).\n");
      ENZO_FAIL("");
    }
 
    if (ReadListOfInts(fptr, GridRank, GridDimension) == FAIL) {
      fprintf(stderr, "Error reading GridDimension(1).\n");
      ENZO_FAIL("");
    }
 
    fscanf(fptr, "GridStartIndex = ");
 
    if (ReadListOfInts(fptr, GridRank, GridStartIndex) == FAIL) {
      fprintf(stderr, "Error reading GridStartIndex.\n");
      ENZO_FAIL("");
    }
 
    fscanf(fptr, "GridEndIndex = ");
 
    if (ReadListOfInts(fptr, GridRank, GridEndIndex) == FAIL) {
      fprintf(stderr, "Error reading GridEndIndex.\n");
      ENZO_FAIL("");
    }
 
    fscanf(fptr, "GridLeftEdge = ");
 
    if (ReadListOfFloats(fptr, GridRank, GridLeftEdge) == FAIL) {
      fprintf(stderr, "Error reading GridLeftEdge.\n");
      ENZO_FAIL("");
    }
 
    fscanf(fptr, "GridRightEdge = ");
 
    if (ReadListOfFloats(fptr, GridRank, GridRightEdge) == FAIL) {
      fprintf(stderr, "Error reading GridRightEdge.\n");
      ENZO_FAIL("");
    }

    if (fscanf(fptr, "Time = %"PSYM"\n", &Time) != 1) {
      fprintf(stderr, "Error reading Time.\n");
      ENZO_FAIL("");
    }
 
    if (fscanf(fptr, "SubgridsAreStatic = %"ISYM"\n", &SubgridsAreStatic) != 1) {
      fprintf(stderr, "Error reading SubgridsAreStatic.\n");
      ENZO_FAIL("");
    }

    /* Read baryon field quantities. */
 
    if (fscanf(fptr, "NumberOfBaryonFields = %"ISYM"\n",
	       &NumberOfBaryonFields) != 1) {
      fprintf(stderr, "Error reading NumberOfBaryonFields.\n");
      ENZO_FAIL("");
    }
    if (NumberOfBaryonFields > 0) {
 
      fscanf(fptr, "FieldType = ");
 
      if (ReadListOfInts(fptr, NumberOfBaryonFields, FieldType) == FAIL) {
	fprintf(stderr, "Error reading FieldType.\n");
	ENZO_FAIL("");
      }
 
      fgetpos(fptr, &BaryonFileNamePosition); //AK
 
      if (fscanf(fptr, "BaryonFileName = %s\n", name) != 1) {
	fprintf(stderr, "Error reading BaryonFileName.\n");
	ENZO_FAIL("");
      }
 
      fscanf(fptr, "CourantSafetyNumber    = %"FSYM"\n", &CourantSafetyNumber);
      fscanf(fptr, "PPMFlatteningParameter = %"ISYM"\n", &PPMFlatteningParameter);
      fscanf(fptr, "PPMDiffusionParameter  = %"ISYM"\n", &PPMDiffusionParameter);
      fscanf(fptr, "PPMSteepeningParameter = %"ISYM"\n", &PPMSteepeningParameter);
    }

    /* 3) Read particle info */
 
    if (fscanf(fptr, "NumberOfParticles = %"ISYM"\n", &NumberOfParticles) != 1) {
      fprintf(stderr, "error reading NumberOfParticles.\n");
      ENZO_FAIL("");
    }
  
    if (NumberOfParticles > 0) {
 
      /* Read particle file name. */

      if (fscanf(fptr, "ParticleFileName = %s\n", name) != 1) {
	fprintf(stderr, "Error reading ParticleFileName.\n");
	ENZO_FAIL("");
      }
    }

    /* 4) Read gravity info */
 
    if (SelfGravity)
      if (fscanf(fptr, "GravityBoundaryType = %"ISYM"\n",&GravityBoundaryType) != 1) {
	fprintf(stderr, "Error reading GravityBoundaryType.\n");
	ENZO_FAIL("");
      }
  }

  /* Compute Flux quantities */

  this->PrepareGridDerivedQuantities();

  if (HydroMethod == MHD_RK) {

    int activesize = 1;
    for (int dim = 0; dim < GridRank; dim++)
      activesize *= (GridDimension[dim]-2*DEFAULT_GHOST_ZONES);
    
    if (divB == NULL) 
      divB = new float[activesize];

    for (int dim = 0; dim < 3; dim++)
      if (gradPhi[dim] == NULL)
	gradPhi[dim] = new float[activesize];

    for (int dim = GridRank; dim < 3; dim++)
      for (int n = 0; n < activesize; n++)
	gradPhi[dim][n] = 0.0;

  } /* if HydroMethod == MHD */

 

  int ii = sizeof(io_type);
 
  switch(ii)
    {
 
    case 4:
      float_type_id = HDF5_R4;
      break;
 
    case 8:
      float_type_id = HDF5_R8;
      break;
 
    default:
      float_type_id = HDF5_R4;
 
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

  if (NumberOfBaryonFields > 0 && ReadData ) {
    if (MyProcessorNumber == ProcessorNumber)
      {
	strcpy(logname, name);
	strcat(logname, ".hdf.log");
	if (io_log) log_fptr = fopen(logname, "a");
      }
 
    if (MyProcessorNumber == ProcessorNumber){

#ifdef USE_HDF4
      if ((sd_id = SDstart(name, DFACC_RDONLY)) == HDF_FAIL) {
	fprintf(stderr, "Error opening file %s.\n", name);
	return FAIL;
      }

      sds_id = SDselect(sd_id, sds_index++);
      while (SDiscoordvar(sds_id)) {
	SDendaccess(sds_id);
	sds_id = SDselect(sd_id, sds_index++);
      }
      if (SDgetinfo(sds_id, dummy, &TempInt, TempIntArray2, &num_type2, 
		    &attributes) == HDF_FAIL) {
	fprintf(stderr, "Error reading dims0 from %s.\n", name);
	return FAIL;
      }
      SDendaccess(sds_id);
      sds_index--;

      if (TempInt != GridRank) {
	fprintf(stderr, "HDF rank (%d) does not match GridRank.\n", TempInt);
	return FAIL;
      }
#else
      if(!ReadText){
	// build filename from grid id
	char id[MAX_GRID_TAG_SIZE];
	sprintf(id, "%"GRID_TAG_FORMAT""ISYM, GridID);
	strcpy(name, PrevParameterFileName);
	strcat(name, ".grid");
	strcat(name, id);
      }
      
      if (io_log) fprintf(log_fptr,"H5Fopen with Name %s\n",name);
 
      file_id = H5Fopen(name,  H5F_ACC_RDONLY, H5P_DEFAULT);
      if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
      if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
#endif

      /* fill in ActiveDim for dims up to 3d */
 
      for (dim = 0; dim < 3; dim++)
	ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;
 
      /* check dimensions of HDF file against this grid
	 (note: we don't bother to check the coordinate arrays)  */
 
      size = 1;
      active_size = 1;
 
      for (dim = 0; dim < GridRank; dim++) {
	size *= GridDimension[dim];
	active_size *= ActiveDim[dim];
      }
 
      //  CAUTION - are the coordinates reversed?
 
      for (dim = 0; dim < GridRank; dim++) {
	OutDims[GridRank-dim-1] = ActiveDim[dim];
	if (io_log) fprintf(log_fptr, "Outdims %"ISYM"\n", (int) OutDims[GridRank-dim-1]);
      }
 
      /* allocate temporary space */
 
#ifdef USE_HDF4
      float *temp = new float[active_size];
#else
      io_type *temp = new io_type[active_size];
#endif
 
      /* loop over fields, reading each one */
 
      for (field = 0; field < NumberOfBaryonFields; field++) {
 
	/* get data into temporary array */

#ifdef USE_HDF4
	if (ReadField(temp, ActiveDim, GridRank, name, 
		      DataLabel[field]) == FAIL) {
	  fprintf(stderr, "Error reading field %d.\n", field);
	  return FAIL;
	}
#else
	file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", DataLabel[field]);
 
	dset_id =  H5Dopen(file_id, DataLabel[field]);
        if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
	h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
#endif 
 
	/* copy active region into whole grid */
 
	BaryonField[field] = new float[size];
 
	for (i = 0; i < size; i++)
	  BaryonField[field][i] = 0;
 
	for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	  for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	    for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	      BaryonField[field][i + j*GridDimension[0] +
				 k*GridDimension[0]*GridDimension[1]] =
		float(temp[(i-GridStartIndex[0])                         +
			   (j-GridStartIndex[1])*ActiveDim[0]            +
			   (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ]);

      } // end: loop over fields
 
      delete [] temp;
 
    }  // end: if (MyProcessorNumber == ProcessorNumber)
  }  // end read baryon fields
 
  /* 3) Read particle info */
 
  if (NumberOfParticles > 0 && ReadData) {
 
    if (MyProcessorNumber == ProcessorNumber) {

#ifdef USE_HDF4
#else
      if(!ReadText){
	// build filename from grid id
	char id[MAX_GRID_TAG_SIZE];
	sprintf(id, "%"GRID_TAG_FORMAT""ISYM, GridID);
	strcpy(name, PrevParameterFileName);
	strcpy(name, ".grid");
	strcat(name, id);
      }
#endif
 
      /* Open file if not already done (note: particle name must = grid name). */
 
      if (NumberOfBaryonFields == 0) {
 
#ifdef USE_HDF4
	if ((sd_id = SDstart(name, DFACC_RDONLY)) == HDF_FAIL) {
	  fprintf(stderr, "Error opening file %s.\n", name);
	  return FAIL;
	}
#else
	if (MyProcessorNumber == ProcessorNumber)
	  {
	    strcpy(logname, name);
	    strcat(logname, ".hdf.log");
	    if (io_log) log_fptr = fopen(logname, "a");
	  }
 
	if (io_log) fprintf(log_fptr, "H5Fopen with Name %s\n", name);
 
	file_id = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
        if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
#endif
 
      } // end: if (NumberOfBaryonFields == 0)
 
      /* Allocate room for particles. */
 
      this->AllocateNewParticles(NumberOfParticles);

      TempIntArray[0] = NumberOfParticles;
      
#ifdef USE_HDF4

      float *temp = new float[active_size];

      int32 HDFDataType = (sizeof(Eflt) == 4) ? DFNT_FLOAT32 : DFNT_FLOAT64;
      if (sizeof(Eflt) == 16) HDFDataType = DFNT_FLOAT128;

      /* Read dims.
	 If Rank != 1, we may have just read some other field SDS.  If so,
	 then try again. */
      
      TempInt = 0;
      while (TempInt != 1) {
	sds_id = SDselect(sd_id, sds_index++);
	while (SDiscoordvar(sds_id)) {
	  SDendaccess(sds_id);
	  sds_id = SDselect(sd_id, sds_index++);
	}
	if (SDgetinfo(sds_id, dummy, &TempInt, TempIntArray2, &num_type2, 
		      &attributes) == HDF_FAIL) {
	  fprintf(stderr, "Error reading dims1 from %s.\n", name);
	  return FAIL;
	}
	SDendaccess(sds_id);
      }
      sds_index--; 

      /* Check dims. */

      if (TempInt != 1 || TempIntArray2[0] != NumberOfParticles) {
	fprintf(stderr, "HDF particle dims do not match NumberOfParticles.\n");
	fprintf(stderr, "  (HDF dim[0] = %d, NumberOfParticles = %d)\n",
		int(TempIntArray2[0]), NumberOfParticles);
	return FAIL;  
      }

      /* Read ParticlePosition (use temporary buffer). */ 
      
      for (dim = 0; dim < GridRank; dim++) {

	if (num_type2 == HDFDataType) {

	  /* same data type: just read. */

	  if (ReadField((float *) ParticlePosition[dim], &NumberOfParticles, 1, name, 
			ParticlePositionLabel[dim]) == FAIL) {
	    fprintf(stderr, "Error reading ParticlePosition %d\n", dim);
	    return FAIL;
	  }
	
	} else {

	  /* convert data: Read into temporary buffer and copy. */

	  if (ReadField(temp, &NumberOfParticles, 1, name, 
			ParticlePositionLabel[dim]) == FAIL) {
	    fprintf(stderr, "Error reading ParticlePosition %d\n", dim);
	    return FAIL;
	  }

	  float64 *temp64 = (float64 *) temp;
	  long_double *temp128 = (long_double *) temp;
	  
	  if (num_type2 == DFNT_FLOAT32)
	    for (i = 0; i < NumberOfParticles; i++)
	      ParticlePosition[dim][i] = Eflt(temp[i]);
	  if (num_type2 == DFNT_FLOAT64)
	    for (i = 0; i < NumberOfParticles; i++)
	      ParticlePosition[dim][i] = Eflt(temp64[i]);
	  if (num_type2 == DFNT_FLOAT128)
	    for (i = 0; i < NumberOfParticles; i++)
	      ParticlePosition[dim][i] = Eflt(temp128[i]);
	}
      } // end: loop over dims

      delete [] temp;

      /* Read ParticleVelocity. */

      for (dim = 0; dim < GridRank; dim++) {
	if (ReadField(ParticleVelocity[dim], &NumberOfParticles, 1, name,
		      ParticleVelocityLabel[dim]) == FAIL) {
	  fprintf(stderr, "Error reading ParticleVelocity %d\n", dim);
	  return FAIL;
	}
      }

      /* Read ParticleMass. */

      if (ReadField(ParticleMass, &NumberOfParticles, 1, name,
		    "particle_mass") == FAIL)
	return FAIL;

      /* Read ParticleNumber */
      
      if (ReadField((float *) ParticleNumber, &NumberOfParticles, 1, name,
		    "particle_index") == FAIL)
	return FAIL;

      /* Read particle type if present */
      
      if (ParticleTypeInFile == TRUE) {
	
	if (ReadField((float *) ParticleType, &NumberOfParticles, 1, name,
		      "particle_type") == FAIL)
	  return FAIL;

#define NO_CHECK_PARTICLE_TYPE
#ifdef CHECK_PARTICLE_TYPE
	for (i = 0; i < NumberOfParticles; i++)
	  if (ParticleType[i] < PARTICLE_TYPE_GAS ||
	      ParticleType[i] > NUM_PARTICLE_TYPES-1) {
	  fprintf(stderr, "file: %s: particle %d has unknown type %d\n",
		  name, i, ParticleType[i]);
	  return FAIL;
	  }
#endif

    }

      /* Read ParticleAttributes. */

#define NO_RESTART_WITH_ATTRIBUTES
      for (j = 0; j < NumberOfParticleAttributes; j++) {
#ifdef RESTART_WITH_ATTRIBUTES
	for (i=0; i < NumberOfParticles; i++)
	  ParticleAttribute[j][i] = 0;
#else
	if (ReadField(ParticleAttribute[j], &NumberOfParticles, 1, name,
		      ParticleAttributeLabel[j]) == FAIL) {
	  fprintf(stderr, "Error reading ParticleAttribute %d\n", j);
	  return FAIL;
	}
#endif
      }
      
      /* If the particle type is not in the file, then set it according
	 to the value of the attributes. */

      if (ParticleTypeInFile != TRUE)
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleType[i] = ReturnParticleType(i);

#else /* USE_HDF4 */

      /* Create a temporary buffer (32 bit or twice the size for 64). */
      
      io_type *temp = NULL;
      
      jj = sizeof(FLOAT);
 
      switch(jj)
	{
	  
	case 4:
	  temp = new io_type[NumberOfParticles];
	  break;
 
	case 8:
	  temp = new io_type[NumberOfParticles*2];
	  break;
 
	case 16:
	  temp = new io_type[NumberOfParticles*4];
	  break;
 
	default:
	  printf("INCORRECT FLOAT DEFINITION\n");
 
	}
 
      if (temp == NULL)
	temp = new io_type[NumberOfParticles];
  
      /* Read ParticlePosition (use temporary buffer). */
 
      for (dim = 0; dim < GridRank; dim++) {
 
	file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	if (io_log) fprintf(log_fptr,"H5Dopen with Name = %s\n", ParticlePositionLabel[dim]);
 
	dset_id =  H5Dopen(file_id, ParticlePositionLabel[dim]);
        if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	num_type = H5Dget_type(dset_id);
	num_size = H5Tget_size(num_type);
 
	if (sizeof(FLOAT) == 16)
	  {
 
	    //                                 NOTE: for 128bits this must be FILE_type_id and NOT FLOAT_type_id!
	    h5_status = H5Dread(dset_id, FILE_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) ParticlePosition[dim]);
	    if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	  }
	else
	  {
 
	    h5_status = H5Dread(dset_id, FLOAT_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) ParticlePosition[dim]);
	    if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	  }
 
	h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      }
 
      /* Read ParticleVelocity. */
 
      for (dim = 0; dim < GridRank; dim++) {
 
	file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", ParticleVelocityLabel[dim]);
 
	dset_id =  H5Dopen(file_id, ParticleVelocityLabel[dim]);
        if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
        if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleVelocity[dim][i] = float(temp[i]);
      }
  
      /* Read ParticleMass into temporary buffer and Copy to ParticleMass. */
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dopen with Name = particle_mass\n");
 
      dset_id =  H5Dopen(file_id, "particle_mass");
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      for (i = 0; i < NumberOfParticles; i++)
	ParticleMass[i] = float(temp[i]);
 
      /* Read ParticleNumber into temporary buffer and Copy to ParticleNumber. */
 
      int *tempint = new int[NumberOfParticles];
 
      file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr,"H5Dopen  with Name = particle_index\n");
 
      dset_id =  H5Dopen(file_id, "particle_index");
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempint);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      for (i = 0; i < NumberOfParticles; i++)
	ParticleNumber[i] = tempint[i];
 
       // Read ParticleType if present
 
      if (ParticleTypeInFile == TRUE) {
 
	/* Read ParticleType into temporary buffer and Copy to ParticleType. */
 
	file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
        if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
        if( file_dsp_id == h5_error){my_exit(EXIT_FAILURE);}
 
	if (io_log) fprintf(log_fptr,"H5Dopen  with Name = particle_type\n");
 
	dset_id =  H5Dopen(file_id, "particle_type");
        if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
        if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) tempint);
        if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	h5_status = H5Sclose(file_dsp_id);
        if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	h5_status = H5Dclose(dset_id);
        if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
        if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleType[i] = tempint[i];
 
	for (i = 0; i < NumberOfParticles; i++)
	  if (ParticleType[i] < PARTICLE_TYPE_GAS ||
	      ParticleType[i] > NUM_PARTICLE_TYPES-1) {
	    fprintf(stderr, "file: %s: particle %"ISYM" has unknown type %"ISYM"\n",
		    name, i, ParticleType[i]);
	    ENZO_FAIL("");
	  }
 
      } else {
 
	/* Otherwise create the type. */
 
	for (i = 0; i < NumberOfParticles; i++)
	  ParticleType[i] = ReturnParticleType(i);
 
      } 
 
      /* Read ParticleAttributes. */
      if (AddParticleAttributes) {
	for (j = 0; j < NumberOfParticleAttributes; j++) {
	  ParticleAttribute[j] = new float[NumberOfParticles];
	  for (i=0; i < NumberOfParticles; i++)
	    ParticleAttribute[j][i] = 0;
	}
      } else {
	for (j = 0; j < NumberOfParticleAttributes; j++) {
 
	  file_dsp_id = H5Screate_simple((Eint32) 1, TempIntArray, NULL);
	  if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
	  if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	  if (io_log) fprintf(log_fptr,"H5Dopen with Name = %s\n",ParticleAttributeLabel[j]);
	  
	  dset_id =  H5Dopen(file_id, ParticleAttributeLabel[j]);
	  if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
	  if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
	  h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
	  if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
	  
	  h5_status = H5Sclose(file_dsp_id);
	  if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
	  
	  h5_status = H5Dclose(dset_id);
	  if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
	  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
	    
	  for (i = 0; i < NumberOfParticles; i++)
	    ParticleAttribute[j][i] = float(temp[i]);
 
	}
      } // ENDELSE AddParticleAttributes 

      delete [] temp;
      delete [] tempint;

#endif /* USE_HDF4 */
 
    } // end: if (MyProcessorNumber == ProcessorNumber)
  } // end: if (NumberOfParticles > 0 && ReadData)
 
  /* Close file. */
 
  if ( (MyProcessorNumber == ProcessorNumber) &&
       (NumberOfParticles > 0 || NumberOfBaryonFields > 0) 
       && ReadData ){

#ifdef USE_HDF4
    SDend(sd_id);
#else
    h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
#endif

  }
 
  if (MyProcessorNumber == ProcessorNumber)
    {
      if (io_log) fclose(log_fptr);
    }
 
  return SUCCESS;
 
}







#ifdef USE_HDF4
/* ----------------------------------------------------------------------                                          
   This routine reads one data field from the file using the appropriate                                           
   data model.  Note that it uses file pointers/handlers that are statically                                       
   declared above. */

int ReadField(float *temp, int Dims[], int Rank, char *name,
	      char *field_name)
{
  int dim;

  int32 sds_id, start[] = {0, 0, 0};
  int32 TempInt, TempIntArray[MAX_DIMENSION], attributes, num_type;
  char dummy[MAX_LINE_LENGTH];

  /* Find the next SDS which is not a coordinate variable. */

  sds_id = SDselect(sd_id, sds_index++);
  while (SDiscoordvar(sds_id)) {
    SDendaccess(sds_id);
    sds_id = SDselect(sd_id, sds_index++);
  }

  if (SDgetinfo(sds_id, dummy, &TempInt, TempIntArray, &num_type, &attributes) == HDF_FAIL) {
    fprintf(stderr, "error getting info from file %s (filed %s)\n", name, field_name);
    return FAIL;
  }
  
  /* check rank against this grid */
  
  if (TempInt != Rank) {
    fprintf(stderr, "HDF rank (%d) does not match GridRank.\n", TempInt);
    return FAIL;
  }

  /* check dimensions of HDF file against this grid */
  
  for (dim = 0; dim < Rank; dim++)
    if (TempIntArray[Rank-dim-1] != Dims[dim]) {
      fprintf(stderr, "HDF file dimensions do not match GridDimensions.\n");
      return FAIL;
    }

  if (SDreaddata(sds_id, start, (int32 *) NULL, TempIntArray, (void *) temp)
      == HDF_FAIL) {
    fprintf(stderr, "Error reading data from file %s (field %s).\n", name, field_name);
    return FAIL;
  }
  SDendaccess(sds_id);

  return SUCCESS;
}
#endif /* USE_HDF4 */

