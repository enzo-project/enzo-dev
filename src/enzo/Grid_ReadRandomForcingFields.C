/***********************************************************************
/
/  GRID CLASS (READ GRID)
/
/  written by: Alexei Kritsuk
/  date:       January, 2004
/  modified  : Robert Harkness
/  date:       May, 2008
/  modified2:  Michael Kuhlen, October 2010, HDF5 hierarchy
/
/  PURPOSE: At restart reads initial velocity fields from the first output
/           data file(s) (e.g. data_0000.grid000?) and stores them as
/           RandomForcingField[].
/
************************************************************************/
 
//  Input a grid from file pointer fpt
 
#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
 


 
//#include "performance.h"
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
 
int grid::ReadRandomForcingFields(FILE *fptr, char DataFilename[])
{
 
  int i, j, k, dim, field, size, active_size;
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
 
  int ActiveDim[MAX_DIMENSION];
  
  fpos_t position; 
  FILE *log_fptr;
 
  hid_t       file_id, dset_id;
  hid_t       float_type_id, FLOAT_type_id;
  hid_t       file_dsp_id;
  hid_t       num_type;
 
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     TempIntArray[MAX_DIMENSION];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int         num_size;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
 
  int ii = sizeof(float32);
 
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
      break;
 
    case 8:
      FLOAT_type_id = HDF5_R8;
      break;
 
    case 16:
      FLOAT_type_id = HDF5_R16;
      break;
 
    default:
      printf("INCORRECT FLOAT DEFINITION\n");
 
  }
 
 
  /* Assume general grid class data are known and
     read velocity fields only (as RandomForcingFields).
     Do not modify the existing current BaryonFields and grid class data. */

  int MinNumberOfFields = GridRank + 1;
  if( EquationOfState == 0 ) MinNumberOfFields++;
  if (NumberOfBaryonFields < MinNumberOfFields) {
    fprintf(stderr, "Error: No.  Baryon Fields (%"ISYM") => Nothing to Force. Right?.\n",
            NumberOfBaryonFields);
    ERROR_MESSAGE;
  }

  if (HierarchyFileInputFormat == 1) {
    /* Store the file position indicator. */
    if (fgetpos(fptr, &position) != 0)
      WARNING_MESSAGE;
    
    /* Read the filename where the current baryon fields are; assume that
       initial fields were sitting in a file with the same name but different
       number; change the current number into '0000' and, thus, prepare the
       required name. */
    
    if (fsetpos(fptr, &BaryonFileNamePosition) != 0)
      ERROR_MESSAGE;
    if (fscanf(fptr, "BaryonFileName = %s\n", name) != 1) {
      fprintf(stderr, "Error reading BaryonFileName.\n");
      ERROR_MESSAGE;
    }
  }

  if (HierarchyFileInputFormat % 2 == 0) {
    strcpy(name, DataFilename);
  }

  /* this wont work if file prefix contains anoter ".grid" and/or
     if restart is invoked for cycle # > 9999 */
  
  char * numberEnd = NULL;
  if ( (numberEnd = strstr(name, ".grid")) == NULL )
    ERROR_MESSAGE;
  *(numberEnd - 1) = '0';
  *(numberEnd - 2) = '0';
  *(numberEnd - 3) = '0';
  *(numberEnd - 4) = '0';
 
  /* check the name. */
 
  if (debug)
    printf("ReadRandomForcingFields from: %s\n", name);
  if ( strstr(name, "0000.grid") == NULL )
    ERROR_MESSAGE;
 
  /* read fields */
 
  if (MyProcessorNumber == ProcessorNumber) {
    char *logname = new char[MAX_NAME_LENGTH];
    strcpy(logname, name);
    strcat(logname, ".hdf.log");
    if (io_log) log_fptr = fopen(logname, "a");
    if (io_log) fprintf(log_fptr,"H5Fopen with Name %s\n",name);
    file_id = H5Fopen(name,  H5F_ACC_RDONLY, H5P_DEFAULT);
    if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
    if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
  }
 
  if (MyProcessorNumber == ProcessorNumber) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                         Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("GRRFF: Error in IdentifyPhysicalQuantities.\n");
    }
    int vel = Vel1Num;
    printf("RandomForcing: Fields %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" \n",
           DensNum, TENum, Vel1Num, Vel2Num, Vel3Num);
 
 
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
 
    float32 *temp = new float32[active_size];
 
    /* skip fields written before the first velocity field;
       dirty trick, sorry; this is done only once at a restart. */
 
    for (i = 0; i < vel; i++) {
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", DataLabel[i]);
 
      dset_id =  H5Dopen(file_id, DataLabel[i]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
      //      LCAPERF_COUNT_READ(dset_id, float_type_id, H5S_ALL);
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
    }
 
    /* loop over fields, reading each one */
 
    for (dim = 0; dim < GridRank; dim++) {
 
      /* get data into temporary array */
 
      file_dsp_id = H5Screate_simple((Eint32) GridRank, OutDims, NULL);
      if (io_log) fprintf(log_fptr, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
      if( file_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", DataLabel[vel+dim]);
 
      dset_id =  H5Dopen(file_id, DataLabel[vel+dim]);
      if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
      if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dread(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
      if (io_log) fprintf(log_fptr, "H5Dread: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
      //JBPERF_COUNT_READ(dset_id, float_type_id, H5S_ALL);
 
      h5_status = H5Sclose(file_dsp_id);
      if (io_log) fprintf(log_fptr, "H5Sclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      h5_status = H5Dclose(dset_id);
      if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
      if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
      /* copy velocity from active region into whole grid forcing field;
         put zeroes into ghost zones. */
 
      if (RandomForcingField[dim] == NULL)
        RandomForcingField[dim] = new float[size];
      for (i = 0; i < size; i++)
        RandomForcingField[dim][i] = 0;
 
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    RandomForcingField[dim][i + j*GridDimension[0] +
				      k*GridDimension[0]*GridDimension[1]] =
	      float(temp[(i-GridStartIndex[0])                         +
	                 (j-GridStartIndex[1])*ActiveDim[0]            +
	                 (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ]);
 
    } // end: loop over fields
 
    delete [] temp;
 
  }  // end: if (MyProcessorNumber == ProcessorNumber)
 
  /* Close file. */
 
  if ( (MyProcessorNumber == ProcessorNumber) &&
       (NumberOfBaryonFields > 0) )
  {
     h5_status = H5Fclose(file_id);
       if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
       if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
  }
 
  if (MyProcessorNumber == ProcessorNumber)
  {
    if (io_log) fclose(log_fptr);
  }
 
  if (HierarchyFileInputFormat == 1) {
    /* Set the file position indicator. */
    if (fsetpos(fptr, &position) != 0)
      WARNING_MESSAGE;
  }

  return SUCCESS;
 
}
