/************************************************************************
 HDF5 INPUT FOR P-GROUPFINDER
 By John Wise

 Created:       02 Jul 2005
 Last Modified: 03 Mar 2009 (modified HDF4 version to read Packed AMR)
 History:
************************************************************************/

#ifdef USE_HDF5
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hdf5.h"
#include "allvars.h"
#include "proto.h"

void ReadParticleFieldHDF5_DOUBLE(hid_t group_id, char *label, int nPart, double **data)
{

  int n;
  // HDF variables
  hid_t data_id, num_type, data_type_id;
  hsize_t freal_size;
  float 	*temp  = NULL;

  if ((data_id = H5Dopen(group_id, label)) < 0) {
    fprintf(stderr, "ReadParticleField: cannot read %s\n", label);
    fprintf(stderr, "GROUP_ID = %d, DATA_ID = %d\n", group_id, data_id);
    MPI_Finalize();
    exit(1);
  }
  
  // Get datatype (float, double, long_double)
  num_type = H5Dget_type(data_id);
  freal_size = H5Tget_size(num_type);
  switch (freal_size) {
  case 4:
    data_type_id = H5T_NATIVE_FLOAT;
    temp = (float *) malloc(sizeof(float) * nPart);
    break;
  case 8:
    data_type_id = H5T_NATIVE_DOUBLE;
    temp = (float *) malloc(2*sizeof(float) * nPart);
    break;
  case 12:
  case 16:
    data_type_id = H5T_NATIVE_LDOUBLE;
    temp = (float *) malloc(4*sizeof(float) * nPart);
    break;
  }

  if (H5Dread(data_id, data_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	      (void*) temp) < 0) {
    fprintf(stderr, "Error reading %s\n", label);
    MPI_Finalize();
    exit(1);
  }
      
  (*data) = (double *) malloc(sizeof(double) * nPart);
  double *temp64 = (double *) temp;
  long double *temp128 = (long double *) temp;
  
  switch (freal_size) {
  case 4:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (double) temp[n];
    break;
  case 8:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (double) temp64[n];
    break;
  case 12:
  case 16:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (double) temp128[n];
    break;
  }  

  free(temp);
  H5Dclose(data_id);

}

/************************************************************************/

void ReadParticleFieldHDF5_FLOAT(hid_t group_id, char *label, int nPart, float **data)
{

  int n;
  // HDF variables
  hid_t data_id, num_type, data_type_id;
  hsize_t freal_size;
  float 	*temp  = NULL;

  if ((data_id = H5Dopen(group_id, label)) < 0) {
    fprintf(stderr, "ReadParticleField: cannot read %s\n", label);
    fprintf(stderr, "GROUP_ID = %d, DATA_ID = %d\n", group_id, data_id);
    MPI_Finalize();
    exit(1);
  }
  
  // Get datatype (float, double, long_double)
  num_type = H5Dget_type(data_id);
  freal_size = H5Tget_size(num_type);
  switch (freal_size) {
  case 4:
    data_type_id = H5T_NATIVE_FLOAT;
    temp = (float *) malloc(sizeof(float) * nPart);
    break;
  case 8:
    data_type_id = H5T_NATIVE_DOUBLE;
    temp = (float *) malloc(2*sizeof(float) * nPart);
    break;
  case 12:
  case 16:
    data_type_id = H5T_NATIVE_LDOUBLE;
    temp = (float *) malloc(4*sizeof(float) * nPart);
    break;
  }

  if (H5Dread(data_id, data_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	      (void*) temp) < 0) {
    fprintf(stderr, "Error reading %s\n", label);
    MPI_Finalize();
    exit(1);
  }
      
  (*data) = (float *) malloc(sizeof(float) * nPart);
  double *temp64 = (double *) temp;
  long double *temp128 = (long double *) temp;
  
  switch (freal_size) {
  case 4:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (float) temp[n];
    break;
  case 8:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (float) temp64[n];
    break;
  case 12:
  case 16:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (float) temp128[n];
    break;
  }  

  free(temp);
  H5Dclose(data_id);

}

/************************************************************************/

void ReadParticleFieldHDF5_INT(hid_t group_id, char *label, int nPart, PINT **data)
{

  int n;
  // HDF variables
  hid_t data_id, data_type_id;
  PINT 	*temp  = NULL;

  if ((data_id = H5Dopen(group_id, label)) < 0) {
    fprintf(stderr, "ReadParticleField: cannot read %s\n", label);
    fprintf(stderr, "GROUP_ID = %d, DATA_ID = %d\n", group_id, data_id);
    MPI_Finalize();
    exit(1);
  }
  
  data_type_id = HDF5_PINT;
  temp = (PINT *) malloc(sizeof(PINT) * nPart);

  if (H5Dread(data_id, data_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	      (void*) temp) < 0) {
    fprintf(stderr, "Error reading %s\n", label);
    MPI_Finalize();
    exit(1);
  }
      
  (*data) = (PINT *) malloc(sizeof(PINT) * nPart);
  for (n = 0; n < nPart; n++)
    (*data)[n] = (PINT) temp[n];

  free(temp);
  H5Dclose(data_id);

}
#endif
