/************************************************************************
 HDF INPUT FOR P-GROUPFINDER
 By John Wise

 Created:       02 Jul 2005
 Last Modified: 
 History:
************************************************************************/

#ifdef USE_HDF4

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "hdf.h"
#include "allvars.h"
#include "proto.h"

void ReadParticleField_FLOAT64 (int32 sd_id, char *label, int nPart, 
				double **data)
{

  char buf1[200];
  int n;

  // HDF variables
  int32 	 sds_id, sds_index;
  int32 	 DummyIntArray[3];
  int32 	 HDFdummyInt, num_type, attributes;
  intn 		 hdf_status;
  float32 	*temp  = NULL;
  int32 start[]        = {0, 0, 0};
  int32 TempIntArray[] = {0, 0, 0};

  if ((sds_index = SDnametoindex(sd_id, label)) == FAIL) {
    fprintf(stderr, "ReadParticleField: cannot read %s\n", label);
    fprintf(stderr, "SD_ID = %d, SDS_INDEX = %d\n", sd_id, sds_index);
    MPI_Finalize();
    exit(1);
  }
  
  sds_id = SDselect(sd_id, sds_index);
    
  // Get datatype (float, double, long_double)
  hdf_status = SDgetinfo(sds_id, buf1, &HDFdummyInt, DummyIntArray, 
			 &num_type, &attributes);
  switch (num_type) {
  case DFNT_FLOAT32:
    temp = (float32 *) malloc(sizeof(float32) * nPart);
    break;
  case DFNT_FLOAT64:
    temp = (float32 *) malloc(sizeof(float64) * nPart);
    break;
  case DFNT_FLOAT128:
    temp = (float32 *) malloc(sizeof(long double) * nPart);
    break;
  }
  
  TempIntArray[0] = (int32) nPart;
  hdf_status = SDreaddata(sds_id, start, (int32 *) NULL, TempIntArray, 
			  (VOIDP) temp);
  
  (*data) = (double *) malloc(sizeof(double) * nPart);
  float64 *temp64 = (float64 *) temp;
  long double *temp128 = (long double *) temp;
  
  switch (num_type) {
  case DFNT_FLOAT32:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (double) temp[n];
    break;
  case DFNT_FLOAT64:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (double) temp64[n];
    break;
  case DFNT_FLOAT128:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (double) temp128[n];
    break;
  }  

  free(temp);

}

/************************************************************************/

void ReadParticleField_FLOAT32 (int32 sd_id, char *label, int nPart, 
				float **data)
{

  char buf1[200];
  int n;

  // HDF variables
  int32 	 sds_id, sds_index;
  int32 	 DummyIntArray[3];
  int32 	 HDFdummyInt, num_type, attributes;
  intn 		 hdf_status;
  float32 	*temp  = NULL;
  int32 start[]        = {0, 0, 0};
  int32 TempIntArray[] = {0, 0, 0};

  if ((sds_index = SDnametoindex(sd_id, label)) == FAIL) {
    fprintf(stderr, "ReadParticleField: cannot read %s\n", label);
    fprintf(stderr, "SD_ID = %d, SDS_INDEX = %d\n", sd_id, sds_index);
    MPI_Finalize();
    exit(1);
  }
  
  sds_id = SDselect(sd_id, sds_index);
    
  // Get datatype (float, double, long_double)
  hdf_status = SDgetinfo(sds_id, buf1, &HDFdummyInt, DummyIntArray, 
			 &num_type, &attributes);
  switch (num_type) {
  case DFNT_FLOAT32:
    temp = (float32 *) malloc(sizeof(float32) * nPart);
    break;
  case DFNT_FLOAT64:
    temp = (float32 *) malloc(sizeof(float64) * nPart);
    break;
  case DFNT_FLOAT128:
    temp = (float32 *) malloc(sizeof(long double) * nPart);
    break;
  }
  
  TempIntArray[0] = (int32) nPart;
  hdf_status = SDreaddata(sds_id, start, (int32 *) NULL, TempIntArray, 
			  (VOIDP) temp);
  
  (*data) = (float *) malloc(sizeof(float) * nPart);
  float64 *temp64 = (float64 *) temp;
  long double *temp128 = (long double *) temp;
  
  switch (num_type) {
  case DFNT_FLOAT32:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (float) temp[n];
    break;
  case DFNT_FLOAT64:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (float) temp64[n];
    break;
  case DFNT_FLOAT128:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (float) temp128[n];
    break;
  }  

  free(temp);

}

/************************************************************************/

void ReadParticleField_INT (int32 sd_id, char *label, int nPart, int **data)
{

  char buf1[200];
  int n;

  // HDF variables
  int32 	 sds_id, sds_index;
  int32 	 DummyIntArray[3];
  int32 	 HDFdummyInt, num_type, attributes;
  intn 		 hdf_status;
  int32 	*temp  = NULL;
  int32 start[]        = {0, 0, 0};
  int32 TempIntArray[] = {0, 0, 0};

  if ((sds_index = SDnametoindex(sd_id, label)) == FAIL) {
    fprintf(stderr, "ReadParticleField: cannot read %s\n", label);
    fprintf(stderr, "SD_ID = %d, SDS_INDEX = %d\n", sd_id, sds_index);
    MPI_Finalize();
    exit(1);
  }
  
  sds_id = SDselect(sd_id, sds_index);
    
  // Get datatype (float, double, long_double)
  hdf_status = SDgetinfo(sds_id, buf1, &HDFdummyInt, DummyIntArray, 
			 &num_type, &attributes);
  switch (num_type) {
  case DFNT_INT32:
    temp = (int32 *) malloc(sizeof(int32) * nPart);
    break;
  case DFNT_INT64:
    fprintf(stderr, "No support for 64-bit integers yet.\n");
    MPI_Finalize();
    exit(1);
    break;
  }
  
  TempIntArray[0] = (int32) nPart;
  hdf_status = SDreaddata(sds_id, start, (int32 *) NULL, TempIntArray, 
			  (VOIDP) temp);
  
  (*data) = (int *) malloc(sizeof(int) * nPart);
  //  int64 *temp64 = (int64 *) temp;
  
  switch (num_type) {
  case DFNT_INT32:
    for (n = 0; n < nPart; n++)
      (*data)[n] = (int) temp[n];
    break;
  case DFNT_INT64:
    break;
  }  

  free(temp);

}

/************************************************************************/

#endif
