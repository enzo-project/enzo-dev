#include <assert.h>
#include "hdf5.h"
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"


#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

//returns the wall time.
extern FILE * wall_ptr;
void wall_time_start(){
  char filename[20];
  sprintf(filename,"data_time_%04d",MyProcessorNumber);
  wall_ptr = fopen(filename,"w");
}
void wall_time_flush(){
  fclose(wall_ptr);
  char filename[20];
  sprintf(filename,"data_time_%04d",MyProcessorNumber);
  wall_ptr = fopen(filename,"a");

}
void wall_time_stop(){
  fclose(wall_ptr);
}  
void wall_time (char * string)
{
  //wall_ptr == null unless wall_time_start is called
  if( wall_ptr ){
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv,&tz);
    fprintf(wall_ptr,"%d TIME %s %f \n", MyProcessorNumber, string, tv.tv_sec + 1e-6*tv.tv_usec);
    fflush(wall_ptr);
  }
}

//This routine checks the list SingleGridDumpList for the given flag value.
int CheckForSingleGridDump(int flag){
  
  int WriteInThisLocal = FALSE;
  for( int SnapperJoe = 0; SnapperJoe < N_DbgWrites; SnapperJoe++){
    if( SingleGridDumpList[SnapperJoe] == flag ) {
      WriteInThisLocal = TRUE;
    }
  }
  
  return WriteInThisLocal;
}

//Writes an HDF5 cube.  Quick and dirty.  
// WriteCube(array, [nx,ny,nz], "ID string", dNum, gNum)
// prints out ID string to std err, along with file name.
// Filename is data111(dNum).grid(gNum)
// datasets are all called "monkey"
// Flow:
// 1.) create filename, print message
// 2.) define size of float 
// 3.) create file
// 3.5) invert array dimensions
// 4.) create dataspace, set
// 5.) close file,space,set.
void WriteSingleCube(float * array, int Dims[], char* string, int dNum, int gNum){
  
  hid_t       file_id, dataset_id, dataspace_id, float_type_id;
  herr_t      status, h5_status, h5_error = -1;
  int FieldRankOut = 3;
  hsize_t     DimsInv[FieldRankOut];
  
  char filename[20];
  
  sprintf(filename, "data111%4.4d.grid%4.4d",dNum,gNum);
  fprintf(stderr,"GPFS WriteCube: %s %s [%d,%d,%d]\n", string, filename, Dims[0],Dims[1],Dims[2]);
  
#define floatdcc double  
  int jj = sizeof(floatdcc);
  switch(jj)
    {
    case 0:
      float_type_id = H5T_NATIVE_INT;
      break;
    case 4:
      float_type_id = HDF5_R4;
      break;
    case 8:
      float_type_id = HDF5_R8;
      break;
    case 16:
      float_type_id = HDF5_R16;
      break;
    default:
      printf("INCORRECT FLOAT DEFINITION\n");
    }
  
  
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  for(int moo=0;moo<FieldRankOut;moo++)
    DimsInv[moo]=Dims[FieldRankOut-moo-1];
  
  
  //Create Dataspace 
  dataspace_id=H5Screate_simple(FieldRankOut, DimsInv, NULL);
  
  //create set
  //                       duh, name,      datatype,  shape of data, Something I dont get
  dataset_id = H5Dcreate(file_id, "monkey", float_type_id, dataspace_id, H5P_DEFAULT);
  
  //Write the Data Set
  //                (set, memory type, mem. space, file space, transfer shit, actual data)
  fprintf(stderr,"Writing set\n");
  status = H5Dwrite(dataset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		    array);
  
  
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);
  
  
}

void commas(char *str, int numin){

  char temp[100]="";
  int num = numin;
  sprintf(str, "");
  for(int i=10; i>=0; i--){
    int dig = (int) pow((double) 10,(double) i);
    if((int) numin/dig < 1 ) continue;
    sprintf(temp, "%i", num/dig);
    strcat(str, temp);
    if( i % 3 == 0 && i != 0 ) strcat(str, ",");
    num=num-(num/dig)*dig;
  }

}
