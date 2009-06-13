#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"

#include "AnalysisBaseClass.h"

void my_exit(int exit_status);

void AnalysisBaseClass::HDF5CreateFile( char *name ){

  hid_t       file_id;
  herr_t      h5_status;
  herr_t      h5_error = -1;

  //Don't create existing files.
  FILE * test = fopen(name,"rb");
  if( test ){
    fclose( test );
    return;
  }

  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if( file_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CreateFile, cannot create file %s\n", name);
    my_exit(EXIT_FAILURE);
  }

  h5_status = H5Fclose(file_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CreateFile, cannot close file %s\n", name);
    my_exit(EXIT_FAILURE);
  }

}

hid_t AnalysisBaseClass::HDF5OpenFile( char *name ){

  hid_t       file_id;
  herr_t      h5_error = -1;

  file_id = H5Fopen(name, H5F_ACC_RDWR, H5P_DEFAULT);

  if( file_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5OpenFile, cannot open file %s\n", name);
    my_exit(EXIT_FAILURE);
  }

  return file_id;
}

void AnalysisBaseClass::HDF5CloseFile( hid_t file_id ){

  herr_t      h5_status;
  herr_t      h5_error = -1;

  h5_status = H5Fclose(file_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CloseFile, cannot close file\n" );
    my_exit(EXIT_FAILURE);
  }

}

hid_t AnalysisBaseClass::HDF5CreateGroup( hid_t loc_id, char *name ){

  hid_t group_id;
  herr_t h5_status;
  herr_t      h5_error = -1;

  group_id = H5Gcreate(loc_id, name, 0);

  if( group_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CreateGroup, cannot create group %s\n", name);
    my_exit(EXIT_FAILURE);
  }
  if( debug )
    fprintf(stderr,"HDF5CreateGroup: created %s\n",name);
  return group_id;
}


hid_t AnalysisBaseClass::HDF5OpenGroup( hid_t loc_id, char *name ){

  hid_t group_id;
  herr_t h5_status;
  herr_t      h5_error = -1;

  group_id = H5Gopen(loc_id, name);

  if( group_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5GetGroup, cannot open group %s\n", name);
    my_exit(EXIT_FAILURE);
  }
  return group_id;
}

void AnalysisBaseClass::HDF5CloseGroup( hid_t group_id ){

  herr_t      h5_status;
  herr_t      h5_error = -1;

  h5_status = H5Gclose(group_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CloseGroup, cannot close group\n" );
    my_exit(EXIT_FAILURE);
  }
}

void AnalysisBaseClass::HDF5ReadDataset( hid_t group_id,
				     char *dset_name,
				     float *data ){

  herr_t h5_status;
  herr_t      h5_error = -1;

  hid_t dset_id =  H5Dopen(group_id, dset_name);

  if( dset_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5ReadDataset, cannot open read dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  h5_status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5ReadDataset, cannot read dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }
  if( debug )
    fprintf(stderr, "AnalysisBaseClass::HDF5ReadDataset, successfully wrote %s.\n",dset_name);
  h5_status = H5Dclose(dset_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5ReadDataset, cannot close dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }
}
 
void AnalysisBaseClass::HDF5MakeDataset( hid_t group_id, 
					 char *dset_name,
					 int rank, hsize_t dims[], 
					 float *data,
					 char *units ){
  herr_t      h5_status;
  herr_t      h5_error = -1;

  // WARN WARN WARN
  // harcoding 32-bit IO
  hid_t float_type_id = HDF5_R4;
  hid_t file_type_id = HDF5_FILE_R4;

  int i, size = 1;
  for(i = 0; i < rank; i++)
    size *= int(dims[i]);

  Eflt32 *buffer = new Eflt32[size];

  for(i = 0; i < size; i++)  
    buffer[i] = Eflt32(data[i]);

  hid_t file_dsp_id = H5Screate_simple(rank, dims, NULL);

  if( file_dsp_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot open dataspace\n");
    my_exit(EXIT_FAILURE);
  }

  hid_t dset_id =  H5Dcreate(group_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);

  if( dset_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot create dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }
  
  /* Write the dataset. */
  h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot write dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  if(units)
    HDF5WriteStringAttr(dset_id, "Units", units);
  
  h5_status = H5Sclose(file_dsp_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot close dataspace\n");
    my_exit(EXIT_FAILURE);
  }

  h5_status = H5Dclose(dset_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot close dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  delete [] buffer;
}

void AnalysisBaseClass::HDF5MakeDataset( hid_t group_id, 
					 char *dset_name,
					 int rank, hsize_t dims[], 
					 int *data,
					 char *units ){
  herr_t      h5_status;
  herr_t      h5_error = -1;

  // WARN WARN WARN
  // harcoding 32-bit IO
  hid_t float_type_id = HDF5_I4;
  hid_t file_type_id = HDF5_FILE_I4;

  int i, size = 1;
  for(i = 0; i < rank; i++)
    size *= int(dims[i]);

  Eint32 *buffer = new Eint32[size];

  for(i = 0; i < size; i++)  
    buffer[i] = Eint32(data[i]);

  hid_t file_dsp_id = H5Screate_simple(rank, dims, NULL);

  if( file_dsp_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot open dataspace\n");
    my_exit(EXIT_FAILURE);
  }

  hid_t dset_id =  H5Dcreate(group_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);

  if( dset_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot create dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }
  
  /* Write the dataset. */
  h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot write dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  if(units)
    HDF5WriteStringAttr(dset_id, "Units", units);
  
  h5_status = H5Sclose(file_dsp_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot close dataspace\n");
    my_exit(EXIT_FAILURE);
  }

  h5_status = H5Dclose(dset_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot close dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  delete [] buffer;
}

void AnalysisBaseClass::HDF5WriteStringAttr(hid_t dset_id, char *Alabel, char *String){

  hid_t       attr_id;
  hid_t       attr_dsp_id;
  hid_t       attr_type_id;
  herr_t      h5_status;
  herr_t      h5_error = -1;

  const char  *NoString = "none";

  attr_dsp_id = H5Screate(H5S_SCALAR);

  if( attr_dsp_id == h5_error ){ 
    fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to open dataspace\n"); 
    my_exit(EXIT_FAILURE);
  }

  attr_type_id = H5Tcopy(H5T_C_S1);
  H5Tset_size(attr_type_id, 80);

  attr_id = H5Acreate(dset_id, Alabel, attr_type_id,  attr_dsp_id, H5P_DEFAULT);

  if(  attr_id == h5_error ){ 
    fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to create attribute %s\n", Alabel ); 
    my_exit(EXIT_FAILURE);
  }

  if( strlen(String) > 0 ){
    h5_status = H5Awrite(attr_id, attr_type_id, (void *) String);
    if(  h5_status == h5_error ){ 
      fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to write string %s\n", String); 
      my_exit(EXIT_FAILURE); 
    }
  }else{
    h5_status = H5Awrite(attr_id, attr_type_id, (void *) NoString);
    if(  h5_status == h5_error ){ 
      fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to write string %s\n", NoString); 
      my_exit(EXIT_FAILURE); 
    }
  }

  h5_status = H5Aclose(attr_id);

  if(  h5_status == h5_error ){ 
    fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to close attribute\n"); 
    my_exit(EXIT_FAILURE); 
  }

  h5_status = H5Sclose(attr_dsp_id);

  if(  h5_status == h5_error ){ 
    fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to close dataspace\n"); 
    my_exit(EXIT_FAILURE); 
  }
}
