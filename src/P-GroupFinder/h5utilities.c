/*********************************************************************
 *
 *  file    : h5utilities.cpp
 *
 *  Project : Visualization of Adaptive Mesh Refinement Data
 *
 *  Company : Zuse Institute Berlin
 *            All rights reserved.
 *
 *  Author  : Ralf Kaehler                           
 *
 *  Date    : 26.01.2006
 *
 *********************************************************************/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <string.h>

#include <hdf5.h>


// need the following ifdef to make standalone version of carpet2amira 
// independent of hxhdf5 and hence of amira kernel libs ...
//#define CARPET2AMIRA_MAKESTANDALONE
//#ifndef CARPET2AMIRA_MAKESTANDALONE
//#include <hxhdf5/hxhdf5.h>
//#endif

#include "h5utilities.h"


int checkErr(int arg, const char* name) {
  int callAssert = 0;
  if (arg<0) {
    fprintf(stderr,"Error accessing %s. \n",name); 
    assert(!callAssert);
    return 1;
  }
  return 0;
}



hid_t openH5File(const char* filename)
{

  hid_t fileId = -1;

  int err = 0;

  // need the following ifdef to make standalone version of carpet2amira 
  // independent of hxhdf5 and hence of amira kernel libs ...

  //#ifdef CARPET2AMIRA_MAKESTANDALONE
#define DOIT
#ifdef DOIT
      
  // try to open HDF5 files
  fileId = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);

#else

  hid_t prop = H5Pcreate (H5P_FILE_ACCESS);

  if (prop < 0) { err = 1; }

  if (!err) {
    if (hx5set_driver (filename, prop) < 0) { err = 2; }
  }

  if (!err) {
    fileId = H5Fopen (filename, H5F_ACC_RDONLY, prop);
    if (fileId<0) {
      err = 3;
    }
  }

  if (prop >= 0) {
    H5Pclose (prop);
  }

#endif
#undef DOIT
  
  return fileId;

}



int readAttribute( const hid_t loc_id, 
		   const hid_t type_id, 
		   const char* attrName, 
		   void* result)
{

  hid_t attr;
  int err = 0;

  err |= checkErr ( attr = H5Aopen_name(loc_id,attrName),attrName );
  err |= checkErr ( H5Aread(attr,type_id,result),attrName );
  err |= checkErr ( H5Aclose(attr),attrName );
  
  return err;

}


int writeScalarAttribute( hid_t       loc_id, 
			  hid_t       type_id, 
			  const char* attrName,
			  const void* buffer)
{

  hid_t attrs, attr;
  
  int err = 0;

  err |= checkErr (attrs = H5Screate (H5S_SCALAR), attrName );                                    
  err |= checkErr (attr  = H5Acreate (loc_id, attrName, type_id, attrs, H5P_DEFAULT), attrName ); 
  err |= checkErr (H5Awrite (attr,  type_id, buffer), attrName  );                                
  err |= checkErr (H5Sclose (attrs), attrName );                                                  
  err |= checkErr (H5Aclose (attr), attrName );                                                   

  return (err);

}


int writeArrayAttribute( hid_t         loc_id, 
			 hid_t         type_id, 
			 const hsize_t dims, 
			 const char*   attrName,
			 const void*   buffer  )
{

  hid_t dataspace_id, attribute_id;
  int   err = 0;

 

  err |= checkErr ( dataspace_id = H5Screate_simple(1, &dims, NULL), attrName );
  err |= checkErr ( attribute_id = H5Acreate(loc_id, attrName, type_id, dataspace_id, H5P_DEFAULT), attrName );
  err |= checkErr ( H5Awrite(attribute_id, type_id, buffer), attrName );
  err |= checkErr ( H5Aclose(attribute_id), attrName );
  err |= checkErr ( H5Sclose(dataspace_id), attrName );

  return err;

}

int writeArrayDataset( hid_t         loc_id, 
			 hid_t         type_id, 
             int ndims,
			 hsize_t *dims, 
			 const char*   attrName,
			 const void*   buffer  )
{

  hid_t dataspace_id, dataset_id;
  int   err = 0;

  err |= checkErr ( dataspace_id = H5Screate_simple(ndims, dims, NULL), attrName );
  err |= checkErr ( dataset_id = H5Dcreate(loc_id, attrName, type_id, dataspace_id, H5P_DEFAULT), attrName );
  err |= checkErr ( H5Dwrite(dataset_id, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer), attrName );
  err |= checkErr ( H5Dclose(dataset_id), attrName );
  err |= checkErr ( H5Sclose(dataspace_id), attrName );

  return err;

}


int writeStringAttribute(hid_t loc_id, const char* name, const char* data)
{

  int err = 0;

  hid_t	  dataspace, attr, datatype;
  int str_length = strlen(data);
  if (str_length == 0) str_length = 1;

  err |= checkErr( dataspace = H5Screate(H5S_SCALAR), data) ;
  err |= checkErr( datatype  = H5Tcopy(H5T_C_S1), data);
  err |= checkErr( H5Tset_size(datatype,str_length), data);
  err |= checkErr( attr = H5Acreate(loc_id, name, datatype, dataspace, H5P_DEFAULT), data);

  err |= checkErr( H5Awrite(attr, datatype, data), data); 
  err |= checkErr( H5Sclose(dataspace), data); 
  err |= checkErr( H5Aclose(attr), data); 

  return err;

}



#ifdef UNUSED
nativeTypeId h5Type2NativeType (const hid_t h5type) 
{

  assert(h5type>0);

  const H5T_class_t h5class = H5Tget_class(h5type);
  const size_t      h5prec  = H5Tget_precision (h5type);
  const H5T_sign_t  h5sign  = H5Tget_sign(h5type);

  if (h5class==H5T_FLOAT) {

    if (h5prec == 32) { return NTID_float; }
    else if (h5prec == 64) { return NTID_double; }
    else { assert(0); }

  }
  else if (h5class == H5T_INTEGER) {

    if (h5sign == H5T_SGN_2) {
      if (h5prec == 16) { return NTID_int16; }
      else if (h5prec == 32) { return NTID_int32; }
      else { assert(0); }
    }
    else if (h5sign == H5T_SGN_NONE) {

      if (h5prec == 8) { return NDTI_uint8; }
      else if (h5prec == 16) { return NTID_uint16; }
      else { assert(0); }

    }
    else {
      assert(0);
    }

  }
  else {
    assert(0);
  }

  assert(0);

  return NTID_unknown;

};
#endif /* UNUSED */

