/*********************************************************************
 *
 *  file    : h5utilities.h
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

#ifndef  _HDF5_UTILITIES_
#define  _HDF5_UTILITIES_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>

typedef int nativeTypeId;
const nativeTypeId NTID_unknown = -1, NTID_float = 0, NTID_double = 1, 
  NTID_int16 = 2, NTID_int32 = 3, NDTI_uint8 = 4, NTID_uint16 = 5;

////////////////////////////////////////////////////
///// HDF5 wrapper functionality
////////////////////////////////////////////////////

int checkErr( int err, const char* name, const bool callAssert=false );

   
int readAttribute( const hid_t loc_id, 
		   const hid_t type_id, 
		   const char* attrName, 
		   void*       result,
		   const bool  callAssert=false);
   

int writeScalarAttribute( hid_t       loc_id, 
			  hid_t       type_id, 
			  const char* attrName,
			  const void* buffer  );


int writeArrayAttribute( hid_t         loc_id, 
			 hid_t         type_id, 
			 const hsize_t dims, 
			 const char*   attrName,
			 const void*   buffer  );

int writeArrayDataset( hid_t         loc_id, 
			 hid_t         type_id, 
             int ndims,
			 hsize_t *dims, 
			 const char*   attrName,
			 const void*   buffer  );


int writeStringAttribute( hid_t loc_id, 
                          const char* datasetname, 
                          const char* data);

hid_t openH5File(const char* filename);

nativeTypeId h5Type2NativeType (const hid_t h5type); 

#endif
