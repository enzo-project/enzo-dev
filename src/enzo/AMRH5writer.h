/*********************************************************************
 *
 *  file    : AMRH5writer.h
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

#ifndef  __AMRHDF5WRITER__
#define  __AMRHDF5WRITER__

// min/max macros conflict with stl
#undef min
#undef max

#include <hdf5.h>

#define _CHECK_CONSISTENCY_

typedef enum staggering  { VERTEX_CENTERED=0, CELL_CENTERED=1} staggering;
typedef enum fieldtype   { SCALAR=1, /*VECTOR=3*/} fieldtype;


class AMRHDF5Writer
{

 public:

  AMRHDF5Writer ();

  void AMRHDF5Create(/* filename of the file that will store the data */
		const char*      fileName, 

		/* relative refinement factor in x,y,z between two levels of refinement;
		   up to now only a factor of 2,2,2 is supported by all vis routines */
		const int*    relativeRefinementFactor,
		
		/* data type of the grid function, as classified by HDF5; 
		   compare: http://hdf.ncsa.uiuc.edu/HDF5/doc/PredefDTypes.html */
		hid_t   dataType,
		
		/* specifies the staggering of the data; currently only cell and 
		   vertex-centered data are supported */
		staggering stag,
		       
		/* type of field; currently only scalar data is supported */
		fieldtype  field_type,
		
		const int        cycle,
		const double     time,
		const double     redshift,
		const float      root_dx,
		const int        overwrite,
		const int        writeTimeFile,
		const int        nFields,
		const int        ParticlesOn,
		const int        nParticleAttr,
		char           **FieldNames,

		/* return true if an error occured */
		bool& error) ;
  
  ~AMRHDF5Writer() ;
  void AMRHDF5Close();

  herr_t WriteTextures(  const int    timeStep,
			 const double physicalTime,
			 const int    levelIndex,
			 const double *delta,
			 
			 const double *physicalOrigin,
			 const double *gridCenter,
			 const Eint64    *integerOrigin,
			 
			 const int    *dims,
			 const int    dim,
			 const int    nFields,
			 char   **names,
			 void   **dataPtr);
  
  herr_t WriteFlat(  const int    timeStep,
		     const double physicalTime,
		     const double redshift,
		     const int    levelIndex,
		     const double *delta,
		     
		     const double *physicalOrigin,
		     const Eint64    *integerOrigin,
		     const int    *bboxflags,
		     const int    *nghostzones,
		     
		     const int    *dims,
		     const int    fieldNum,
		     const int    nFields,
		     const int    nParticles,
		     const char   *name,
		     const void   *dataPtr);

  herr_t writeParticles( const int nPart,
			 const int nAttributes,
			 const int nBaryonFields,
			 const int Rank,
			 void **pos,
			 void **vel,
			 void *type,
			 void *ID,
			 void *mass,
			 void **attr );

  void IncreaseGridCount();

 protected:

  hid_t      h5DataType;
  staggering staggerType;
  fieldtype  fieldType;

  hid_t  fileId;
  FILE  *index_file;
  int    relRef[3];
  int    gridId;

  double rootDelta[3];
  
};

// Redefine them
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

#endif


