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

// int ifdef screws around with hdf5.  Undefine it.
#ifdef LARGE_INTS
#undef int
#endif

#include <hdf5.h>

// Redefine them
#ifdef LARGE_INTS
#define int long_int
#endif

#define _CHECK_CONSISTENCY_

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
  void AMRHDF5CloseSeparateParticles();

  void AMRHDF5CreateSeparateParticles( const char*      fileName, 
				       const int        ParticlesOn,
				       const int        nParticleAttr,
				       bool&            error) ;
  
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

  herr_t writeParticles2( const int nPart,
			  const int nAttributes,
			  const int nBaryonFields,
			  const int Rank,
			  void **pos,
			  void **vel,
			  void *type,
			  void *ID,
			  void *mass,
			  void **attr,
			  int& alreadyopenedentry,
			  int& NumberOfStarParticlesOnProcOnLvlEntry,
			  
			  const int    timeStep,
			  const double physicalTime,
			  const double redshift,
			  const int    levelIndex,
			  const double *delta,
			  
			  const double *physicalOrigin,
			  const Eint64    *integerOrigin,
			  const int    *bboxflags,
			  const int    *nghostzones ) ;
  
  herr_t writeSeparateParticles ( const int nPart,
				  const int nAttributes,
				  const int Rank,
				  void **pos,
				  void **vel,
				  void *type,
				  void *ID,
				  void *mass,
				  void **attr,
				  double doubleTime,
				  double doubleRedshift,
				  int& alreadyopenedentry,
				  int& NumberOfStarParticlesOnProcEntry ) ;

  void IncreaseGridCount();
  void IncreaseParticleGridCount();
  void IncreaseOutputParticleCount();

 protected:

  hid_t      h5DataType;
  staggering staggerType;
  fieldtype  fieldType;

  hid_t  fileId, fileId_particle;
  FILE  *index_file;
  int    relRef[3];
  int    gridId, particlegridId, output_particle;

  double rootDelta[3];
  
};

#endif


