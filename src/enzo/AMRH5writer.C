/***********************************************************************
 *
 *  file    : AMRH5writer.cpp
 *
 *  Project : Visualization of Adaptive Mesh Refinement Data
 *
 *  Company : Zuse Institute Berlin
 *            All rights reserved.
 *
 *  Author  : Ralf Kaehler                           
 *
 *  Date    : 14.01.2006
 *
 *********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "h5utilities.h"

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "AMRH5writer.h"
#include "StarParticleData.h"

AMRHDF5Writer::AMRHDF5Writer() :
  h5DataType(H5T_NATIVE_FLOAT),
  staggerType(CELL_CENTERED),
  fieldType(SCALAR)
{
  
};

void AMRHDF5Writer::AMRHDF5Create( const char*      fileName, 
				   const int*       relativeRefinement,
				   const hid_t      dataType,
				   const staggering stag,
				   const fieldtype  field_type,
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
				   bool&            error)
{
  error=false;

  const char *ParticlePositionLabel[] = 
    {"particle_position_x", "particle_position_y", "particle_position_z"};
  const char *ParticleVelocityLabel[] = 
    {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
  const char *ParticleOtherLabel[] =
    {"particle_type", "particle_index", "particle_mass"};
#ifdef WINDS
  const char *ParticleAttributeLabel[] =
    {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x", 
     "particle_jet_y", "particle_jet_z", "typeia_fraction"};
#else
  const char *ParticleAttributeLabel[] = 
    {"creation_time", "dynamical_time", "metallicity_fraction", "typeia_fraction"};
#endif

  int i;
    
  for (i=0; i<3; i++) { 
    rootDelta[i] = -1.0;
    relRef[i] = relativeRefinement[i];
    if (relRef[i]!=2) {
      fprintf(stderr,"AMRHDF5Writer: Only relative refinement factors of 2,2,2 are supported yet. \n");
      error=true;
    }
  }  

  fileId = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if (fileId<0) {
    fprintf(stderr,"Can't open file %s for writing.\n", fileName);
    error = true;
  }

  /* Write time data to a text file */

  /* Create global attributes */

  hid_t groupId;
  int err = 0;

  gridId = 0;
  particlegridId = 0;

  /* Create a string array that has the field names */

  char buf[80];
  int NumberOfParticleFields = (ParticlesOn == TRUE) ? 9+nParticleAttr : 0;
  int NumberOfAllFields = nFields + NumberOfParticleFields;
  int iField = 0;
  char **HDF5_FieldNames = new char*[NumberOfAllFields];
  for (i = 0; i < NumberOfAllFields; i++)
    HDF5_FieldNames[i] = new char[64];

  for (i = 0; i < nFields; i++)
    strcpy(HDF5_FieldNames[iField++], FieldNames[i]);

  if (ParticlesOn == TRUE) {
    for (i = 0; i < 3; i++)
      strcpy(HDF5_FieldNames[iField++], ParticlePositionLabel[i]);
    for (i = 0; i < 3; i++)
      strcpy(HDF5_FieldNames[iField++], ParticleVelocityLabel[i]);
    for (i = 0; i < 3; i++)
      strcpy(HDF5_FieldNames[iField++], ParticleOtherLabel[i]);
    for (i = 0; i < nParticleAttr; i++)
      strcpy(HDF5_FieldNames[iField++], ParticleAttributeLabel[i]);
  } // ENDIF ParticlesOn

  /* Write global attributes */

  char str[100];
  hid_t atype;
  strcpy(str, "/Parameters and Global Attributes");
  err |= checkErr(groupId = H5Gcreate(fileId, str, 0), str);
  err |= writeScalarAttribute(groupId, H5T_NATIVE_INT, "staggering", &stag);
  err |= writeScalarAttribute(groupId, H5T_NATIVE_INT, "NumberOfFields", &nFields);
  err |= writeScalarAttribute(groupId, H5T_NATIVE_INT, "NumberOfParticleFields", 
			      &NumberOfParticleFields);

  atype = H5Tcopy(H5T_C_S1);
  err |= H5Tset_size(atype, H5T_VARIABLE);
  err |= writeArrayAttribute(groupId, atype, NumberOfAllFields, "FieldNames", 
			     HDF5_FieldNames);

//  err |= writeScalarAttribute(groupId, H5T_NATIVE_INT, "ParticlesPresent", 
//			      &ParticlesOn);

  H5Gclose(groupId);

  H5Fflush(fileId, H5F_SCOPE_LOCAL);

  /* Open the accompanying binary index file */

  strcpy(str, fileName);
  strcat(str, ".idx");
  if ((index_file = fopen(str, "w")) == NULL) {
    fprintf(stderr, "Error opening movie index file %s\n", str);
    error = true;
  }

  // The 3.14159 header is to test the endian since this is a binary
  // file.  The root cell spacing is necessary to calculate the dx for
  // all other grids.
  float myPi = 3.14159;
  char nFields8bit = (char) NumberOfAllFields;
  char stag8bit = (char) stag;
  int field;
  fwrite(&myPi, sizeof(float), 1, index_file);
  fwrite(&root_dx, sizeof(float), 1, index_file);
  fwrite(&nFields8bit, sizeof(char), 1, index_file);
  fwrite(&stag8bit, sizeof(char), 1, index_file);
  for (field = 0; field < NumberOfAllFields; field++)
    fwrite(HDF5_FieldNames[field], sizeof(char), 64, index_file);
  
  /* Write a file with all of the topgrid times*/

  FILE *fptr;
  if (writeTimeFile == 1) {
    if (overwrite)
      fptr = fopen("MovieTimes.dat", "w");
    else
      fptr = fopen("MovieTimes.dat", "a");

    assert(fptr != NULL);

    if (overwrite)
      fprintf(fptr, "# %10s %20s %20s\n", "Cycle", "Time", "Redshift");

    fprintf(fptr, "  %10d %20.15f %20.15f\n", cycle, time, redshift);
    fclose(fptr);

  } // ENDIF writeTimeFile

};


AMRHDF5Writer::~AMRHDF5Writer() 
{

};

herr_t AMRHDF5Writer::WriteTextures(  const int    timeStep,
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
				      void   **dataPtr)

/* Write a file with a flat, continuous structure (i.e. not grouped) */

{
  
  char gridDataName[100], fieldName[100];
  sprintf(gridDataName, "/grid-%d", gridId);

  hid_t gridGrp, dataspace, dataset;
  hsize_t hdims[2] = { dims[0], dims[1] };
  int err = 0;
  int size = dims[0]*dims[1];
  int totalSize = 0;
  int i;

  /* Write index entry */

  char level8bit = (char) levelIndex;
  short shortdims[3] = { dims[0], dims[1], dims[2] };

  if (dim == 0) {
    fwrite(&gridId, sizeof(int), 1, index_file);
    fwrite(&physicalTime, sizeof(double), 1, index_file);
    fwrite(&level8bit, sizeof(char), 1, index_file);
    fwrite(integerOrigin, sizeof(Eint64), 3, index_file);
    fwrite(shortdims, sizeof(short), 3, index_file);
  } // ENDIF dim == 0

  /* Write Data */

  // Create or open grid group
  if (dim == 0) {
    gridGrp = H5Gcreate(fileId, gridDataName, 0);

    /* Write attributes */

    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "level", &levelIndex);

    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "timestep", 
				&timeStep);

    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_DOUBLE, "time", 
				&physicalTime);

    err |= writeArrayAttribute(gridGrp, H5T_NATIVE_DOUBLE, 3, "center", 
			       gridCenter);

    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_DOUBLE, 3, "delta",        
				delta);

    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_DOUBLE, 3, "origin",        
				physicalOrigin);

    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_LLONG,  3, "iorigin",  
				integerOrigin);

  }
  else
    err |= checkErr(gridGrp = H5Gopen(fileId, gridDataName), gridDataName);

  assert(err==0);

  for (i = 0; i < nFields; i++) {

    sprintf(fieldName, "%s/%s-dim%d", gridDataName, names[i], dim);

    dataspace = H5Screate_simple(2, hdims, NULL);
    dataset = H5Dcreate(fileId, fieldName, h5DataType, dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataPtr[i]);
    
    H5Dclose(dataset);
    H5Sclose(dataspace);

  } /* ENDFOR: field */

  H5Gclose(gridGrp);
  H5Fflush(fileId, H5F_SCOPE_LOCAL);

  return 0;

}

/***********************************************************************/

void AMRHDF5Writer::AMRHDF5Close() 
{

  //  writeRemainingMetaData();
  
  //  H5Tclose (h5DataType);
  
  // Number of total grids
  fwrite(&gridId, sizeof(int), 1, index_file);
  fclose(index_file);

  H5Fclose (fileId);

};

herr_t AMRHDF5Writer::WriteFlat(  const int    timeStep,
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
				  const void   *dataPtr)

/* Write a file with a flat, continuous structure (i.e. not grouped) */

{
  
  char gridDataName[100], fieldName[100];
  sprintf(gridDataName, "/grid-%d", gridId);
  sprintf(fieldName, "%s/%s", gridDataName, name);

  hid_t gridGrp, dataspace, dataset;
  hsize_t hdims[3] = { dims[2], dims[1], dims[0] };
  int err = 0;
  int size = dims[0] * dims[1] * dims[2];
  char level8bit = (char) levelIndex;
  short shortdims[3] = { dims[0], dims[1], dims[2] };

  /* Write index entry */

  if (fieldNum == 0) {
    fwrite(&gridId, sizeof(int), 1, index_file);
    fwrite(&physicalTime, sizeof(double), 1, index_file);
    fwrite(&timeStep, sizeof(int), 1, index_file);
    fwrite(&redshift, sizeof(double), 1, index_file);
    fwrite(&level8bit, sizeof(char), 1, index_file);
    fwrite(integerOrigin, sizeof(Eint64), 3, index_file);
    fwrite(shortdims, sizeof(short), 3, index_file);
    fwrite(&nParticles, sizeof(int), 1, index_file);
  }

  if (nFields == 0)
    return 0;

  /* Write Data */

  // Create or open grid group
  if (fieldNum == 0) {
    gridGrp = H5Gcreate(fileId, gridDataName, 0);

    /* Write attributes */

    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "level", &levelIndex);

    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "timestep", 
				&timeStep);

    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_DOUBLE, "time", 
				&physicalTime);

    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_DOUBLE, "redshift", 
				&redshift);

    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_INT,    6, "cctk_bbox", 
				bboxflags);

    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_INT,    3, 
				"cctk_nghostzones", nghostzones);

    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_DOUBLE, 3, "delta",        
				delta);

    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_DOUBLE, 3, "origin",        
				physicalOrigin);

    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_LLONG,  3, "iorigin",  
				integerOrigin);

  }
  else
    err |= checkErr(gridGrp = H5Gopen(fileId, gridDataName), gridDataName);

  assert(err==0);

  dataspace = H5Screate_simple(3, hdims, hdims);
  dataset = H5Dcreate(fileId, fieldName, h5DataType, dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataPtr);

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Gclose(gridGrp);
  H5Fflush(fileId, H5F_SCOPE_LOCAL);

  return 0;

}

/***********************************************************************/

herr_t AMRHDF5Writer::writeParticles ( const int nPart,
				       const int nAttributes,
				       const int nBaryonFields,
				       const int Rank,
				       void **pos,
				       void **vel,
				       void *type,
				       void *ID,
				       void *mass,
				       void **attr )
{

  int dim, i;
  int err = 0;
  hid_t gridGrp, dataspace, dataset, PositionDatatype;
  hsize_t hdims = nPart;
  char gridDataName[100];

  const char *ParticlePositionLabel[] = 
     {"particle_position_x", "particle_position_y", "particle_position_z"};
  const char *ParticleVelocityLabel[] = 
     {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
#ifdef WINDS
  const char *ParticleAttributeLabel[] =
    {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x", 
     "particle_jet_y", "particle_jet_z", "typeia_fraction"};
#else
  const char *ParticleAttributeLabel[] = 
    {"creation_time", "dynamical_time", "metallicity_fraction", "typeia_fraction"};
#endif

  sprintf(gridDataName, "/grid-%d", gridId);
  if (nBaryonFields > 0) 
    gridGrp = H5Gopen(fileId, gridDataName);
  else
    gridGrp = H5Gcreate(fileId, gridDataName, 0);

  err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "NumberOfParticles", 
			      &nPart);
  if (nPart == 0) {
    H5Gclose(gridGrp);
    return 0;
  }

  assert(err==0);

  /* Write particle data */

  int FLOAT_Size = sizeof(FLOAT);
  switch (FLOAT_Size) {
  case 4:
    PositionDatatype = H5T_NATIVE_FLOAT;
    break;
  case 8:
    PositionDatatype = H5T_NATIVE_DOUBLE;
    break;
  case 12:
  case 16:
    PositionDatatype = H5T_NATIVE_LDOUBLE;
    break;
  default:
    PositionDatatype = H5T_NATIVE_FLOAT;
  } // ENDSWITCH

  for (dim = 0; dim < Rank; dim++) {

    // Position
    dataspace = H5Screate_simple(1, &hdims, &hdims);
    dataset = H5Dcreate(gridGrp, ParticlePositionLabel[dim], PositionDatatype,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, PositionDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos[dim]);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    // Velocity
    dataspace = H5Screate_simple(1, &hdims, &hdims);
    dataset = H5Dcreate(gridGrp, ParticleVelocityLabel[dim], h5DataType,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel[dim]);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    
  } // ENDFOR dim

  // Type
  dataspace = H5Screate_simple(1, &hdims, &hdims);
  dataset = H5Dcreate(gridGrp, "particle_type", H5T_NATIVE_INT,
		      dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, type);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  // ID
  dataspace = H5Screate_simple(1, &hdims, &hdims);
  dataset = H5Dcreate(gridGrp, "particle_index", HDF5_FILE_PINT,
		      dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, HDF5_PINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ID);
  H5Dclose(dataset);
  H5Sclose(dataspace);


  // Mass
  dataspace = H5Screate_simple(1, &hdims, &hdims);
  dataset = H5Dcreate(gridGrp, "particle_mass", h5DataType,
		      dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  // Attributes
  for (i = 0; i < nAttributes; i++) {
    dataspace = H5Screate_simple(1, &hdims, &hdims);
    dataset = H5Dcreate(gridGrp, ParticleAttributeLabel[i], h5DataType,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, attr[i]);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  } // ENDFOR attributes

  H5Gclose(gridGrp);
  H5Fflush(fileId, H5F_SCOPE_LOCAL);
  
  return 0;

}

/***********************************************************************/

void AMRHDF5Writer::IncreaseGridCount() { gridId++; }













/* The following 2 functions were made for NON_DM_PARTICLES_MERGED_LEVEL
  - Ji-hoon Kim, Apr. 2010 */   

/***********************************************************************/

void AMRHDF5Writer::IncreaseParticleGridCount() { particlegridId++; }



/***********************************************************************/

herr_t AMRHDF5Writer::writeParticles2( const int nPart,
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
				       const int    *nghostzones )
{

  int dim, i, nPart_recorded_here, nPart_recorded_here_new;
  int err = 0;
  hid_t gridGrp, dataspace, dataspace2, dataset, PositionDatatype, attrname;
  herr_t ret;
  hsize_t hdims, hdims2;  
  char gridDataName[100];

  const char *ParticlePositionLabel[] = 
     {"particle_position_x", "particle_position_y", "particle_position_z"};
  const char *ParticleVelocityLabel[] = 
     {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
#ifdef WINDS
  const char *ParticleAttributeLabel[] =
    {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x", 
     "particle_jet_y", "particle_jet_z", "typeia_fraction"};
#else
  const char *ParticleAttributeLabel[] = 
    {"creation_time", "dynamical_time", "metallicity_fraction", "typeia_fraction"};
#endif

  /* if there's no particle, don't bother,
     but if this is a root-level grid, print them anyway --> for Ralf's visualization purpose! */
  if (nPart == 0 && levelIndex != 0)
    return 0;

  sprintf(gridDataName, "/particlegrid-%d", particlegridId);
  fprintf(stdout, "AMRH5writer: nPart = %d (NofSP_OnProcOnLvl = %d), alreadyopenedentry = %d, gridDataName = %s\n",
	  nPart, NumberOfStarParticlesOnProcOnLvlEntry, alreadyopenedentry, gridDataName); //#####

  if (alreadyopenedentry == FALSE) {
    gridGrp = H5Gcreate(fileId_particle, gridDataName, 0);

    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "level", &levelIndex);
    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "timestep", 
				&timeStep);
    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_DOUBLE, "time", 
				&physicalTime);
    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_DOUBLE, "redshift", 
				&redshift);
    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_INT,    6, "cctk_bbox", 
				bboxflags);
    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_INT,    3, 
				"cctk_nghostzones", nghostzones);
    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_DOUBLE, 3, "delta",        
				delta);
    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_DOUBLE, 3, "origin",        
				physicalOrigin);
    err |= writeArrayAttribute( gridGrp, H5T_NATIVE_LLONG,  3, "iorigin",  
				integerOrigin);
    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "NumberOfStarParticles", 
				&NumberOfStarParticles);
    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "NumberOfStarParticlesOnProcOnLvlEntry", 
				&NumberOfStarParticlesOnProcOnLvlEntry);
    err |= writeScalarAttribute(gridGrp, H5T_NATIVE_INT, "nPart_recorded_here", 
				&nPart);
  } else {
    gridGrp = H5Gopen(fileId_particle, gridDataName);

    attrname = H5Aopen_name(gridGrp, "nPart_recorded_here");
    ret  = H5Aread(attrname, H5T_NATIVE_INT, &nPart_recorded_here);

    // increase the recorded number of star particles
    nPart_recorded_here_new = nPart_recorded_here + nPart;
    ret  = H5Awrite(attrname, H5T_NATIVE_INT, &nPart_recorded_here_new);
    H5Aclose(attrname);

    fprintf(stdout, "nPart_recorded_here was %d, is now %d \n",
	    nPart_recorded_here, nPart_recorded_here_new); //#####
  }

  if (nPart == 0) {
    H5Gclose(gridGrp);
    return 0;
  }

  assert(err==0);

  /* Write particle data */

  int FLOAT_Size = sizeof(FLOAT);
  switch (FLOAT_Size) {
  case 4:
    PositionDatatype = H5T_NATIVE_FLOAT;
    break;
  case 8:
    PositionDatatype = H5T_NATIVE_DOUBLE;
    break;
  case 12:
  case 16:
    PositionDatatype = H5T_NATIVE_LDOUBLE;
    break;
  default:
    PositionDatatype = H5T_NATIVE_FLOAT;
  } // ENDSWITCH

  if (alreadyopenedentry == FALSE) {
    
    //When creating the group, the size is NumberOfStarParticlesOnProcEntry
    hdims = NumberOfStarParticlesOnProcOnLvlEntry;  
  
    // setting the maximum size to be NULL is critical, 
    // in order to extend the size later with H5Dextend
    dataspace = H5Screate_simple(1, &hdims, NULL);

    for (dim = 0; dim < Rank; dim++) {
      
      // Position
      dataset = H5Dcreate(gridGrp, ParticlePositionLabel[dim], PositionDatatype,
			  dataspace, H5P_DEFAULT);
      H5Dwrite(dataset, PositionDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos[dim]);
      H5Dclose(dataset);
      
      // Velocity
      dataset = H5Dcreate(gridGrp, ParticleVelocityLabel[dim], h5DataType,
			  dataspace, H5P_DEFAULT);
      H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel[dim]);
      H5Dclose(dataset);
      
    } // ENDFOR dim

    // Type
    dataset = H5Dcreate(gridGrp, "particle_type", H5T_NATIVE_INT,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, type);
    H5Dclose(dataset);
    
    // ID
    dataset = H5Dcreate(gridGrp, "particle_index", HDF5_FILE_PINT,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, HDF5_PINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ID);
    H5Dclose(dataset);
    
    // Mass
    dataset = H5Dcreate(gridGrp, "particle_mass", h5DataType,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass);
    H5Dclose(dataset);
    
    // Attributes
    for (i = 0; i < nAttributes; i++) {
      dataset = H5Dcreate(gridGrp, ParticleAttributeLabel[i], h5DataType,
			  dataspace, H5P_DEFAULT);
      H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, attr[i]);
      H5Dclose(dataset);
    } // ENDFOR attributes

    H5Sclose(dataspace);

//    alreadyopenedentry = TRUE;  //now it is fully open   //#####

  } else {

    // When group alreadyopened, 
    // first extend the filespace, and then append the new particles at the end

    hsize_t file_count[1], file_offset[1];

    file_offset[0] = nPart_recorded_here;
    file_count[0] = nPart;

    hdims = NumberOfStarParticlesOnProcOnLvlEntry;  //file
    hdims2 = nPart;                                 //memory

    dataspace = H5Screate_simple(1, &hdims, NULL);
    dataspace2 = H5Screate_simple(1, &hdims2, NULL);
    ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, file_offset, NULL, file_count, NULL);

    hsize_t a1[1], a2[1];
    int a3;
    H5Sget_simple_extent_dims(dataspace, a1, a2);
    a3 = H5Sget_select_npoints(dataspace);

    for (dim = 0; dim < Rank; dim++) {
      
      // Position
      dataset = H5Dopen(gridGrp, ParticlePositionLabel[dim]);
      H5Dwrite(dataset, PositionDatatype, dataspace2, dataspace, H5P_DEFAULT, pos[dim]);
      H5Dclose(dataset);
      
      // Velocity
      dataset = H5Dopen(gridGrp, ParticleVelocityLabel[dim]);
      H5Dwrite(dataset, h5DataType, dataspace2, dataspace, H5P_DEFAULT, vel[dim]);
      H5Dclose(dataset);
      
    } // ENDFOR dim

    // Type
    dataset = H5Dopen(gridGrp, "particle_type");
    H5Dwrite(dataset, H5T_NATIVE_INT, dataspace2, dataspace, H5P_DEFAULT, type);
    H5Dclose(dataset);
    
    // ID
    dataset = H5Dopen(gridGrp, "particle_index");
    H5Dwrite(dataset, HDF5_PINT, dataspace2, dataspace, H5P_DEFAULT, ID);
    H5Dclose(dataset);
    
    // Mass
    dataset = H5Dopen(gridGrp, "particle_mass");
    H5Dwrite(dataset, h5DataType, dataspace2, dataspace, H5P_DEFAULT, mass);
    H5Dclose(dataset);
    
    // Attributes
    for (i = 0; i < nAttributes; i++) {
      dataset = H5Dopen(gridGrp, ParticleAttributeLabel[i]);
      H5Dwrite(dataset, h5DataType, dataspace2, dataspace, H5P_DEFAULT, attr[i]);
      H5Dclose(dataset);
    } // ENDFOR attributes

    H5Sclose(dataspace);
    H5Sclose(dataspace2);

  } // ENDIF alreadyopenedentry

  H5Gclose(gridGrp);
  H5Fflush(fileId_particle, H5F_SCOPE_LOCAL);

  return 0;

}









/* The following 4 functions were made for NON_DM_PARTICLES_IN_MERGED_ALL
  - Ji-hoon Kim, Apr. 2010 */   

/***********************************************************************/

void AMRHDF5Writer::IncreaseOutputParticleCount() { output_particle++; }



/***********************************************************************/

void AMRHDF5Writer::AMRHDF5CloseSeparateParticles() { H5Fclose (fileId_particle); }



/***********************************************************************/

void AMRHDF5Writer::AMRHDF5CreateSeparateParticles( const char*      fileName, 
						    const int        ParticlesOn,
						    const int        nParticleAttr,
						    bool&            error)
{
  error=false;

  const char *ParticlePositionLabel[] = 
    {"particle_position_x", "particle_position_y", "particle_position_z"};
  const char *ParticleVelocityLabel[] = 
    {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
  const char *ParticleOtherLabel[] =
    {"particle_type", "particle_index", "particle_mass"};
#ifdef WINDS
  const char *ParticleAttributeLabel[] =
    {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x", 
     "particle_jet_y", "particle_jet_z", "typeia_fraction"};
#else
  const char *ParticleAttributeLabel[] = 
    {"creation_time", "dynamical_time", "metallicity_fraction", "typeia_fraction"};
#endif

  int i;
    
  fileId_particle = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if (fileId_particle<0) {
    fprintf(stderr,"Can't open file %s for writing.\n", fileName);
    error = true;
  }

  /* Write time data to a text file */

  /* Create global attributes */

  hid_t groupId;
  int err = 0;

  gridId = 0;
  output_particle = 0;  //This is the counter when recording particles in separate hdf5 file

  /* Create a string array that has the field names */

  char buf[80];
  int NumberOfParticleFields = (ParticlesOn == TRUE) ? 9+nParticleAttr : 0;
  int iField = 0;
  char **HDF5_FieldNames = new char*[NumberOfParticleFields];
  for (i = 0; i < NumberOfParticleFields; i++)
    HDF5_FieldNames[i] = new char[64];

  if (ParticlesOn == TRUE) {
    for (i = 0; i < 3; i++)
      strcpy(HDF5_FieldNames[iField++], ParticlePositionLabel[i]);
    for (i = 0; i < 3; i++)
      strcpy(HDF5_FieldNames[iField++], ParticleVelocityLabel[i]);
    for (i = 0; i < 3; i++)
      strcpy(HDF5_FieldNames[iField++], ParticleOtherLabel[i]);
    for (i = 0; i < nParticleAttr; i++)
      strcpy(HDF5_FieldNames[iField++], ParticleAttributeLabel[i]);
  } // ENDIF ParticlesOn

  /* Write global attributes */

  char str[100];
  hid_t atype;
  strcpy(str, "/Parameters and Global Attributes");
  err |= checkErr(groupId = H5Gcreate(fileId_particle, str, 0), str);
  err |= writeScalarAttribute(groupId, H5T_NATIVE_INT, "NumberOfParticleFields", 
			      &NumberOfParticleFields);

  atype = H5Tcopy(H5T_C_S1);
  err |= H5Tset_size(atype, H5T_VARIABLE);
  err |= writeArrayAttribute(groupId, atype, NumberOfParticleFields, "FieldNames", 
			     HDF5_FieldNames);

  H5Gclose(groupId);

  H5Fflush(fileId_particle, H5F_SCOPE_LOCAL);

};


/***********************************************************************/

herr_t AMRHDF5Writer::writeSeparateParticles ( const int nPart,
					       const int nAttributes,
					       const int Rank,
					       void **pos,
					       void **vel,
					       void *type,
					       void *ID,
					       void *mass,
					       void **attr,
					       double physicalTime,
					       double redshift,
					       int& alreadyopenedentry,
					       int& NumberOfStarParticlesOnProcEntry )
{

  int dim, i, nPart_recorded_here, nPart_recorded_here_new;
  int err = 0;
  herr_t ret;
  hsize_t hdims, hdims2;  
  
  hid_t partGrp, dataspace, dataspace2, dataset, PositionDatatype, attrname;
  char partDataName[100];

  const char *ParticlePositionLabel[] = 
     {"particle_position_x", "particle_position_y", "particle_position_z"};
  const char *ParticleVelocityLabel[] = 
     {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
#ifdef WINDS
  const char *ParticleAttributeLabel[] =
    {"creation_time", "dynamical_time", "metallicity_fraction", "particle_jet_x", 
     "particle_jet_y", "particle_jet_z", "typeia_fraction"};
#else
  const char *ParticleAttributeLabel[] = 
    {"creation_time", "dynamical_time", "metallicity_fraction", "typeia_fraction"};
#endif

  if (nPart == 0) 
    return 0;

  sprintf(partDataName, "/Timestep-%d", output_particle);
  //  fprintf(stdout, "AMRH5writer: nPart = %d, alreadyopenedentry = %d, partDataName = %s\n",
  //	  nPart, alreadyopenedentry, partDataName); 
  //  fprintf(stdout, "fileId_particle = %g\n", fileId_particle);

  if (alreadyopenedentry == FALSE) {
    partGrp = H5Gcreate(fileId_particle, partDataName, 0);
    err |= writeScalarAttribute(partGrp, H5T_NATIVE_INT, "NumberOfStarParticles", 
				&NumberOfStarParticles);
    err |= writeScalarAttribute(partGrp, H5T_NATIVE_INT, "NumberOfStarParticlesOnProcEntry", 
				&NumberOfStarParticlesOnProcEntry);
    err |= writeScalarAttribute(partGrp, H5T_NATIVE_INT, "nPart_recorded_here", 
				&nPart);
    err |= writeScalarAttribute(partGrp, H5T_NATIVE_DOUBLE, "time", 
				&physicalTime);
    err |= writeScalarAttribute(partGrp, H5T_NATIVE_DOUBLE, "redshift", 
				&redshift);    
  } else {
    partGrp = H5Gopen(fileId_particle, partDataName);

    attrname = H5Aopen_name(partGrp, "nPart_recorded_here");
    ret  = H5Aread(attrname, H5T_NATIVE_INT, &nPart_recorded_here);

    // increase the recorded number of star particles
    nPart_recorded_here_new = nPart_recorded_here + nPart;
    ret  = H5Awrite(attrname, H5T_NATIVE_INT, &nPart_recorded_here_new);
    H5Aclose(attrname);

  //  fprintf(stdout, "AMRH5writer: nPart = %d, alreadyopenedentry = %d, partDataName = %s\n"
  //	  "nPart_recorded_here was %d, is now %d \n",
  //	  nPart, alreadyopenedentry, partDataName, 
  //	  nPart_recorded_here, nPart_recorded_here_new); 
  }

  assert(err==0);

  /* Write particle data */

  int FLOAT_Size = sizeof(FLOAT);
  switch (FLOAT_Size) {
  case 4:
    PositionDatatype = H5T_NATIVE_FLOAT;
    break;
  case 8:
    PositionDatatype = H5T_NATIVE_DOUBLE;
    break;
  case 12:
  case 16:
    PositionDatatype = H5T_NATIVE_LDOUBLE;
    break;
  default:
    PositionDatatype = H5T_NATIVE_FLOAT;
  } // ENDSWITCH


  if (alreadyopenedentry == FALSE) {
    
    //When creating the group, the size is NumberOfStarParticlesOnProcEntry
    hdims = NumberOfStarParticlesOnProcEntry;  
  
    //setting the maximum size to be NULL is critical, 
    // in order to extend the size later with H5Dextend
    dataspace = H5Screate_simple(1, &hdims, NULL);

    for (dim = 0; dim < Rank; dim++) {
      
      // Position
      dataset = H5Dcreate(partGrp, ParticlePositionLabel[dim], PositionDatatype,
			  dataspace, H5P_DEFAULT);
      H5Dwrite(dataset, PositionDatatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos[dim]);
      H5Dclose(dataset);
      
      // Velocity
      dataset = H5Dcreate(partGrp, ParticleVelocityLabel[dim], h5DataType,
			  dataspace, H5P_DEFAULT);
      H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel[dim]);
      H5Dclose(dataset);
      
    } // ENDFOR dim

    // Type
    dataset = H5Dcreate(partGrp, "particle_type", H5T_NATIVE_INT,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, type);
    H5Dclose(dataset);
    
    // ID
    dataset = H5Dcreate(partGrp, "particle_index", HDF5_FILE_PINT,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, HDF5_PINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ID);
    H5Dclose(dataset);
    
    // Mass
    dataset = H5Dcreate(partGrp, "particle_mass", h5DataType,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass);
    H5Dclose(dataset);
    
    // Attributes
    for (i = 0; i < nAttributes; i++) {
      dataset = H5Dcreate(partGrp, ParticleAttributeLabel[i], h5DataType,
			  dataspace, H5P_DEFAULT);
      H5Dwrite(dataset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, attr[i]);
      H5Dclose(dataset);
    } // ENDFOR attributes

    H5Sclose(dataspace);

    alreadyopenedentry = TRUE;  //now it is fully open

  } else {

    // When group alreadyopened, 
    // first extend the filespace, and then append the new particles at the end

    hsize_t file_count[1], file_offset[1];

    file_offset[0] = nPart_recorded_here;
    file_count[0] = nPart;

    hdims = NumberOfStarParticlesOnProcEntry;  //file
    hdims2 = nPart;                            //memory

    dataspace = H5Screate_simple(1, &hdims, NULL);
    dataspace2 = H5Screate_simple(1, &hdims2, NULL);
    ret = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, file_offset, NULL, file_count, NULL);

    hsize_t a1[1], a2[1];
    int a3;
    H5Sget_simple_extent_dims(dataspace, a1, a2);
    a3 = H5Sget_select_npoints(dataspace);

    for (dim = 0; dim < Rank; dim++) {
      
      // Position
      dataset = H5Dopen(partGrp, ParticlePositionLabel[dim]);
      H5Dwrite(dataset, PositionDatatype, dataspace2, dataspace, H5P_DEFAULT, pos[dim]);
      H5Dclose(dataset);
      
      // Velocity
      dataset = H5Dopen(partGrp, ParticleVelocityLabel[dim]);
      H5Dwrite(dataset, h5DataType, dataspace2, dataspace, H5P_DEFAULT, vel[dim]);
      H5Dclose(dataset);
      
    } // ENDFOR dim

    // Type
    dataset = H5Dopen(partGrp, "particle_type");
    H5Dwrite(dataset, H5T_NATIVE_INT, dataspace2, dataspace, H5P_DEFAULT, type);
    H5Dclose(dataset);
    
    // ID
    dataset = H5Dopen(partGrp, "particle_index");
    H5Dwrite(dataset, HDF5_PINT, dataspace2, dataspace, H5P_DEFAULT, ID);
    H5Dclose(dataset);
    
    // Mass
    dataset = H5Dopen(partGrp, "particle_mass");
    H5Dwrite(dataset, h5DataType, dataspace2, dataspace, H5P_DEFAULT, mass);
    H5Dclose(dataset);
    
    // Attributes
    for (i = 0; i < nAttributes; i++) {
      dataset = H5Dopen(partGrp, ParticleAttributeLabel[i]);
      H5Dwrite(dataset, h5DataType, dataspace2, dataspace, H5P_DEFAULT, attr[i]);
      H5Dclose(dataset);
    } // ENDFOR attributes

    H5Sclose(dataspace);
    H5Sclose(dataspace2);

  } // ENDIF alreadyopenedentry

  H5Gclose(partGrp);
  H5Fflush(fileId_particle, H5F_SCOPE_LOCAL);

  return 0;

}

