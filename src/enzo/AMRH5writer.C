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
#include "AMRH5writer.h"


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
    {"particle_type", "particle_number", "particle_mass"};
  /*  const char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
      "metallicity_fraction", "alpha_fraction", "p5", "p6"}; */
  char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
				    "metallicity_fraction", "particle_jet_x", "particle_jet_y", "particle_jet_z", "alpha_fraction"};

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
  char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
				    "metallicity_fraction", "particle_jet_x", "particle_jet_y", "particle_jet_z", "alpha_fraction"};

  /*  const char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
      "metallicity_fraction", "alpha_fraction", "p5", "p6"}; */

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
  dataset = H5Dcreate(gridGrp, "particle_number", H5T_NATIVE_INT,
		      dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ID);
  H5Dclose(dataset);
  H5Sclose(dataspace);


  // Mass
  dataspace = H5Screate_simple(1, &hdims, &hdims);
  dataset = H5Dcreate(gridGrp, "particle_mass", H5T_NATIVE_FLOAT,
		      dataspace, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  
  // Attributes
  for (i = 0; i < nAttributes; i++) {
    dataspace = H5Screate_simple(1, &hdims, &hdims);
    dataset = H5Dcreate(gridGrp, ParticleAttributeLabel[i], H5T_NATIVE_FLOAT,
			dataspace, H5P_DEFAULT);
    H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, attr[i]);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  } // ENDFOR attributes

  H5Gclose(gridGrp);
  H5Fflush(fileId, H5F_SCOPE_LOCAL);
  
  return 0;

}

/***********************************************************************/

void AMRHDF5Writer::IncreaseGridCount() { gridId++; }

