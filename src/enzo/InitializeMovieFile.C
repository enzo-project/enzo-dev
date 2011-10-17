/***********************************************************************
/
/  INITIALIZE STREAMING DATA FILES
/
/  written by: John Wise
/  date:       July, 2008
/  modified1:
/
/  PURPOSE:  Open and initialize streaming data files.
/
/  RETURNS: ENZO_SUCCESS or FAIL
/
************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "fortran.def"
#include "AMRH5writer.h"

int InitializeMovieFile(TopGridData &MetaData, HierarchyEntry &TopGrid)
{

  if (MovieSkipTimestep == INT_UNDEFINED)
    return SUCCESS;

  /* If streaming movie output, write header file. */

  FILE *header;
  char *headerName = (char*) "movieHeader.dat";
  int sizeOfFLOAT = sizeof(FLOAT);
  int sizeOfRecord = (7+MAX_MOVIE_FIELDS)*sizeof(int) + sizeof(float) + 
    6*sizeof(FLOAT);
  char *movieVersion = (char*) "1.3";
  int nMovieFields = 0, dim;
  while (MovieDataField[nMovieFields] != INT_UNDEFINED &&
 	 nMovieFields < MAX_MOVIE_FIELDS) nMovieFields++;
     
  if ((header = fopen(headerName, "w")) == NULL) {
    fprintf(stderr, "Error in opening movie header.\n");
    return FAIL;
  }
  fprintf(header, "MovieVersion = %s\n", movieVersion);
  fprintf(header, "RootReso = %d\n", MetaData.TopGridDims[0]);
  fprintf(header, "FLOATSize = %d\n", sizeOfFLOAT);
  fprintf(header, "RecordSize = %d\n", sizeOfRecord);
  fprintf(header, "NumFields = %d\n", nMovieFields);
  fprintf(header, "NumCPUs = %d\n", NumberOfProcessors);
  fprintf(header, "FileStem = %s\n", NewMovieName);
  fclose(header);
 
  /* Open Amira Data file, if requested */
 
  char *AmiraFileName = new char[80];
  int *RefineByArray = new int[3];
  bool error = FALSE;
  hid_t DataType = H5T_NATIVE_FLOAT;
  char fileID[4], pid[6];
  float root_dx = 1.0 / MetaData.TopGridDims[0];
   
  staggering stag = CELL_CENTERED;
  fieldtype field_type = SCALAR;
  for (dim = 0; dim < MAX_DIMENSION; dim++) RefineByArray[dim] = RefineBy;
 
  sprintf(pid, "_P%3.3d", MyProcessorNumber);
  sprintf(fileID, "%4.4d", NewMovieDumpNumber);
 
  strcpy(AmiraFileName, NewMovieName);
  strcat(AmiraFileName, fileID);
  strcat(AmiraFileName, pid);
  strcat(AmiraFileName, ".hdf5");

  int field, nBaryonFields;
  int nFields = 0;
  while (MovieDataField[nFields] != INT_UNDEFINED)
    nFields++;

  nBaryonFields = TopGrid.GridData->ReturnNumberOfBaryonFields();

  char **FieldNames = new char*[nFields];
  for (field = 0; field < nFields; field++) {
    FieldNames[field] = new char[64];
      
    if (MovieDataField[field] != TEMPERATURE_FIELD)
      strcpy(FieldNames[field], DataLabel[MovieDataField[field]]);
    else
      strcpy(FieldNames[field], "Temperature");
  }
  
  if (Movie3DVolumes > 0) {
    MetaData.AmiraGrid.AMRHDF5Create(AmiraFileName, RefineByArray, 
				     DataType, stag, field_type, 
				     MetaData.CycleNumber, MetaData.Time, 
				     InitialRedshift, root_dx, 1,
				     (MyProcessorNumber == ROOT_PROCESSOR),
				     nFields, 
				     (NewMovieParticleOn > 0 &&
				      NewMovieParticleOn < 3),
				     NumberOfParticleAttributes, FieldNames, 
				     error);
    if (error) {
      fprintf(stderr, "Error in AMRHDF5Writer.\n");
      return FAIL;
    }

    if (NewMovieParticleOn == NON_DM_PARTICLES_MERGED_LEVEL ||
	NewMovieParticleOn == NON_DM_PARTICLES_MERGED_ALL) {   

      char *AmiraParticleFileName = new char[80];
      strcpy(AmiraParticleFileName, NewMovieName);
      strcat(AmiraParticleFileName, "Particle");
      strcat(AmiraParticleFileName, fileID);
      strcat(AmiraParticleFileName, pid);
      strcat(AmiraParticleFileName, ".hdf5");

      MetaData.AmiraGrid.AMRHDF5CreateSeparateParticles(AmiraParticleFileName, 
							(NewMovieParticleOn > 0),
							NumberOfParticleAttributes,  
							error);
      if (error) {
	fprintf(stderr, "Error in AMRHDF5Writer.\n");
	return FAIL;
      }
      delete [] AmiraParticleFileName;
    }
  }

  delete [] AmiraFileName;
  delete [] RefineByArray;
  for (field = 0; field < nFields; field++)
    delete [] FieldNames[field];

  return SUCCESS;

}
