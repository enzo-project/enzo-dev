/*------------------------------------------------------------------------
  WRITE MOVIE DATAFILE (Binary)
  By John Wise
  
  Created : 20 Jan 2004
  Modified: 27 Feb 2005

  Purpose : To output grid data into a stacked datafile and to record
  its exisitence in an index file.

  Input   : instance = 0 : Created
                       1 : Update
		       2 : Destroyed

  Index file format:
            1. Grid number (int)
	    2. Time (FLOAT)
	    3. dt (float)
	    4. Redshift (FLOAT)
	    5. Grid Dimensions (3 x int)
	    6. Cell Width (FLOAT)
	    7. LeftEdge (3 x FLOAT)
	    8. Data fields (MAX_MOVIE_FIELDS x int)
	    9. Number of particles (int)
	   10. Instance (int)
	   11. Data written? (int)

  History : ( 9 Feb 2004) Attempting to parallelize this routine.
            (17 Feb 2004) Converting to binary format
	    (30 Jul 2004) Writes Particle Data
	    (09 Sep 2004) Writes all child grids when written
	    (01 Oct 2004) Index files are now binary.  Writes all
	                  grids' index files.
	    (27 Feb 2005) Now writes every N-th level L timestep.  No
	                  longer recursive.
------------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);

int grid::WriteNewMovieData(FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[],
			    FLOAT StopTime, int lastMovieStep, int &cycle)
{

  static char suffix[] = ".mdat";
  static char partSuffix[] = ".part";
  static char indexSuffix[] = ".idx";
  static long maxEntriesAll = 40000;    // Maximum entries for all processors
  long maxEntries = maxEntriesAll / NumberOfProcessors;

  /* Declarations */
  int i, j, k, dim, field, size=1, ret, ActiveDim[MAX_DIMENSION], instance;
  int dataWritten;
  FLOAT CurrentRedshift, Left[MAX_DIMENSION], Right[MAX_DIMENSION], a = 1, dadt;
  float *temp;
  char fileID[3], fieldID[2];
  FILE *index = NULL;
  FILE *movie = NULL;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Exit if last timestep, but not called from WriteDataHierarchy */
  if (!lastMovieStep && Time >= StopTime)
    return SUCCESS;

  /* Determine whether to output, its instance, and if to increment
     the counter */
  instance = (TimestepsSinceCreation++ != 0);
  if (lastMovieStep) instance = 2;

  // Flag to write data if it's the last timestep or the n-th timestep
  dataWritten = (lastMovieStep || cycle % MovieSkipTimestep == 0 || 
		 instance == 0);
  
  /* If outside the region, skip the grid */
  
  for (dim = 0; dim < GridRank; dim++) {
    Left[dim] = max(RegionLeftEdge[dim], GridLeftEdge[dim]);
    Right[dim] = min(RegionRightEdge[dim], GridRightEdge[dim]);
    if (Left[dim] >= Right[dim])
      return SUCCESS;
  }

  char pid[MAX_TASK_TAG_SIZE];
  sprintf(pid, "_%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);

  /* Check if we need to start a new moviefiles. */
  if (MovieEntriesPP[ProcessorNumber] > maxEntries) {
    NewMovieDumpNumber++;
    NewMovieEntries = 0;
    for (i=0; i<NumberOfProcessors; i++)
      MovieEntriesPP[i] = 0;
    MaxMovieFilenum = max(MaxMovieFilenum, NewMovieDumpNumber);
  }

  /* Get expansion factor */
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      return FAIL;
    }
  CurrentRedshift = (1 + InitialRedshift)/a - 1;

  /* Subtract ghost cells from dimensions */

  sprintf(fileID, "%3.3d", NewMovieDumpNumber);
  for (dim = 0; dim < GridRank; dim++)
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;

  for (field = 0; field < MAX_MOVIE_FIELDS; field++) {

    /* Break if undefined */
    if (MovieDataField[field] == INT_UNDEFINED) break;
    
    /* Write index only on first field */
    if (field == 0) {

      /* Create index filename */
      char *iname = new char[MAX_NAME_LENGTH];
      strcpy(iname, NewMovieName);
      strcat(iname, fileID);
      strcat(iname, indexSuffix);
      strcat(iname, pid);

      /* Open index */
      if (NewMovieEntries == 0) {
	if ((index = fopen(iname, "w")) == NULL) {
	  fprintf(stderr, "Error opening movie index file %s\n", iname);
	  return FAIL;
	} 
      } else {
	if ((index = fopen(iname, "a")) == NULL) {
	  fprintf(stderr, "Error opening movie index file %s\n", iname);
	  return FAIL;
	} 
      }

      delete [] iname;

      fwrite(&NewMovieEntries, 	 	sizeof(int), 1, index);
      fwrite(&Time, 		 	sizeof(FLOAT), 1, index);
      fwrite(&dtFixed, 		 	sizeof(float), 1, index);
      fwrite(&CurrentRedshift, 	 	sizeof(FLOAT), 1, index);
      fwrite(&ActiveDim, 	 	sizeof(int), 3, index);
      fwrite((void *) CellWidth[0], 	sizeof(FLOAT), 1, index);
      fwrite(&GridLeftEdge, 	 	sizeof(FLOAT), 3, index);
      fwrite(&MovieDataField, 	 	sizeof(int), MAX_MOVIE_FIELDS, index);
      fwrite(&NumberOfParticles, 	sizeof(int), 1, index);
      fwrite(&instance, 	 	sizeof(int), 1, index);
      fwrite(&dataWritten, 		sizeof(int), 1, index);

      fclose(index);

      // If not n-th timestep, exit
      if (!dataWritten) return SUCCESS;

    } /* END: Create index (root_processor) */

    /* Now output the field to the movie file */

    // Determine grid size
    for (dim=0; dim<GridRank; dim++) size *= ActiveDim[dim];

    temp = new float[size];

    /* Prepare to write the field to file */

    // Check if the field isn't temperature
    if (MovieDataField[field] != NumberOfBaryonFields) {
      // Copy non-ghost grid points into temp array
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           + 
		 (j-GridStartIndex[1])*ActiveDim[0]              + 
		 (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
	      BaryonField[MovieDataField[field]]
	      [i + j*GridDimension[0] + 
	       k*GridDimension[0]*GridDimension[1]];
    } else {

    // Separate but very similar procedure for temperature
      float *temperature = new float[size];
      if (this->ComputeTemperatureField(temperature) == FAIL) {
	fprintf(stderr, "Error in grid->ComputeTemperatureField.\n");
	return FAIL;
      }

      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
	    temp[(i-GridStartIndex[0])                           + 
	         (j-GridStartIndex[1])*ActiveDim[0]              + 
	         (k-GridStartIndex[2])*ActiveDim[0]*ActiveDim[1] ] =
	      temperature[(k*GridDimension[1] + j)*GridDimension[0] + i];

      delete [] temperature;
            
    } /* END: temperature */

    /* Now write the grid data */

    // Create movie filename
    char *fname = new char[MAX_NAME_LENGTH];
    sprintf(fieldID, ".%1.1d", field);
    strcpy(fname, NewMovieName);
    strcat(fname, fileID);
    strcat(fname, suffix);
    strcat(fname, fieldID);
    strcat(fname, pid);

    // Open movie file
    if (MovieEntriesPP[ProcessorNumber] == 0) {
      if ((movie = fopen(fname, "wb")) == NULL) {
	fprintf(stderr, "Error opening movie file %s\n", fname);
	return FAIL;
      }
    } else {
      if ((movie = fopen(fname, "ab")) == NULL) {
	fprintf(stderr, "Error opening movie file %s\n", fname);
	return FAIL;
      }
    }

    // Write it and clean up
    fwrite((void*) temp, sizeof(FLOAT), size, movie);
    fclose(movie);
    delete [] temp;
    delete [] fname;

} /* END: Field output */

  /*********** Output Particle Data ***********/

  if (NumberOfParticles > 0 && NewMovieParticleOn > 0) {

    // Create moviefile name
    char *fname = new char[MAX_NAME_LENGTH];
    strcpy(fname, NewMovieName);
    strcat(fname, fileID);
    strcat(fname, partSuffix);
    strcat(fname, pid);

    // Open movie file
    if (MovieEntriesPP[ProcessorNumber] == 0) {

      if ((movie = fopen(fname, "wb")) == NULL) {
	fprintf(stderr, "Error opening movie file %s\n", fname);
	return FAIL;
      }

    } else {
      
      if ((movie = fopen(fname, "ab")) == NULL) {
	fprintf(stderr, "Error opening movie file %s\n", fname);
	return FAIL;
      }

    }
      
    /*** Write the particle position, velocity, and attributes ***/

    // Position, Type, and Particle number
    for (dim=0; dim < GridRank; dim++)
      fwrite((float*) ParticlePosition[dim], sizeof(float), NumberOfParticles, 
	     movie);
    //fwrite((short*) ParticleType, sizeof(short), NumberOfParticles, movie);
    fwrite((int*) ParticleNumber, sizeof(int), NumberOfParticles, movie);
    fwrite((void*) ParticleMass, sizeof(float), NumberOfParticles, movie);

    if (NewMovieParticleOn > 1) {

      // Velocity, ID, & Mass
      for (dim=0; dim < GridRank; dim++)
	fwrite((void*) ParticleVelocity[dim], sizeof(float), NumberOfParticles, 
	       movie);

    }  /* END NewMovieParticleOn > 1 */

    if (NewMovieParticleOn > 2) {

      // Particle Attributes
      for (i=0; i < NumberOfParticleAttributes; i++)
	fwrite((void*) ParticleAttribute[i], sizeof(float), NumberOfParticles, 
	       movie);

    }  /* END NewMovieParticleOn > 2 */

    /*** END Write Particle Data ***/

    // Close the particle movie file
    fclose(movie);
    delete [] fname;

  } /* END: output particles */

  // Reset timestep if written
  //  TimestepsSinceCreation = 1;

  NewMovieEntries++;
  MovieEntriesPP[ProcessorNumber]++;
    
  return SUCCESS;

}
