/*------------------------------------------------------------------------
  WRITE MOVIE DATAFILE (Binary)
  By John Wise
  
  Created : 20 Jan 2004
  Modified: 25 Oct 2006

  Purpose : To output grid data into a stacked datafile and to record
  its existence in an index file.

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
            (Feb-Mar 2006) Ji-hoon Kim and I modified this routine to 
	                   write data in an Amira-readable format.
	    (Oct 2006)    Cleaned up, and changed index file back to 
	                  Feb 2005 version.
------------------------------------------------------------------------*/

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
#include "CosmologyParameters.h"
#include "AMRH5writer.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);

//  --JHK notes--
//  import RootResolution from TopGridDimension[0] to calculate levelIdx
//  also had to change : WriteStreamData.C, WriteDataHierarchy.C, Grid.h
//  import MovieTimestepCounter to calculate the timeStep
//  also had to look at : TopGridData.C, WriteStreamData.C,
//  WriteDataHierarchy.C, EvolveLevel.C, SetDefaultGlobalValues.C,
//  Grid.h

int grid::WriteNewMovieData(FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[], 
			    int RootResolution, FLOAT StopTime, 
			    AMRHDF5Writer &AmiraGrid,
			    int lastMovieStep, int TopGridCycle, 
			    int WriteMe, int MovieTimestepCounter, int open, 
			    FLOAT WriteTime, 
			    int alreadyopened[][MAX_DEPTH_OF_HIERARCHY], 
			    int NumberOfStarParticlesOnProcOnLvl[][MAX_DEPTH_OF_HIERARCHY])
{

  static char suffix[] = ".mdat";
  static char partSuffix[] = ".part";

  static char indexSuffix[] = ".idx";

  static long maxEntriesAll = 40000;    // Maximum entries for all processors
  long maxEntries = maxEntriesAll;      // / NumberOfProcessors;

  /* Declarations */
  int i, j, k, dim, field, size=1, allsize=1, vcsize = 1;
  int ret, gridindex, tempindex;
  int ActiveDim[MAX_DIMENSION], instance;
  int dataWritten, StartNewFile = FALSE;
  FLOAT CurrentRedshift, Left[MAX_DIMENSION], Right[MAX_DIMENSION], a = 1, dadt;
  float *temp, *ThisField;
  float *temperature;
  char fileID[3], pid[6];
  const double ln2 = 0.69314718;
  int TemperatureField = NumberOfBaryonFields+1;
  int thislevel;   
  float root_dx = 1.0 / RootResolution;
  
  double doubleTime, doubleRedshift, *doubletemp;
  char *referenceFileName = NULL, *referenceDataPath = NULL;

  FILE *index = NULL; 
  FILE *index2 = NULL;
  FILE *movie = NULL;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Exit if last timestep, but not called from WriteDataHierarchy */
  if (!lastMovieStep && Time >= StopTime)
    return SUCCESS;

  /* Determine whether to output, its instance, and if to increment
     the counter */
  instance = (lastMovieStep) ? 2 : 0;

  // Flag to write data if it's the last timestep or the n-th timestep
  dataWritten = (lastMovieStep || WriteMe);

  // If not n-th timestep, exit
  if (!dataWritten) return SUCCESS;

  /* If outside the region, skip the grid */
  
  for (dim = 0; dim < GridRank; dim++) {
    Left[dim] = max(RegionLeftEdge[dim], GridLeftEdge[dim]);
    Right[dim] = min(RegionRightEdge[dim], GridRightEdge[dim]);
    if (Left[dim] >= Right[dim])
      return SUCCESS;
  }

  int nFields = 0;
  int NeedTemperature = FALSE;
  while (MovieDataField[nFields] != INT_UNDEFINED) {
    if (MovieDataField[nFields] == TEMPERATURE_FIELD) 
      NeedTemperature = TRUE;
    nFields++;
  }
  if (NumberOfBaryonFields <= 0) nFields = 0;

  /* Find the density field. */

  int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);

  /* Get expansion factor */
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      return FAIL;
    }
    CurrentRedshift = (1 + InitialRedshift)/a - 1;
  } else
    CurrentRedshift = -1.0;

  /* Check if we need to start a new moviefile (at top grid timesteps). */

  if (NewMovieDumpNumber < TopGridCycle && !open) {
    //printf("Inside: %d %d %d\n", NewMovieDumpNumber, TopGridCycle, open);
    StartNewFile = TRUE;
    NewMovieDumpNumber = TopGridCycle;

    char **FieldNames = new char*[nFields];
    for (field = 0; field < nFields; field++) {
      FieldNames[field] = new char[64];
      if (MovieDataField[field] != TEMPERATURE_FIELD)
	strcpy(FieldNames[field], DataLabel[MovieDataField[field]]);
      else
	strcpy(FieldNames[field], "Temperature");
    }

    char *AmiraFileName = new char[80];
    int *RefineByArray = new int[3];
    bool error = FALSE;
    hid_t DataType = H5T_NATIVE_FLOAT;
    char Amira_fileID[5], Amira_pid[6];
  
    staggering stag = (MovieVertexCentered) ? VERTEX_CENTERED : CELL_CENTERED;
    fieldtype field_type = SCALAR;
    for (dim = 0; dim < MAX_DIMENSION; dim++) RefineByArray[dim] = RefineBy;

    sprintf(Amira_pid, "_P%3.3d", MyProcessorNumber);
    sprintf(Amira_fileID, "%4.4d", NewMovieDumpNumber);

    strcpy(AmiraFileName, NewMovieName);
    strcat(AmiraFileName, Amira_fileID);
    strcat(AmiraFileName, Amira_pid);
    strcat(AmiraFileName, ".hdf5");
    
    if (Movie3DVolumes > 0) {
      AmiraGrid.AMRHDF5Close();
      AmiraGrid.AMRHDF5Create(AmiraFileName, RefineByArray, 
			      DataType, stag, field_type, TopGridCycle,
			      Time, CurrentRedshift, root_dx, 0, 
			      (MyProcessorNumber == ROOT_PROCESSOR),
			      nFields, 
			      (NewMovieParticleOn > 0 &&
			       NewMovieParticleOn < 3),
			      NumberOfParticleAttributes, FieldNames, 
			      error);
      //printf("NewMovie: Opened movie data file %s\n", AmiraFileName);

      if (error) {
	fprintf(stderr, "Error in AMRHDF5Writer %s.\n", AmiraFileName);
	return FAIL;
      }

      if (NewMovieParticleOn == NON_DM_PARTICLES_MERGED_LEVEL ||
	  NewMovieParticleOn == NON_DM_PARTICLES_MERGED_ALL) {   

	char *AmiraParticleFileName = new char[80];
	strcpy(AmiraParticleFileName, NewMovieName);
	strcat(AmiraParticleFileName, "Particle");
	strcat(AmiraParticleFileName, Amira_fileID);
	strcat(AmiraParticleFileName, Amira_pid);
	strcat(AmiraParticleFileName, ".hdf5");
	
	AmiraGrid.AMRHDF5CloseSeparateParticles();
	AmiraGrid.AMRHDF5CreateSeparateParticles(AmiraParticleFileName, 
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

    for (field = 0; field < nFields; field++)
      delete [] FieldNames[field];
    delete [] FieldNames;
    delete [] AmiraFileName;
    delete [] RefineByArray;

  } /* ENDIF start new movie file */

  /* Subtract ghost cells from dimensions */

  for (dim = 0; dim < GridRank; dim++)
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;

  int ghostzoneFlags[6];
  int numGhostzones[3];
  int AmiraDims[3];
  Eint64 integerOrigin[3];
  double delta[3], doubleGridLeftEdge[3];

  //PART 1

  for (dim=0; dim<GridRank; dim++) {
    size *= ActiveDim[dim];
    vcsize *= ActiveDim[dim]+1;
    allsize *= GridDimension[dim];
  }

  /* Create attributes for this grid's entry */

  thislevel = (int) nint( -log( *CellWidth[0] * RootResolution) / ln2);
  doubleTime = (double) WriteTime;
  doubleRedshift = (double) CurrentRedshift;

  for (i = 0; i < 6; i++)
    ghostzoneFlags[i] = 0;

  int DimensionCorrection = (MovieVertexCentered) ? 1 : 0;

  for (dim = 0; dim < 3; dim++) {
    numGhostzones[dim] = 0;
    integerOrigin[dim] = (Eint64) (GridLeftEdge[dim] / *CellWidth[0]);
    delta[dim] = (double) *CellWidth[0];
    //AmiraDims[2-dim] = ActiveDim[dim];
    AmiraDims[dim] = ActiveDim[dim] + DimensionCorrection;
    doubleGridLeftEdge[dim] = (double) GridLeftEdge[dim];
  }

  if (NeedTemperature == TRUE) {
    temperature = new float[allsize];

    if (this->ComputeTemperatureField(temperature) == FAIL) {
      fprintf(stderr, "Error in grid->ComputeTemperatureField.\n");
      return FAIL;
    }

    BaryonField[TemperatureField] = new float[allsize];
    for (i = 0; i < allsize; i++)
      BaryonField[TemperatureField][i] = temperature[i];

  }

  int field_num;

  if (Movie3DVolumes > 0) {
    for (field = 0; field < nFields; field++) {

      //PART 2

      /* Now output the field to the movie file */

      /* Prepare to write the field to file */

      // Check if the field isn't temperature                
      // Create temp[] : "number of baryon" or "temperature" (?)
      // Copy non-ghost grid points into temp array

      field_num = (MovieDataField[field] == TEMPERATURE_FIELD) ? TemperatureField :
	MovieDataField[field];
      
      if (MovieVertexCentered) {

	/* VERTEX-CENTERED DATA */

	if (this->ComputeVertexCenteredField(field_num) == FAIL) {
	  fprintf(stderr, "Error in grid->ComputeVertexCenteredField.\n");
	  return FAIL;
	}

	temp = new float[vcsize];
	for (i = 0; i < vcsize; i++)
	  temp[i] = InterpolatedField[field_num][i];

      } else {

	/* CELL-CENTERED DATA */

	temp = new float[size];
	ThisField = BaryonField[field_num];
      
	for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
	  for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	    gridindex = (j + k*GridDimension[1]) * GridDimension[0] + 
	      GridStartIndex[0];
	    tempindex = ((k-GridStartIndex[2])*ActiveDim[1] + (j-GridStartIndex[1])) *
	      ActiveDim[0];
	    for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, gridindex++, 
		   tempindex++)
	      temp[tempindex] = ThisField[gridindex];
	  } /* ENDFOR: j */
      } // ENDELSE vertex centered
      
      //PART 3

      char FieldName[100];
      if (MovieDataField[field] != TEMPERATURE_FIELD)
	strcpy(FieldName, DataLabel[MovieDataField[field]]);
      else
	strcpy(FieldName, "Temperature");

      int NumberOfWrittenParticles = 0;
      if (NewMovieParticleOn == ALL_PARTICLES)
	NumberOfWrittenParticles = NumberOfParticles;
      else if (NewMovieParticleOn == NON_DM_PARTICLES)
	for (i = 0; i < NumberOfParticles; i++)
	  if (ParticleType[i] != PARTICLE_TYPE_DARK_MATTER)
	    NumberOfWrittenParticles++;

      if (AmiraGrid.WriteFlat(MovieTimestepCounter, doubleTime, doubleRedshift, 
			      thislevel, delta, doubleGridLeftEdge, 
			      integerOrigin, ghostzoneFlags, numGhostzones, 
			      AmiraDims, field, nFields, NumberOfWrittenParticles, 
			      FieldName, (void*) temp) != 0) {
	fprintf(stderr, "Error in WriteFlat.\n");fflush(stdout);
	return FAIL;
      }
   
      if (MovieVertexCentered) {
	delete [] InterpolatedField[field_num];
	InterpolatedField[field_num] = NULL;
      }

    } /* END: Field output */

    // If no baryon data exists, just write the index file
    if (nFields == 0)
      if (AmiraGrid.WriteFlat(MovieTimestepCounter, doubleTime, doubleRedshift, 
			      thislevel, delta, doubleGridLeftEdge, 
			      integerOrigin, ghostzoneFlags, numGhostzones, 
			      AmiraDims, field, nFields, NumberOfParticles, 
			      NULL, NULL) != 0) {
	fprintf(stderr, "Error in WriteFlat.\n");
	return FAIL;
      }

    delete [] temp;

  } /* ENDIF: 3D Volumes */

  //PART 4

  /*********** Output Particle Data ***********/

  if (NewMovieParticleOn == ALL_PARTICLES) {
    if (AmiraGrid.writeParticles(NumberOfParticles, NumberOfParticleAttributes,
				 NumberOfBaryonFields, GridRank, 
				 (void **) ParticlePosition, 
				 (void **) ParticleVelocity,
				 (void *) ParticleType, (void *) ParticleNumber, 
				 (void *) ParticleMass,
				 (void **) ParticleAttribute) != 0) {
      fprintf(stderr, "Error in AMRHDF5Writer->writeParticles\n");
      return FAIL;
    }
//  fprintf(stdout, "grid::WriteNewMovieData: NumberOfParticles = %d\n", NumberOfParticles); 
//    fprintf(stdout, "ParticleNumber[j][0] = %d", ParticleNumber[j][0]);
  } /* ENDIF: output all particles */

  if (NewMovieParticleOn == NON_DM_PARTICLES) {

    /* Search for non dark matter particles and record their array
       element.  We don't store the data first because we need the
       particle number first. */

    int *NonDMParticleIndices = new int[NumberOfParticles];
    int NumberOfNonDMParticles = 0;
    int ii, iattr;

    FLOAT *TempPosition[3];
    float *TempVelocity[3], *TempMass;
    float *TempAttr[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
    int *TempType;
    PINT *TempNumber;
    
    for (i = 0; i < NumberOfParticles; i++)
      NonDMParticleIndices[i] = -1;
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleType[i] != PARTICLE_TYPE_DARK_MATTER)
	NonDMParticleIndices[NumberOfNonDMParticles++] = i;

    /* Allocate memory */

    if (NumberOfNonDMParticles > 0) {
      for (dim = 0; dim < GridRank; dim++) {
	TempPosition[dim] = new FLOAT[NumberOfNonDMParticles];
	TempVelocity[dim] = new float[NumberOfNonDMParticles];
      }
      TempMass = new float[NumberOfNonDMParticles];
      for (i = 0; i < NumberOfParticleAttributes; i++)
	TempAttr[i] = new float[NumberOfNonDMParticles];
      TempType = new int[NumberOfNonDMParticles];
      TempNumber = new PINT[NumberOfNonDMParticles];
    } // ENDIF non-DM particles > 0

    /* Move non-DM particles into temp arrays */

    for (i = 0; i < NumberOfNonDMParticles; i++) {
      j = NonDMParticleIndices[i];
      for (dim = 0; dim < GridRank; dim++) {
	TempPosition[dim][i] = ParticlePosition[dim][j];
	TempVelocity[dim][i] = ParticleVelocity[dim][j];
      }
      TempMass[i] = ParticleMass[j];
      for (iattr = 0; iattr < NumberOfParticleAttributes; iattr++)
	TempAttr[iattr][i] = ParticleAttribute[iattr][j];
      TempType[i] = ParticleType[j];
      TempNumber[i] = ParticleNumber[j];
    } // ENDFOR non-DM particles

    /* Write */

    if (AmiraGrid.writeParticles(NumberOfNonDMParticles, NumberOfParticleAttributes,
				 NumberOfBaryonFields, GridRank, 
				 (void **) TempPosition, 
				 (void **) TempVelocity,
				 (void *) TempType, (void *) TempNumber, 
				 (void *) TempMass,
				 (void **) TempAttr) != 0) {
      fprintf(stderr, "Error in AMRHDF5Writer->writeParticles\n");
      return FAIL;
    }

    /* Free memory */

    if (NumberOfNonDMParticles > 0) {
      for (dim = 0; dim < GridRank; dim++) {
	delete [] TempPosition[dim];
	delete [] TempVelocity[dim];
      }
      for (i = 0; i < NumberOfParticleAttributes; i++)
	delete [] TempAttr[i];
      delete [] TempMass;
      delete [] TempType;
      delete [] TempNumber;
    } // ENDIF non-DM particles > 0

    delete [] NonDMParticleIndices;

  } // ENDIF Non-DM particles

  AmiraGrid.IncreaseGridCount();
















  /* Print out particles merged by level if requested.
     - Ji-hoon Kim, Apr.2010 */

  if (NewMovieParticleOn == NON_DM_PARTICLES_MERGED_LEVEL) {

    int *NonDMParticleIndices = new int[NumberOfParticles];
    int NumberOfNonDMParticles = 0, filled_upto_here = 0;
    int ii, iattr;

    for (i = 0; i < NumberOfParticles; i++)
      NonDMParticleIndices[i] = -1;
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleType[i] == PARTICLE_TYPE_STAR)
	NonDMParticleIndices[NumberOfNonDMParticles++] = i;


    /* Figure out upto which index the particles are filled by previous grids */

    if (NumberOfNonDMParticles > 0) {
      filled_upto_here = 0;
      while (StarParticlesOnProcOnLvl_Type[ProcessorNumber][filled_upto_here] != 0 &&
	     filled_upto_here < NumberOfStarParticlesOnProcOnLvl[ProcessorNumber][thislevel])
	filled_upto_here++;
    }

    /* Move star particles into temp arrays */

    for (i = 0; i < NumberOfNonDMParticles; i++) {
      j = NonDMParticleIndices[i];
      for (dim = 0; dim < GridRank; dim++) {
	StarParticlesOnProcOnLvl_Position[ProcessorNumber][dim][filled_upto_here+i] 
	  = ParticlePosition[dim][j];
	StarParticlesOnProcOnLvl_Velocity[ProcessorNumber][dim][filled_upto_here+i] 
	  = ParticleVelocity[dim][j];
      }
      StarParticlesOnProcOnLvl_Mass[ProcessorNumber][filled_upto_here+i] = ParticleMass[j];
      for (iattr = 0; iattr < NumberOfParticleAttributes; iattr++)
	StarParticlesOnProcOnLvl_Attr[ProcessorNumber][iattr][filled_upto_here+i] 
	  = ParticleAttribute[iattr][j];
      StarParticlesOnProcOnLvl_Type[ProcessorNumber][filled_upto_here+i] = ParticleType[j];
      StarParticlesOnProcOnLvl_Number[ProcessorNumber][filled_upto_here+i] = ParticleNumber[j];
    } // ENDFOR non-DM particles




    /* Check whether we are doen with collecting the particles */

    if (NumberOfNonDMParticles > 0) {
      filled_upto_here = 0;
      while (StarParticlesOnProcOnLvl_Type[ProcessorNumber][filled_upto_here] != 0 &&
	     filled_upto_here < NumberOfStarParticlesOnProcOnLvl[ProcessorNumber][thislevel])
	filled_upto_here++;
    }



    /*
    int NumStar = 0;
    for (i = 0; i < NumberOfStarParticlesOnProcOnLvl[ProcessorNumber][thislevel]; i++)
      NumStar += StarParticlesOnProcOnLvl_Type[ProcessorNumber][i]/2;

    if (NumStar == NumberOfStarParticlesOnProcOnLvl[ProcessorNumber][thislevel]) {
    */



    /* Write (printing star particles when one level is done on a processor) */

    if (filled_upto_here == NumberOfStarParticlesOnProcOnLvl[ProcessorNumber][thislevel]) {

      alreadyopened[ProcessorNumber][thislevel] = FALSE;  //##### nothing hasn't been written

      if (AmiraGrid.writeParticles2(NumberOfNonDMParticles, NumberOfParticleAttributes,
				    NumberOfBaryonFields, GridRank, 
				    (void **) StarParticlesOnProcOnLvl_Position[ProcessorNumber], 
				    (void **) StarParticlesOnProcOnLvl_Velocity[ProcessorNumber],
				    (void *) StarParticlesOnProcOnLvl_Type[ProcessorNumber], 
				    (void *) StarParticlesOnProcOnLvl_Number[ProcessorNumber], 
				    (void *) StarParticlesOnProcOnLvl_Mass[ProcessorNumber],
				    (void **) StarParticlesOnProcOnLvl_Attr[ProcessorNumber],
				    alreadyopened[ProcessorNumber][thislevel],
				    NumberOfStarParticlesOnProcOnLvl[ProcessorNumber][thislevel],
				    MovieTimestepCounter, doubleTime, doubleRedshift, 
				    thislevel, delta, doubleGridLeftEdge,
				    integerOrigin, ghostzoneFlags, numGhostzones) != 0) {
	fprintf(stderr, "Error in AMRHDF5Writer->writeParticles2\n");
	return FAIL;
      }

    }





    delete [] NonDMParticleIndices;

  } // ENDIF NON_DM_PARTICLES_MERGED_LEVEL

  /* Clean up */

  if (NeedTemperature == TRUE) {
    delete [] temperature;
    if (BaryonField[TemperatureField] != NULL) {
      delete [] BaryonField[TemperatureField];
      BaryonField[TemperatureField] = NULL;
    }
  }
  
  return SUCCESS;

}
