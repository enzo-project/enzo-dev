/*------------------------------------------------------------------------
  WRITE STREAMING DATA (Binary)
  By John Wise
  
  Created : 09 Sep 2004
  Modified: 

  Purpose : To walk through the hierarchy and call WriteNewMovieData
            for each grid.

  History : 
------------------------------------------------------------------------*/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
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
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#include "CommunicationUtilities.h"

/****************************** Prototypes ******************************/
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
/************************************************************************/

int WriteStreamData(LevelHierarchyEntry *LevelArray[], int level,
		    TopGridData *MetaData, int *CycleCount, int open=FALSE) 
{

  float MaxDensity = -1e20;
  float root_dx = 1.0 / MetaData->TopGridDims[0];
  FLOAT *pos;
  FLOAT lbbox[] = { huge_number,  huge_number,  huge_number};
  FLOAT rbbox[] = {-huge_number, -huge_number, -huge_number};
  FLOAT Left[3], Right[3];
  int Dims[3], Rank, i, j;

  int ilevel, Zero = FALSE;
  LevelHierarchyEntry *Temp;

  int count;

//  if (debug)
//    printf("Movie: %d %d\n", level, CycleCount[level]);

  if (MovieSkipTimestep == INT_UNDEFINED)
    return SUCCESS;

  if (!MetaData->FirstTimestepAfterRestart)
    if (CycleCount[level] != MovieSkipTimestep) {
      CycleCount[level]++;
      return SUCCESS;
    }

  pos = new FLOAT[MAX_DIMENSION];

  // Set flag to FALSE if no radiative transfer.  If RT, then it'll
  // get reset in RestartPhotons.
  if (MetaData->FirstTimestepAfterRestart) {
    NewMovieDumpNumber++;
    open = TRUE;
//    if (!RadiativeTransfer)
//      MetaData->FirstTimestepAfterRestart = FALSE;
  }

  /* Open file if it's the first timestep */

  if (open == TRUE) {

    FILE *header;
    char *headerName = (char*) "movieHeader.dat";
    int sizeOfFLOAT = sizeof(FLOAT);
    int sizeOfRecord = (7+MAX_MOVIE_FIELDS)*sizeof(int) + sizeof(float) + 
      6*sizeof(FLOAT);
    char *movieVersion = (char*) "1.3";
    int nMovieFields = 0;
    while (MovieDataField[nMovieFields] != INT_UNDEFINED &&
	   nMovieFields < MAX_MOVIE_FIELDS) nMovieFields++;
    
    if (MovieSkipTimestep != INT_UNDEFINED) {
      if ((header = fopen(headerName, "w")) == NULL)
	ENZO_FAIL("Error in opening movie header.\n");
      fprintf(header, "MovieVersion = %s\n", movieVersion);
      fprintf(header, "RootReso = %d\n", MetaData->TopGridDims[0]);
      fprintf(header, "FLOATSize = %d\n", sizeOfFLOAT);
      fprintf(header, "RecordSize = %d\n", sizeOfRecord);
      fprintf(header, "NumFields = %d\n", nMovieFields);
      fprintf(header, "NumCPUs = %d\n", NumberOfProcessors);
      fprintf(header, "FileStem = %s\n", NewMovieName);
      fclose(header);
    } /* END: write movie header file */

    int field, nBaryonFields;
    int nFields = 0;
    while (MovieDataField[nFields] != INT_UNDEFINED)
      nFields++;
    nBaryonFields = LevelArray[0]->GridData->ReturnNumberOfBaryonFields();

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
    char fileID[5], pid[6];
  
    staggering stag = (MovieVertexCentered) ? VERTEX_CENTERED : CELL_CENTERED;
    fieldtype field_type = SCALAR;
    for (int dim = 0; dim < MAX_DIMENSION; dim++) RefineByArray[dim] = RefineBy;

    /* Get expansion factor */

    FLOAT CurrentRedshift, a = 1, dadt;
    if (ComovingCoordinates)
      CosmologyComputeExpansionFactor(MetaData->Time, &a, &dadt);
    CurrentRedshift = (1 + InitialRedshift)/a - 1;

    // Restarts will have a number of -001, not to overwrite any
    // previous data
    //    int RestartFileNumber = -1;

    sprintf(pid, "_P%3.3d", MyProcessorNumber);
    sprintf(fileID, "%4.4d", NewMovieDumpNumber);
    
    strcpy(AmiraFileName, NewMovieName);
    strcat(AmiraFileName, fileID);
    strcat(AmiraFileName, pid);
    strcat(AmiraFileName, ".hdf5");

    if (Movie3DVolumes > 0) {
      MetaData->AmiraGrid.AMRHDF5Create(AmiraFileName, RefineByArray, 
					DataType, stag, field_type, 
					MetaData->CycleNumber, MetaData->Time, 
					CurrentRedshift, root_dx, TRUE, 
					(MyProcessorNumber == ROOT_PROCESSOR),
					nFields, 
					(NewMovieParticleOn > 0 &&
					 NewMovieParticleOn < 3),
					NumberOfParticleAttributes, FieldNames, 
					error);
      if (error)
	ENZO_FAIL("Error in AMRHDF5Writer.\n");

      if (NewMovieParticleOn == NON_DM_PARTICLES_MERGED_LEVEL ||
	  NewMovieParticleOn == NON_DM_PARTICLES_MERGED_ALL) {   

	char *AmiraParticleFileName = new char[80];
	strcpy(AmiraParticleFileName, NewMovieName);
	strcat(AmiraParticleFileName, "Particle");
	strcat(AmiraParticleFileName, fileID);
	strcat(AmiraParticleFileName, pid);
	strcat(AmiraParticleFileName, ".hdf5");
	
	MetaData->AmiraGrid.AMRHDF5CreateSeparateParticles(AmiraParticleFileName, 
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

  } // ENDIF open == TRUE

  /* Get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, MetaData->Time);

  /* Write all grids at (1) l<level with radiation (their temperatures
     can change quickly!) or (2) l=level and finer */

  int FineGrainWriting, ilvl, StartLevel = level, WriteMe = TRUE;
  FLOAT WriteTime = LevelArray[level]->GridData->ReturnTime();
  int IonFrontPresent[MAX_DEPTH_OF_HIERARCHY];

  // Only write coarse grids with ionization fronts if the movie data
  // has something other than density or metallicity
  FineGrainWriting = FALSE;
  for (i = 0; i < MAX_MOVIE_FIELDS; i++) {
    if (MovieDataField[i] == INT_UNDEFINED) break;
    if (!(MovieDataField[i] == Density || MovieDataField[i] == SNColour ||
	  MovieDataField[i] == Metallicity)) {
      FineGrainWriting = TRUE;
      break;
    } // ENDIF not Density or Metallicity
  } // ENDFOR movie fields

  // First search grids with l<level for ionization fronts
#ifdef TRANSFER
  if (RadiativeTransfer && level > 0 && FineGrainWriting) {
    for (ilvl = 0; ilvl < level; ilvl++) {
      IonFrontPresent[ilvl] = FALSE;
      Temp = LevelArray[ilvl];
      while (Temp != NULL) {
	if (Temp->GridData->DetectIonizationFrontApprox(TemperatureUnits) == TRUE) {
	  IonFrontPresent[ilvl] = TRUE;
	  break;
	} else
	  Temp = Temp->NextGridThisLevel;
      } // ENDWHILE grids
    } // ENDFOR level

    for (ilvl = 0; ilvl < level; ilvl++)
      if (IonFrontPresent[ilvl] == TRUE) {
	StartLevel = ilvl;
	break;
      }
    
#ifdef USE_MPI
    CommunicationAllReduceValues(&StartLevel, 1, MPI_MIN);
#endif

  } // ENDIF RadiativeTransfer  
#endif /* TRANSFER */

  int NumberOfStarParticlesOnProcOnLvl[NumberOfProcessors][MAX_DEPTH_OF_HIERARCHY];
  int alreadyopened[NumberOfProcessors][MAX_DEPTH_OF_HIERARCHY];
  
  for (i = 0; i < NumberOfProcessors; i++) 
    for (j = 0; j < MAX_DEPTH_OF_HIERARCHY; j++) {
      alreadyopened[i][j] = FALSE;  //whether the group for a certain timestep is already opened
      NumberOfStarParticlesOnProcOnLvl[i][j] = 0;
    }

  /* Print out particles merged by level if requested.  
     This is to speed up the access to the particle data 
     when only particles are being drawn.  - Ji-hoon Kim, Apr.2010 */  
    
  if (NewMovieParticleOn == NON_DM_PARTICLES_MERGED_LEVEL) {
    
    // collect the number of star particles on each processor
    for (ilvl = StartLevel; ilvl < MAX_DEPTH_OF_HIERARCHY; ilvl++) {
      
      Temp = LevelArray[ilvl];
      
      while (Temp != NULL) {
	
	/* Write data */
	if (Temp->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
	  NumberOfStarParticlesOnProcOnLvl[Temp->GridData->ReturnProcessorNumber()][ilvl]
	    += Temp->GridData->ReturnNumberOfStarParticles();
	}

	Temp = Temp->NextGridThisLevel;

      } 
    } 

#ifdef USE_MPI
    CommunicationBarrier();
#endif

    /*
    for (i = 0; i < NumberOfProcessors; i++)  
      for (j = 0; j < 8; j++) 
 	fprintf(stdout, "WriteStreamData: NumberOfStarParticlesOnProcOnLvl[%d][%d] = %d\n", 
 		i, j, NumberOfStarParticlesOnProcOnLvl[i][j]); 
    */
    
  } // ENDIF NON_DM_PARTICLES_MERGED_LEVEL

  if (debug)
    printf("WriteStreamData: level = %d, StartLevel = %d, timestep = %d\n", 
	   level, StartLevel, MetaData->MovieTimestepCounter);

  for (ilvl = StartLevel; ilvl < MAX_DEPTH_OF_HIERARCHY; ilvl++) {

#ifdef USE_MPI  
    CommunicationBarrier();
#endif

    Temp = LevelArray[ilvl];

    /* Allocate memory */

    for (int l = 0; l < NumberOfProcessors; l++) {
      if (NumberOfStarParticlesOnProcOnLvl[l][ilvl] > 0) {
	for (int dim = 0; dim < 3; dim++) {
	  StarParticlesOnProcOnLvl_Position[l][dim] = new FLOAT[NumberOfStarParticlesOnProcOnLvl[l][ilvl]];
	  StarParticlesOnProcOnLvl_Velocity[l][dim] = new float[NumberOfStarParticlesOnProcOnLvl[l][ilvl]];
	}
	StarParticlesOnProcOnLvl_Mass[l] = new float[NumberOfStarParticlesOnProcOnLvl[l][ilvl]];
	for (int nattr = 0; nattr < NumberOfParticleAttributes; nattr++)
	  StarParticlesOnProcOnLvl_Attr[l][nattr] = new float[NumberOfStarParticlesOnProcOnLvl[l][ilvl]];
	StarParticlesOnProcOnLvl_Type[l] = new int[NumberOfStarParticlesOnProcOnLvl[l][ilvl]];
	StarParticlesOnProcOnLvl_Number[l] = new PINT[NumberOfStarParticlesOnProcOnLvl[l][ilvl]];

	// to use it later to determine whether all particles are collected
	for (i = 0; i < NumberOfStarParticlesOnProcOnLvl[l][ilvl]; i++)
	  StarParticlesOnProcOnLvl_Type[l][i] = 0;
      }
    }



    while (Temp != NULL) {

      /* Write data */
      Temp->GridData->WriteNewMovieData
	(MetaData->NewMovieLeftEdge, MetaData->NewMovieRightEdge, 
	 MetaData->TopGridDims[0], MetaData->StopTime, MetaData->AmiraGrid, 
	 Zero, MetaData->CycleNumber, WriteMe, 
	 MetaData->MovieTimestepCounter, open, WriteTime,
	 alreadyopened, NumberOfStarParticlesOnProcOnLvl);

#define NOFIND_DENSEST
#ifdef FIND_DENSEST
      Temp->GridData->FindMaximumBaryonDensity(&MaxDensity, pos);

      if (MyProcessorNumber == ROOT_PROCESSOR) {
	Temp->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
	for (i = 0; i < MAX_DIMENSION; i++) {
	  lbbox[i] = MIN(lbbox[i], Left[i]);
	  rbbox[i] = MAX(rbbox[i], Right[i]);
	}
      }

#endif /* FIND_DENSEST */

      Temp = Temp->NextGridThisLevel;

    } /* ENDWHILE: grid loop */

#ifdef USE_MPI
    CommunicationBarrier();
#endif


    /* Free memory */

    for (int l = 0; l < NumberOfProcessors; l++) {
      if (NumberOfStarParticlesOnProcOnLvl[l][ilvl] > 0) {
	for (int dim = 0; dim < 3; dim++) {
	  delete [] StarParticlesOnProcOnLvl_Position[l][dim];
	  delete [] StarParticlesOnProcOnLvl_Velocity[l][dim];
	}
	for (int nattr = 0; nattr < NumberOfParticleAttributes; nattr++)
	  delete [] StarParticlesOnProcOnLvl_Attr[l][nattr];
	delete [] StarParticlesOnProcOnLvl_Mass[l];
	delete [] StarParticlesOnProcOnLvl_Type[l];
	delete [] StarParticlesOnProcOnLvl_Number[l];
      } 
    }



    /* Increase particlegridId for NON_DM_PARTICLES_MERGED_LEVEL
       when a level is done with nonzero star particles 
       (the group name should be the same on a level) 
       now this is done AmiraGrid.writeParticles2 */

    if (NewMovieParticleOn == NON_DM_PARTICLES_MERGED_LEVEL) 
      MetaData->AmiraGrid.IncreaseParticleGridCount();  

#ifdef USE_MPI 
    CommunicationBarrier();
#endif

  } /* ENDFOR: level loop */











  /* Print out particles all merged to a separate file if requested.  
     This is to speed up the access to the particle data 
     when only particles are being drawn.  - Ji-hoon Kim, Apr.2010 */  

  if (NewMovieParticleOn == NON_DM_PARTICLES_MERGED_ALL) {
    
    int NumberOfStarParticlesOnProc[NumberOfProcessors];
    int alreadyopened[NumberOfProcessors];

    for (i = 0; i < NumberOfProcessors; i++) {
      alreadyopened[i] = FALSE;  //whether the group for a certain timestep is already opened
      NumberOfStarParticlesOnProc[i] = 0;
    }

    // collect the number of star particles on each processor
    for (ilvl = 0; ilvl < MAX_DEPTH_OF_HIERARCHY; ilvl++) {
      
      Temp = LevelArray[ilvl];
      
      while (Temp != NULL) {
	
	/* Write data */
	if (Temp->GridData->ReturnProcessorNumber() == MyProcessorNumber) {
	  NumberOfStarParticlesOnProc[Temp->GridData->ReturnProcessorNumber()]
	    += Temp->GridData->ReturnNumberOfStarParticles();
	}

	Temp = Temp->NextGridThisLevel;

      } 
    } 

#ifdef USE_MPI
    CommunicationBarrier();
#endif

//    for (i = 0; i < NumberOfProcessors; i++)
//      fprintf(stdout, "WriteStreamData: NumberOfStarParticlesOnProc[%d] = %d\n", 
//	      i, NumberOfStarParticlesOnProc[i]);  

    // starting from ilvl=0, to print all the particles in the dataset
    for (ilvl = 0; ilvl < MAX_DEPTH_OF_HIERARCHY; ilvl++) {
      
      Temp = LevelArray[ilvl];
      
      while (Temp != NULL) {
	
	/* Write data */
	Temp->GridData->WriteNewMovieDataSeparateParticles
	  (MetaData->NewMovieLeftEdge, MetaData->NewMovieRightEdge, 
	   MetaData->StopTime, MetaData->AmiraGrid, 
	   Zero, WriteMe, WriteTime, 
	   alreadyopened, NumberOfStarParticlesOnProc);

	Temp = Temp->NextGridThisLevel;

      } 
    } 

#ifdef USE_MPI
    CommunicationBarrier();
#endif

    MetaData->AmiraGrid.IncreaseOutputParticleCount();

#ifdef USE_MPI  
    CommunicationBarrier();
#endif

  }

#ifdef FIND_DENSEST
  float value = MaxDensity;
  int haveMax = FALSE;
  int procWithMax = -1;
  int tempProcNum = 0;

#ifdef USE_MPI
  CommunicationAllReduceValues(&value, 1, MPI_MAX);
  if (fabs(value-MaxDensity) < 1e-5) haveMax = TRUE;
  MaxDensity = value;

  if (haveMax == FALSE)
    for (i = 0; i < MAX_DIMENSION; i++)
      pos[i] = -huge_number;

  CommunicationReduceValues(pos, 3, MPI_MAX);
#endif /* USE_MPI */    

  FILE *fptr;
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if (open == TRUE) {
      if ((fptr = fopen("densestPoint.dat", "w")) == NULL)
	ENZO_FAIL("Error in opening file densestPoint.dat\n");
      fprintf(fptr, "# Time   Density   Densest_Point    Left_Bounding_Box    Right_Bounding_Box\n");
    } else
      if ((fptr = fopen("densestPoint.dat", "a")) == NULL)
	ENZO_FAIL("Error in opening file densestPoint.dat\n");

    fprintf(fptr, "%"GOUTSYM" %15.6g %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
	    MetaData->Time, MaxDensity, pos[0], pos[1], pos[2],
	    lbbox[0], lbbox[1], lbbox[2], rbbox[0], rbbox[1], rbbox[2]);
    fclose(fptr);
  } // ENDIF ROOT PROCESSOR

#endif  /* FIND_DENSEST */

  delete [] pos;

  // Increase movie cycle count
  for (i = level; i < MAX_DEPTH_OF_HIERARCHY; i++)
    CycleCount[i] = 0;
  CycleCount[level]++;

  return SUCCESS;

}
