//
// OutputFromEvolveLevel.C
// 
// Written by: David Collins 
// date      : June 10, 2009, 3:37 pm
// 
// Purpose   : Control various outputs from the EvolveLevel routine.
#include "preincludes.h"
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "ImplicitProblemABC.h"


int WriteTracerParticleData(char *basename, int filenumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime);
//#ifdef USE_HDF5_GROUPS
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		       TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
		       ImplicitProblemABC *ImplicitSolver,
#endif
		       FLOAT WriteTime = -1);
// #else
// int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
//                  TopGridData &MetaData, ExternalBoundary *Exterior,
//#ifdef TRANSFER
//	            ImplicitProblemABC *ImplicitSolver,
//#endif
//                  FLOAT WriteTime = -1);
// #endif
void my_exit(int status);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);

#define TIME_MESSAGING 

EXTERN int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];

int OutputFromEvolveLevel(LevelHierarchyEntry *LevelArray[],TopGridData *MetaData,
			  int level, ExternalBoundary *Exterior
#ifdef TRANSFER
			  , ImplicitProblemABC *ImplicitSolver
#endif
			  ){

  int Write = FALSE, ExitEnzo = FALSE, NumberOfGrids;

  //Do all "bottom of hierarchy" checks
  if (LevelArray[level+1] == NULL){
    
    /* Check for tracer particle output */
    
    if (LevelArray[level]->GridData->ReturnTime() >=
	MetaData->TimeLastTracerParticleDump +
	MetaData->dtTracerParticleDump &&
	MetaData->dtTracerParticleDump > 0.0) {
      MetaData->TimeLastTracerParticleDump += MetaData->dtTracerParticleDump;
      if (WriteTracerParticleData(MetaData->TracerParticleDumpName,
				  MetaData->TracerParticleDumpNumber++,
				  LevelArray, MetaData,
				  LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in WriteTracerParticleData.\n");
	return FAIL;
      }
    }
    
    /* Check for new level output */

    if (MetaData->OutputFirstTimeAtLevel > 0 &&
	level >= MetaData->OutputFirstTimeAtLevel){
      MetaData->OutputFirstTimeAtLevel = level+1;
      Write = TRUE;
    }
 
    // File directed output:
    // Existence of the file outputNow will cause enzo to output the next time the bottom
    //    of the hierarchy is reached.
    // a file subcycleCount will change the number of subcycle skip output
    // a file stopNow will output and then exit enzo.
    
    int outputNow = -1, stopNow = -1, subcycleCount=-1;
    if( FileDirectedOutput == TRUE){
      
      CommunicationBarrier();
      outputNow = access("outputNow", F_OK);
      subcycleCount = access("subcycleCount", F_OK);
      stopNow = access("stopNow", F_OK) ;

      if ( outputNow != -1 ){
	printf("Detected outputNow\n");
	Write = TRUE;
      }

      if( stopNow != -1 ) {
	printf("Detected stopNow\n");
	ExitEnzo = TRUE;
	Write = TRUE;
      }

      /* Check to see if new subcycle information has been given to us */
      
      if ( subcycleCount != -1 ){
	printf("Detected subcycleCount\n");
	
	FILE *fptr; 
	if ((fptr = fopen("subcycleCount", "r")) == NULL) {
	  fprintf(stderr, "Error opening subcycle file subcycleCount.  Continuing.\n");
	}
	else {
	  /* Grab the number of cycles to dump on */
	  char line[MAX_LINE_LENGTH];
	  if (fgets(line, MAX_LINE_LENGTH, fptr) == NULL) {
	    fprintf(stderr, "Error reading subcycle file subcycleCount.  Skipping.\n");
	  } else {
	    sscanf(line, "%"ISYM, &MetaData->SubcycleSkipDataDump);
	    MetaData->SubcycleLastDataDump = MetaData->SubcycleNumber;
	  }
	  fclose(fptr);
	}
      }
      
      
      CommunicationBarrier();
      if (MyProcessorNumber == ROOT_PROCESSOR){
	if( outputNow != -1 )
	  if (unlink("outputNow")) {
	    fprintf(stderr, "Error deleting 'outputNow'\n");
	    return FAIL;
	  }
	if( subcycleCount != -1 )
	  if (unlink("subcycleCount")) {
	    fprintf(stderr, "Error deleting subcycleCount.\n");
	  }
	if( stopNow != -1 )
	  if (unlink("stopNow")) {
	    fprintf(stderr, "Error deleting stopNow\n");
	    return FAIL;
	  } 
      } 
      
      CommunicationBarrier();
      
    }//File Directed Output
    
    /* Check to see if we should start outputting interpolated data based on
       the cycles of the highest level */
    
    if (MetaData->SubcycleNumber >= MetaData->SubcycleLastDataDump +
	MetaData->SubcycleSkipDataDump   &&
	MetaData->SubcycleSkipDataDump > 0) {
      printf("Writing data based on SubcycleDumpSkipping (%"ISYM" %"ISYM" %"ISYM")\n",
	     MetaData->SubcycleNumber, MetaData->SubcycleLastDataDump,
	     MetaData->SubcycleSkipDataDump);
      MetaData->SubcycleLastDataDump += MetaData->SubcycleSkipDataDump;
      Write= TRUE;
    } 
    
    if (MetaData->StopFirstTimeAtLevel > 0 &&
	level >= MetaData->StopFirstTimeAtLevel){
      ExitEnzo = TRUE;
      Write = TRUE;
    }
  }//Finest Level

  FILE *Exit_fptr;

  if( ExitEnzo == TRUE ){
    if (MovieSkipTimestep != INT_UNDEFINED) {
      fprintf(stderr, "Closing movie file.\n");
      MetaData->AmiraGrid.AMRHDF5Close();
    }
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      fprintf(stderr, "Stopping due to request on level %"ISYM"\n", level);
      Exit_fptr = fopen("RunFinished", "w");
      fclose(Exit_fptr);
    }
    my_exit(EXIT_SUCCESS);
  }
  
  if( Write == TRUE ){
    
    LevelHierarchyEntry *Temp2 = LevelArray[0];
    while (Temp2->NextGridThisLevel != NULL)
      Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
    //#ifdef USE_HDF5_GROUPS
    if (Group_WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
			   Temp2->GridHierarchyEntry, *MetaData, Exterior,
#ifdef TRANSFER
			   ImplicitSolver,
#endif
			   LevelArray[level]->GridData->ReturnTime()) == FAIL) {
      fprintf(stderr, "Error in Group_WriteAllData.\n");
      return FAIL;
    }
// #else
//     if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
// 		     Temp2->GridHierarchyEntry, *MetaData, Exterior, 
//#ifdef TRANSFER
//		     ImplicitSolver,
//#endif
// 		     LevelArray[level]->GridData->ReturnTime()) == FAIL) {
//       fprintf(stderr, "Error in WriteAllData.\n");
//       return FAIL;
//     }
// #endif
  }//Write == TRUE

  
  
  return SUCCESS;
}
