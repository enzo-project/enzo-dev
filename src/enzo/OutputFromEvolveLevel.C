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

#include "ErrorExceptions.h"
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
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif


int WriteTracerParticleData(char *basename, int filenumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime);
//#ifdef USE_HDF5_GROUPS
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		       TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
		       ImplicitProblemABC *ImplicitSolver,
#endif
		       FLOAT WriteTime = -1, int CheckpointDump = FALSE);
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
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
void CommunicationBroadcastValues(int *Value, int Number, int BroadcastProcessor);

#define TIME_MESSAGING 

EXTERN int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];

int OutputFromEvolveLevel(LevelHierarchyEntry *LevelArray[],TopGridData *MetaData,
			  int level, ExternalBoundary *Exterior, int OutputNow
#ifdef TRANSFER
			  , ImplicitProblemABC *ImplicitSolver
#endif
			  ){

  int WriteOutput = FALSE, ExitEnzo = FALSE, NumberOfGrids;
  int PackedStatus = 0;
  int CheckpointDump = FALSE;
  WriteOutput = OutputNow;

  //Do all "bottom of hierarchy" checks
  if (LevelArray[level+1] == NULL){
    
    /* Check for tracer particle output */
    
    if (LevelArray[level]->GridData->ReturnTime() >=
	MetaData->TimeLastTracerParticleDump + MetaData->dtTracerParticleDump &&
	MetaData->dtTracerParticleDump > 0.0) {
      MetaData->TimeLastTracerParticleDump += MetaData->dtTracerParticleDump;
      if (WriteTracerParticleData(MetaData->TracerParticleDumpName,
				  MetaData->TracerParticleDumpNumber++,
				  LevelArray, MetaData,
				  LevelArray[level]->GridData->ReturnTime()) == FAIL) {
		ENZO_FAIL("Error in WriteTracerParticleData.");
      }
    }
    
    /* Check for new level output */

    if (MetaData->OutputFirstTimeAtLevel > 0 &&
	level >= MetaData->OutputFirstTimeAtLevel){
      MetaData->OutputFirstTimeAtLevel = level+1;
      WriteOutput = TRUE;
    }

    if (OutputOnDensity == 1 ||
        StopFirstTimeAtDensity > 0. ||
        StopFirstTimeAtMetalEnrichedDensity > 0.) {

      /* Get our units, but only if we need to. */
      float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
            TimeUnits = 1, VelocityUnits = 1;
      FLOAT Time = LevelArray[level]->GridData->ReturnTime();
      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
      }

      /* Make sure we are all synced up across processors. */
      CurrentMaximumDensity = CommunicationMaxValue(CurrentMaximumDensity);
      if (MyProcessorNumber == ROOT_PROCESSOR) {
        fprintf(stderr, "Current maximum density is %"GSYM" g/cm^3.\n",
                (CurrentMaximumDensity*DensityUnits));
      }

      if (StopFirstTimeAtMetalEnrichedDensity > 0.) {
        CurrentMaximumMetalEnrichedDensity = CommunicationMaxValue(CurrentMaximumMetalEnrichedDensity);
        if (MyProcessorNumber == ROOT_PROCESSOR) {
          fprintf(stderr, "Current maximum metal enriched density is %"GSYM" g/cm^3.\n",
                  (CurrentMaximumMetalEnrichedDensity*DensityUnits));
        }
      }

      if (OutputOnDensity == 1) {
        if(log10(CurrentMaximumDensity*DensityUnits) > CurrentDensityOutput) {
          while (log10(CurrentMaximumDensity*DensityUnits) > CurrentDensityOutput) {
            CurrentDensityOutput += IncrementDensityOutput;
          }
          fprintf(stderr, "Outputting based on DensMax == %"FSYM" (now set to %"FSYM")\n",
                  log10(CurrentMaximumDensity*DensityUnits), CurrentDensityOutput);
          WriteOutput = TRUE;
        }
      }

      if (StopFirstTimeAtDensity > 0. && 
          CurrentMaximumDensity*DensityUnits >= StopFirstTimeAtDensity) {
        if (MyProcessorNumber == ROOT_PROCESSOR) {
          fprintf(stderr, "Exiting after reaching max density of %"GSYM" g/cm^3.\n",
                  (StopFirstTimeAtDensity));
        }
        ExitEnzo = TRUE;
        WriteOutput = TRUE;
      }

      if (StopFirstTimeAtMetalEnrichedDensity > 0. && 
          CurrentMaximumMetalEnrichedDensity*DensityUnits >= StopFirstTimeAtMetalEnrichedDensity) {
        if (MyProcessorNumber == ROOT_PROCESSOR) {
          fprintf(stderr, "Exiting after reaching max density of %"GSYM" g/cm^3.\n",
                  (StopFirstTimeAtMetalEnrichedDensity));
        }
        ExitEnzo = TRUE;
        WriteOutput = TRUE;
      }

    } // end: if (OutputOnDensity == 1 || ...

    // File directed output:
    // Existence of the file outputNow will cause enzo to output the next time the bottom
    //    of the hierarchy is reached.
    // a file subcycleCount will change the number of subcycle skip output
    // a file stopNow will output and then exit enzo.
    
    int outputNow = -1, stopNow = -1, subcycleCount=-1, checkpointDumpNow=-1;
    if( FileDirectedOutput == TRUE){

    int OutputFlagArray[5] = {-1, -1, -1, -1, -1};

    if(MyProcessorNumber == ROOT_PROCESSOR) {
      
      outputNow = access("outputNow", F_OK);
      subcycleCount = access("subcycleCount", F_OK);
      stopNow = access("stopNow", F_OK) ;
      
      checkpointDumpNow = access("checkpointDump", F_OK);

      if ( outputNow != -1 ){
	printf("Detected outputNow\n");
	WriteOutput = TRUE;
      }

      if( stopNow != -1 ) {
	printf("Detected stopNow\n");
	ExitEnzo = TRUE;
	WriteOutput = TRUE;
      }

      if( checkpointDumpNow != -1 ) {
	//ExitEnzo = TRUE;
	WriteOutput = TRUE;
	CheckpointDump = TRUE;
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
      
      if( outputNow != -1 )
	if (unlink("outputNow")) {
	  ENZO_FAIL("Error deleting 'outputNow'");
	}
      if( subcycleCount != -1 )
	if (unlink("subcycleCount")) {
	  fprintf(stderr, "Error deleting subcycleCount.\n");
	}
      if( stopNow != -1 )
	if (unlink("stopNow")) {
	  ENZO_FAIL("Error deleting stopNow");
	} 
      if( checkpointDumpNow != -1 )
	if (unlink("checkpointDump")) {
	  ENZO_FAIL("Error deleting checkpointDump");
	} 

    // Now we fill the array we're going to broadcast
    OutputFlagArray[0] = WriteOutput;
    OutputFlagArray[1] = ExitEnzo;
    OutputFlagArray[2] = CheckpointDump;
    OutputFlagArray[3] = MetaData->SubcycleSkipDataDump;
    OutputFlagArray[4] = MetaData->SubcycleLastDataDump;
    
    }//Root Processor only

    /* This packs up a broadcast for all the flags */

    CommunicationBroadcastValues(OutputFlagArray, 5, ROOT_PROCESSOR);
    WriteOutput = OutputFlagArray[0];
    ExitEnzo = OutputFlagArray[1];
    CheckpointDump = OutputFlagArray[2];
    MetaData->SubcycleSkipDataDump = OutputFlagArray[3];
    MetaData->SubcycleLastDataDump = OutputFlagArray[4];

    }//File Directed Output
    
    /* We also reset checkpoint state here */
    if (CheckpointRestart == TRUE) CheckpointRestart = FALSE;


    /* Check to see if we should start outputting interpolated data based on
       the time passed (dtInterpolatedDataDump < dtDataDump).
       This is mostly for making movies or looking at the interim data where TopGrid dt is too long.
       In principle, this output shouldn't be used for restart. */

    if (LevelArray[level]->GridData->ReturnTime() >= 
	MetaData->TimeLastInterpolatedDataDump + MetaData->dtInterpolatedDataDump   && 
	MetaData->dtInterpolatedDataDump > 0.0) {
      printf("Writing data based on dtInterpolatedDataDump (%"FSYM" %"FSYM" %"FSYM")\n",
	     LevelArray[level]->GridData->ReturnTime(), MetaData->TimeLastInterpolatedDataDump,
	     MetaData->dtInterpolatedDataDump);
      MetaData->TimeLastInterpolatedDataDump += MetaData->dtInterpolatedDataDump;
      WriteOutput = TRUE;
    }

    /* Check to see if we should start outputting interpolated data based on
       the cycles of the highest level */
    
    if (MetaData->SubcycleNumber >= MetaData->SubcycleLastDataDump +
	MetaData->SubcycleSkipDataDump   &&
	MetaData->SubcycleSkipDataDump > 0) {
      printf("Writing data based on SubcycleDumpSkipping (%"ISYM" %"ISYM" %"ISYM")\n",
	     MetaData->SubcycleNumber, MetaData->SubcycleLastDataDump,
	     MetaData->SubcycleSkipDataDump);
      MetaData->SubcycleLastDataDump += MetaData->SubcycleSkipDataDump;
      WriteOutput= TRUE;
    } 
    
    if (MetaData->StopFirstTimeAtLevel > 0 &&
	level >= MetaData->StopFirstTimeAtLevel){
      ExitEnzo = TRUE;
      WriteOutput = TRUE;
    }


  }//Finest Level

  FILE *Exit_fptr;

  if( WriteOutput == TRUE ){    
    LevelHierarchyEntry *Temp2 = LevelArray[0];
    while (Temp2->NextGridThisLevel != NULL)
      Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
    //#ifdef USE_HDF5_GROUPS
    if (Group_WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
			   Temp2->GridHierarchyEntry, *MetaData, Exterior,
#ifdef TRANSFER
			   ImplicitSolver,
#endif
			   LevelArray[level]->GridData->ReturnTime(), CheckpointDump) == FAIL) {
            ENZO_FAIL("Error in Group_WriteAllData.");
    }
// #else
//     if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
// 		     Temp2->GridHierarchyEntry, *MetaData, Exterior, 
// #ifdef TRANSFER
// 		     ImplicitSolver,
// #endif
// 		     LevelArray[level]->GridData->ReturnTime()) == FAIL) {
//       ENZO_FAIL("Error in WriteAllData.\n");
//     }
// #endif
  }//WriteOutput == TRUE

  if( ExitEnzo == TRUE ){
    if (MovieSkipTimestep != INT_UNDEFINED) {
      fprintf(stderr, "Closing movie file.\n");
      MetaData->AmiraGrid.AMRHDF5Close();
      MetaData->AmiraGrid.AMRHDF5CloseSeparateParticles();
    }
    if (MyProcessorNumber == ROOT_PROCESSOR) {

      fprintf(stderr, "Stopping due to request on level %"ISYM"\n", level);
      Exit_fptr = fopen("RunFinished", "w");
      fclose(Exit_fptr);
    }
    my_exit(EXIT_SUCCESS);
  }
  
  return SUCCESS;
}
