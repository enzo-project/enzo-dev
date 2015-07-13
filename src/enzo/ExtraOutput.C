
//
// ExtraOutput
// 
// Written by: David Collins 
// date      : May 25, 2011.  3:36 pm.  Cold in my office.
// 
// Purpose   : To provide a horrific amount of data.  Calls to ExtraOutput(), thus WriteAllData, 
//             are sprinkled throughout the code,which helps debug events throughout the evolution.
//             Output positions are indicated by an integer, and outputs are of the form 
//             EDXX_YYYY/ExtraOutputXX_YYYY, where XX is the identifier of the dump and YYYY is a counter
//             for that dump location.  
//
//             Output is controlled by the parameter ExtraOutput, a list of integers.  Each integer
//             corresponds to a location in the code.  Each dump is performed every time the location is reached,
//             which will result in qute a lot of data if one is not careful.
//
//             ExtraOuput is, by design, not re-written into the parameter file when data dumps are made,
//             to prevent writing too much data.  
//
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
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		       TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
		       ImplicitProblemABC *ImplicitSolver,
#endif
		       FLOAT WriteTime = -1, int CheckpointDump = FALSE);
 int *output_number = NULL;
int ExtraOutput(int output_flag, LevelHierarchyEntry *LevelArray[],TopGridData *MetaData, int level, ExternalBoundary *Exterior
#ifdef TRANSFER
			  , ImplicitProblemABC *ImplicitSolver
#endif
       ,char * message ){
    //initialize output_number array.
    int n_outputs=100;
    if( output_number == NULL){
        output_number = new int[n_outputs];
        for( int nout=0;nout<n_outputs;nout++){
            output_number[nout]=0;
        }
    }

    //Check for ouptut
    int WriteOut = FALSE;
    for( int i=0; i<MAX_EXTRA_OUTPUTS; i++){
        if( output_flag == ExtraOutputs[i]){
            WriteOut = TRUE;
            break;
        }
    }

    if( WriteOut ){
        fflush(stdout);
        fprintf(stderr,"Extra Output ED%02"ISYM"_%04d %s\n",output_flag,output_number[output_flag],message);
        LevelHierarchyEntry *Temp2 = LevelArray[0];
        while (Temp2->NextGridThisLevel != NULL)
          Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
        //#ifdef USE_HDF5_GROUPS
        sprintf(MetaData->ExtraDumpName,"Extra%02"ISYM"_",output_flag);
        sprintf(MetaData->ExtraDumpDir,"ED%02"ISYM"_",output_flag);
        if (Group_WriteAllData(MetaData->ExtraDumpName, output_number[output_flag]++,
                   Temp2->GridHierarchyEntry, *MetaData, Exterior,
#ifdef TRANSFER
                   ImplicitSolver,
#endif
                   -1, FALSE) == FAIL) {
                ENZO_FAIL("Error in Group_WriteAllData.");
        }
        fflush(stdout);
    }
  return SUCCESS;
}
