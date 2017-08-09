/***********************************************************************
/
/  DISK GRAVITY UPDATE PARTICLE COM
/
/  written by: 
/  date:       
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h> 
#include <math.h>
#include <stdio.h>
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
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"


int CommunicationAllSumValues(float *Values, int Number);


int DiskGravityUpdateParticleCOM(LevelHierarchyEntry *LevelArray[],
                                   TopGridData &MetaData){


  if (DiskGravityDarkMatterUpdateCOM == 0)
    return SUCCESS;


  // grid-by-grid and local processor COM
  FLOAT gridCOM[MAX_DIMENSION], localCOM[MAX_DIMENSION];
  float gridMass = 0.0, localMass = 0.0;

  // global COM and global mass

  for (int i = 0; i < MAX_DIMENSION; i ++){
    gridCOM[i] = 0.0;
    localCOM[i] = 0.0;
  }

  for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++){
    LevelHierarchyEntry *Temp = LevelArray[level];

    while (Temp != NULL){

      if (Temp->GridData->isLocal()){

          if (Temp->GridData->DiskGravityComputeParticleCOM(gridCOM, gridMass) == FAIL){
            ENZO_FAIL("Error in grid->DiskGravityComputeParticleCOM\n");
          }

          // running weighted sum for local COM
          for (int i = 0; i < MAX_DIMENSION; i ++){
            localCOM[i] += gridCOM[i] * gridMass;
          }

          localMass += gridMass;
      }

      Temp = Temp->NextGridThisLevel;
    } // end while
  } // end hierarchy loop

  // now we communicate
  // currently have total mass on this processor and weighted positions on this processor
  // need to do an all sum of the weighted positions in each dimension and total mass
  // then divide by total mass to get global COM
  CommunicationAllSumValues(localCOM, MAX_DIMENSION);
  //
  CommunicationAllSumValues(&localMass, 1);

  float inv_mass = 1.0 / localMass;
  for (int i = 0; i < MAX_DIMENSION; i ++){
    DiskGravityDarkMatterCOM[i] = localCOM[i] * inv_mass;
  }

  // DONE

  return SUCCESS;
}
