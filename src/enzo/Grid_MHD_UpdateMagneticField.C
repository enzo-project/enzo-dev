/***********************************************************************
/
/  GRID CLASS (UPDATE MAGNETIC FIELDS)
/
/  written by: David Collins
/  date:       2005?
/  modified1:
/
/  PURPOSE: Takes the curl a second time, after UpdateFromFinerGrids
/           replaces ElectricField with the fine grid data.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "ErrorExceptions.h"
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int grid::MHD_UpdateMagneticField(int level, LevelHierarchyEntry * NextLevel){

   if(MyProcessorNumber != ProcessorNumber || UseMHDCT != TRUE)
    return SUCCESS;

  //if we're projecting the magnetic field, this routine isn't necessary.
  //Note that one should never project the magnetic field.
  if( MHD_ProjectE != TRUE )
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, field;
  float dtUsed;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
  
  //dtUsed is set to 1 here because the field ElectricField is NOT E, but E/dt. (This is to ensure
  // that the field is properly projected to parents over subgrid timesteps of potentially uneven lengths.)
  // I left in the option of dt != 1 in order to maintain the versitility of the code: might come out later.

  dtUsed = 1.0;

   int CurlStart[3] = {0,0,0}, CurlEnd[3] = {GridDimension[0]-1,GridDimension[1]-1,GridDimension[2]-1};
   MHD_Curl( CurlStart,CurlEnd, 2);

  CenterMagneticField();
      
    int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  
  //Also update AvgElectricField, if this is a subgrid.
  if( level > 0 ){
    
    
    for(int field=0;field<3;field++){
      
      if(AvgElectricField[field] == NULL)
        ENZO_FAIL("AvgElectricField is null!");

      for(int i=0;i<ElectricSize[field];i++){
        AvgElectricField[field][i] += ElectricField[field][i];
      }
    }//field
  }//level>0

  return SUCCESS;
}

