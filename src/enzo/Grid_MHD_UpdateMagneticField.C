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
/  Behavior: In order to keep DivB=0 in a cosmological simulation, 
/            the magnetic field that is advanced is 
/            B_{half comoving} = B_{comoving}*sqrt(a).
/            The rest of enzo stores B_{comoving} in MagneticField and BaryonField[BiNum].
/            Thus in this routine, the following takes place:
/
/            B_{half comoving}^{N} = B_{comoving}^{N}*sqrt{a^{N}}
/            B_{half}^{N+1} = B_{half}^N - Curl(E_{half comoving})
/            B_{comoving}^{N+1} = B_{half}^{N+1}/sqrt{a^{N+1}}

/  Inputs:  level; if level > 0, the Average Electric field is updated.
/                  THis is unnecessary on the root grid.
/           TimeIsBeforeSetLevelTimestep: The second call to this
/                  routine is after SetLevelTimestep, so
/                  grid->Time = t^{N+1}
/                  but the update still needs t^{N}, t^{N+1},
/                  so t^{N} = grid->Time - grid->dtFixed.
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

int grid::MHD_UpdateMagneticField(int level, LevelHierarchyEntry * NextLevel, 
                                  int TimeIsBeforeSetLevelTimestep){

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
  FLOAT aN = 1, dadtN=0, aNp1=1, dadtNp1=0, sqrt_aN, inv_sqrt_aN, inv_sqrt_aNp1;
  FLOAT tN, tNp1; //t^{N}, t^{N+1}
  if ( TimeIsBeforeSetLevelTimestep ){
      tN   = Time;
      tNp1 = Time+dtFixed;
  }else{
      tN   = Time-dtFixed;
      tNp1 = Time;
  }
  int i;
  if( ComovingCoordinates ){
      CosmologyComputeExpansionFactor(tN  , &aN, &dadtN);
      CosmologyComputeExpansionFactor(tNp1, &aNp1, &dadtNp1);
      sqrt_aN = sqrt(aN); 
      inv_sqrt_aNp1 =  1./sqrt(aNp1);
      inv_sqrt_aN = 1./sqrt_aN;
      //Since MagneticField = OldMagneticField - Curl,
      //we only need to update OldMagneticField to half-comoving
      for(field=0;field<3;field++){
          for( i=0; i<MagneticSize[field];i++){
              OldMagneticField[field][i] *= sqrt_aN;
          }
      }
  }

  
  //dtUsed is set to 1 here because the field ElectricField is NOT E, but E/dt. (This is to ensure
  // that the field is properly projected to parents over subgrid timesteps of potentially uneven lengths.)
  // I left in the option of dt != 1 in order to maintain the versitility of the code: might come out later.

  dtUsed = 1.0;

   int CurlStart[3] = {0,0,0}, CurlEnd[3] = {GridDimension[0]-1,GridDimension[1]-1,GridDimension[2]-1};
   MHD_Curl( CurlStart,CurlEnd, 2);
  if ( ComovingCoordinates ){
      for(field=0;field<3;field++){
          for( i=0; i<MagneticSize[field];i++){
              MagneticField[field][i] *= inv_sqrt_aNp1;
          }
      }
      for(field=0;field<3;field++){
          for( i=0; i<MagneticSize[field];i++){
              OldMagneticField[field][i] *= inv_sqrt_aN;
          }
      }
  }

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

