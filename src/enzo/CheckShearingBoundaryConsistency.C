/***********************************************************************
/
/  CHECK INPUT SHEARING BOUNDARIES ARE CONSISTENT
/
/  written by: Fen Zhao
/  date:       June, 2009
/  modified1:
/
/  PURPOSE:
/         Check that input parameters for Shearing Boundaries are sensical 
/         and sets helper variables properlyk 
/
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"

int CheckShearingBoundaryConsistency(TopGridData &MetaData){

  //Check only one direction has a shearing boundary condition and both faces shearing
  int numberOfShearingBoundariesL = 0 ,  numberOfShearingBoundariesR = 0 ;
  for (int i=0; i<MetaData.TopGridRank; i++){
    if (MetaData.LeftFaceBoundaryCondition[i]== shearing) { numberOfShearingBoundariesL++; ShearingBoundaryDirection=i;}
    if (MetaData.RightFaceBoundaryCondition[i]== shearing) { numberOfShearingBoundariesR++; ShearingBoundaryDirection=i;}
  }

  if (numberOfShearingBoundariesL> 1 || numberOfShearingBoundariesR > 1) 
    ENZO_FAIL("Too Many Shearing Boundaries");
  if (numberOfShearingBoundariesL == 0 && numberOfShearingBoundariesR == 0) 
    return SUCCESS;
 

  //Check that we have one direction at least where there are periodic BC on both sides, fail otherwise
  int numberOfPeriodicBoundariesBoth = 0;
  int tempShearingVelocityDirection;
  for (int i=0; i<MetaData.TopGridRank; i++)
      if (MetaData.LeftFaceBoundaryCondition[i] == periodic && MetaData.RightFaceBoundaryCondition[i] == periodic){ 
      numberOfPeriodicBoundariesBoth++;
      tempShearingVelocityDirection=i;
    }
  

  if (numberOfPeriodicBoundariesBoth == 0) {
   
     ENZO_FAIL("Need a Periodic Boundary (on both faces) for a Shearing Boundary");
  }
   //if only one direction is periodic, then we autoset ShearingVelocityDirection to the right value
   if (numberOfPeriodicBoundariesBoth == 1){
     ShearingVelocityDirection=tempShearingVelocityDirection;
     return SUCCESS;
   }
   
   //if both directions periodic, then either we use the user specific ShearingVelocityDirection
   //or use the high index of the two directions
  
  if (numberOfPeriodicBoundariesBoth == 2){
     
     if (ShearingVelocityDirection==-1) ShearingVelocityDirection=tempShearingVelocityDirection;
     else if (MetaData.LeftFaceBoundaryCondition[ShearingVelocityDirection]!=periodic 
	      || MetaData.RightFaceBoundaryCondition[ShearingVelocityDirection]!=periodic)
       ENZO_FAIL("ShearingVelocityDirection is not Periodic on Both Faces");
   }
  
 

  if (ShearingBoundaryDirection==0){
    if (ShearingVelocityDirection==1) ShearingOtherDirection=2;
    else if (ShearingVelocityDirection==2) ShearingOtherDirection=1;}
  else if (ShearingBoundaryDirection==1){
    if (ShearingVelocityDirection==0) ShearingOtherDirection=2;
    else if (ShearingVelocityDirection==2) ShearingOtherDirection=0;}
  else if (ShearingBoundaryDirection==2){
    if (ShearingVelocityDirection==0) ShearingOtherDirection=1;
    else if (ShearingVelocityDirection==1) ShearingOtherDirection=0;} 


  return SUCCESS;
  
}
