#ifdef SAB

//  grid::AttachAcceleration  and grid::DetachAcceleration().
// Pointer juggling for the boundary set of the acceleration field.

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

#ifdef FAST_SIB
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#else
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          int level, TopGridData *MetaData,
                          ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#endif




// Begin the pointer juggle to set the boundary on the acceleration field.
// Save all BaryonField pointers in temporary array, and set them to be Acceleration Field
// pointers.  This lets the SetBoundary condition machenery operate without heft code rewrite.


int grid::AttachAcceleration(){


  //This redundancy check is for the parent grid.  Multiple subgrids will have the same 
  //parent grid.

  int field;

  if( AccelerationHack == TRUE )
    return SUCCESS;
  else
    AccelerationHack = TRUE;


  ActualNumberOfBaryonFields = NumberOfBaryonFields;
  NumberOfBaryonFields = GridRank; 

  for(field = 0; field < ActualNumberOfBaryonFields; field++){

    ActualBaryonField[field] = BaryonField[field];
    ActualOldBaryonField[field] = OldBaryonField[field];
    ActualFieldType[field] = FieldType[field];


    if( field < GridRank ){

      BaryonField[field] = AccelerationField[field];
      OldBaryonField[field] = OldAccelerationField[field];

    }else{

      BaryonField[field] = NULL;
      OldBaryonField[field] = NULL;

    }

    FieldType[field]=FieldUndefined;

  }

  FieldType[0] = ((GridRank >= 1 ) ? Acceleration0 : FieldUndefined );
  FieldType[1] = ((GridRank >= 2 ) ? Acceleration1 : FieldUndefined );
  FieldType[2] = ((GridRank >= 3 ) ? Acceleration2 : FieldUndefined );

 
  
  return SUCCESS;
}


// end pointer juggle for Boundary Set of AccelerationField.
// Return saved BaryonField pointers to their rightful position.

int grid::DetachAcceleration(){

  int field;

  if( AccelerationHack == FALSE )
    return SUCCESS;  // the detachment has already been done.
  else
    AccelerationHack = FALSE;
    
  NumberOfBaryonFields = ActualNumberOfBaryonFields;

  for( field = 0; field < NumberOfBaryonFields; field++){
    
    BaryonField[field] = ActualBaryonField[field];
    OldBaryonField[field] = ActualOldBaryonField[field];
    FieldType[field] = ActualFieldType[field];

  }


  return SUCCESS;
}

//SetAccelerationBoundary ensures that all subgrids agree in the boundary.
//Not a big deal for hydro, but essential for DivB = 0 in MHD runs.
//Only called on level > 0 because the root grid is dealt with differently than SG's.

int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
			    SiblingGridList SiblingList[],
			    int level, TopGridData *MetaData,
			    ExternalBoundary *Exterior, LevelHierarchyEntry * Level,
			    int CycleNumber)
{

  if ( ! ( (SelfGravity || UniformGravity || PointSourceGravity) && level > 0 ))
    return SUCCESS;

  if ( Grids[0]->GridData->ReturnNumberOfBaryonFields() == 0 ){
      return SUCCESS;
  }

  //Set the boundary on the Acceleration field.  Reuse SetBoundaryConditions.  
  //Juggle pointers around.

  int grid1, ConservativeTruth;
  char basename[30];  

  //We don't want conservative interpolation actually being done for the acceleration field.
  ConservativeTruth = ConservativeInterpolation;
  ConservativeInterpolation = FALSE;

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
    if( Grids[grid1]->GridData->AttachAcceleration() == FAIL ) {
      ENZO_FAIL("Error in AttachAcceleration \n");
    }
    if( Grids[grid1]->ParentGrid->GridData->AttachAcceleration() ==FAIL ){
      ENZO_FAIL("Error in AttachAcceleration, Parent \n");
    }
  }

#ifdef FAST_SIB
  if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList, level, MetaData,
			    NULL, NULL) == FAIL)
    ENZO_FAIL("SetBoundaryConditions() failed!\n");
#else
  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData, 
			    NULL, NULL) == FAIL)
    ENZO_FAIL("SetBoundaryConditions() failed!\n");
#endif
  

  
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

    if( Grids[grid1]->GridData->DetachAcceleration() == FAIL ) {
      ENZO_FAIL("Error in DetachAcceleration\n");
    }
    if( Grids[grid1]->ParentGrid->GridData->DetachAcceleration() == FAIL ) {
      ENZO_FAIL("Error in DetachAcceleration, parent\n");

    }

  }


  ConservativeInterpolation = ConservativeTruth;

  return SUCCESS;

}

#endif /* SAB */
