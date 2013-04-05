/***********************************************************************
/
/  GRID CLASS ()
/
/  written by: David Collins
/  date:       2005?
/  modified1:
/
/  PURPOSE:    Averages the face centered magnetic field to the cell centered one.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::CenterMagneticField(int * Start, int * End){

    //Only if CT is on.
    if(  MHD_CT_Method == CT_None )
        return SUCCESS;

  //Setting default start & end indicies.  If there's a slicker way to do this, I'm all ears, but I think
  // default setting in C++ must be a static variable.

  if( Start == NULL ) Start = this->GridStartIndex;
  if( End   == NULL )   End = this->GridEndIndex;

  //If this processor (MyProcessorNumber) doesn't have the data belonging to this grid (ProcessorNumber) 
  //leave quietly.
  if( ProcessorNumber != MyProcessorNumber )
    return SUCCESS;

  //CenterMagneticField is called from Set Boundary Conditions.  Set Boundary Conditions is sometimes called from
  //SetAccelerationBoundary, which does some pointer juggling in order to apply the proper boundary conditions 
  //to the acceleration field.

  if( AccelerationHack == TRUE )
    return SUCCESS;

  //
  // Ensure that the centering method is appropriate:
  //  For non-CT runs, it needs to be OFF (because you evolve CenteredB
  //  For GridRank < 3, use directy copy on the flat dimensions, simple averaging on non-flat.

  int i,j,k,dim, indexC, indexB1, indexB2;
  int Offset[3]={1,MagneticDims[1][0], MagneticDims[2][1]*MagneticDims[2][0]};

  if( GridRank < 3 ) Offset[2] = 0;
  if( GridRank < 2 ) Offset[1] = 0;

  for( k=0; k<=GridDimension[2]-1; k++)
    for( j=0; j<=GridDimension[1]-1; j++)
      for( i=0; i<=GridDimension[0]-1; i++)
	for(dim=0;dim<3;dim++){
	  indexC = i + GridDimension[0]*(j + GridDimension[1]*k);
	  indexB1= i + MagneticDims[dim][0]*(j+MagneticDims[dim][1]*k);
	  indexB2= i + MagneticDims[dim][0]*(j+MagneticDims[dim][1]*k) + Offset[dim];
	  CenteredB[dim][indexC] =  0.5 *(MagneticField[dim][indexB1] + MagneticField[dim][indexB2]);
	  
	}
  return SUCCESS;
}

