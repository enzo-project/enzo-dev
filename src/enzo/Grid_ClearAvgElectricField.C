/***********************************************************************
/
/  GRID CLASS (ZERO AND ALLOCATE AvgElectricField)
/
/  written by: David Collins
/  date:       2005
/  modified1:
/
/  PURPOSE:  This is an analogous routine to the grid::ClearBoundaryFluxes 
/            routine, but acts on the AvgElectricField.  AvgElectricField
/            acts similarly to the BoundaryFluxes, but in order to more 
/            robustly preserve Div(B) this quantity is used for the projection
/            of the electric field.  See Collins et al 2010 for a description.
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

int grid::ClearAvgElectricField(){
  
  //Also clear and allocate the Average Electric Field.
  //This has a similar function the Fluxes, but not the limited spatial extent that 
  //the Fluxes object does.

  if(MyProcessorNumber != ProcessorNumber )
    return SUCCESS;
  if ( UseMHDCT != TRUE )
      return SUCCESS;

  for(int field=0;field<3;field++){
    
    if( AvgElectricField[field] != NULL )
      delete AvgElectricField[field];
    
    AvgElectricField[field] = new float[ ElectricSize[field] ];
    
    for(int i=0; i<ElectricSize[field]; i++)
      AvgElectricField[field][i] = 0.0;
  }

  return SUCCESS;
}
