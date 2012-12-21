
#include <stdio.h>
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

  for(int field=0;field<3;field++){
    
    if( AvgElectricField[field] != NULL )
      delete AvgElectricField[field];
    
    AvgElectricField[field] = new float[ ElectricSize[field] ];
    
    for(int i=0; i<ElectricSize[field]; i++)
      AvgElectricField[field][i] = 0.0;
  }

  return SUCCESS;
}
