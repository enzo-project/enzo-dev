/***********************************************************************
/
/  GRID CLASS (Setup MHDCT Dimensions)
/
/  written by: Sam Skillman
/  date:       September, 2012
/  modified1: 
/
/  PURPOSE: Helper function to set up MHDCT dimensions.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void grid::MHD_SetupDims(void){
  if (!UseMHDCT)
    return;

  for(int field=0; field<3; field++){
    this->MagneticSize[field] = 1;
    this->ElectricSize[field] = 1;

    for(int dim=0; dim<3; dim++){
      this->MagneticDims[field][dim] = this->GridDimension[dim];
      this->ElectricDims[field][dim] = this->GridDimension[dim] +1;


      this->MHDStartIndex[field][dim] = this->GridStartIndex[dim];
      this->MHDEndIndex[field][dim] = this->GridEndIndex[dim];

      this->MHDeStartIndex[field][dim] = this->GridStartIndex[dim];
      this->MHDeEndIndex[field][dim] = this->GridEndIndex[dim]+1;

      this->MHDAdd[field][dim]=0;
      if( field == dim ){
        this->MagneticDims[field][dim]++;
        this->ElectricDims[field][dim]--;
        this->MHDEndIndex[field][dim]++;
        this->MHDeEndIndex[field][dim]--;
        this->MHDAdd[field][dim]=1;
      }

      this->MagneticSize[field] *= this->MagneticDims[field][dim];
      this->ElectricSize[field] *= this->ElectricDims[field][dim];
    }

  }
  return;
}
