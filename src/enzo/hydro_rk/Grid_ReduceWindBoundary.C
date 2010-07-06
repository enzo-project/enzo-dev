/***********************************************************************
/
/  GRID CLASS (UPDATE MHD VARIABLES)
/
/  written by: Peng Wang
/  date:       August, 2008
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::ReduceWindBoundary()
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, TimeUnits = 1.0,
    VelocityUnits = 1.0;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

  int igrid;
  float Vmax = 5e5/VelocityUnits;
  float WindReduceFactor = 10.0;

  /* Left x boundary */
  
  if (GridLeftEdge[0] == DomainLeftEdge[0]) {
    
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridStartIndex[0]; i++) {
	  igrid = i + (j + k*GridDimension[1])*GridDimension[0];
	  if (BaryonField[Vel1Num][igrid] > Vmax) {
	    BaryonField[Vel1Num][igrid] /= WindReduceFactor;
	  }
	}
      }
    }

  }

  /* Right x boundary */
  
  if (GridRightEdge[0] == DomainRightEdge[0]) {
    
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = GridEndIndex[0]+1; i < GridDimension[0]; i++) {
	  igrid = i + (j + k*GridDimension[1])*GridDimension[0];
	  if (BaryonField[Vel1Num][igrid] < -Vmax) {
	    BaryonField[Vel1Num][igrid] /= WindReduceFactor;
	  }
	}
      }
    }

  }

  /* Left y boundary */
  
  if (GridLeftEdge[1] == DomainLeftEdge[1]) {
    
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridStartIndex[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++) {
	  igrid = i + (j + k*GridDimension[1])*GridDimension[0];
	  if (BaryonField[Vel2Num][igrid] > Vmax) {
	    BaryonField[Vel2Num][igrid] /= WindReduceFactor;
	  }
	}
      }
    }

  }

  /* Right y boundary */
  
  if (GridRightEdge[1] == DomainRightEdge[1]) {
    
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = GridEndIndex[1]+1; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++) {
	  igrid = i + (j + k*GridDimension[1])*GridDimension[0];
	  if (BaryonField[Vel2Num][igrid] < -Vmax) {
	    BaryonField[Vel2Num][igrid] /= WindReduceFactor;
	  }
	}
      }
    }

  }

  /* Left z boundary */
  
  if (GridLeftEdge[2] == DomainLeftEdge[2]) {
    
    for (int k = 0; k < GridStartIndex[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++) {
	  igrid = i + (j + k*GridDimension[1])*GridDimension[0];
	  if (BaryonField[Vel3Num][igrid] > Vmax) {
	    BaryonField[Vel3Num][igrid] /= WindReduceFactor;
	  }
	}
      }
    }

  }

  /* Right z boundary */
  
  if (GridRightEdge[2] == DomainRightEdge[2]) {
    
    for (int k = GridEndIndex[2]+1; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++) {
	  igrid = i + (j + k*GridDimension[1])*GridDimension[0];
	  if (BaryonField[Vel3Num][igrid] < -Vmax) {
	    BaryonField[Vel3Num][igrid] /= WindReduceFactor;
	  }
	}
      }
    }

  }
  
  return SUCCESS;
}
