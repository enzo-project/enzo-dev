/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY SLOPE)
/
/  written by: Tom Abel & Fen Zhao
/  date:       July, 2008
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::FlagCellsToBeRefinedByResistiveLength()
{
   if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }


  float denu = 1.0, lenu = 1.0, velu = 1.0, tu, tempu;
    
  int igrid;
  
  float rho, eint, etot, vx, vy, vz, v2, p, h, cs, dpdrho, dpde, l_res, 
    curlBx, curlBy, curlBz, absB2, curlB2;
  int xmo, xpo, ymo, ypo, zmo, zpo;
  FLOAT x, y, z, r;

    int iBx, iBy, iBz;
    if (useMHD){
      iBx=FindField(Bfield1, FieldType, NumberOfBaryonFields);
      iBy=FindField(Bfield2, FieldType, NumberOfBaryonFields);
      if (GridRank==3) iBz=FindField(Bfield3, FieldType, NumberOfBaryonFields);
      
    }

      
  for (int k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {	
	igrid = (j + k*GridDimension[1])*GridDimension[0]+i;
	xpo = (j + k*GridDimension[1])*GridDimension[0]+(i+1);
	xmo = (j + k*GridDimension[1])*GridDimension[0]+(i-1);
	ypo = (j+1 + k*GridDimension[1])*GridDimension[0]+i;
	ymo = (j-1 + k*GridDimension[1])*GridDimension[0]+i;
	zpo = (j + (k+1)*GridDimension[1])*GridDimension[0]+i;
	zmo = (j + (k-1)*GridDimension[1])*GridDimension[0]+i;
	rho = BaryonField[DensNum][igrid];

	curlBx   = ((BaryonField[iBz  ][ypo] - BaryonField[iBz  ][ymo]) -
		    (BaryonField[iBy  ][zpo] - BaryonField[iBy  ][zmo]))/2.;	
	curlBy   =  ((BaryonField[iBz  ][xpo] - BaryonField[iBz  ][xmo]) -
		     (BaryonField[iBx  ][zpo] - BaryonField[iBx  ][zmo]))/2. ;
	curlBz   = ((BaryonField[iBy  ][zpo] - BaryonField[iBy  ][zmo]) -
		    (BaryonField[iBz  ][ypo] - BaryonField[iBz  ][ymo]))/2.  ;


	absB2 = BaryonField[iBx  ][igrid]*BaryonField[iBx  ][igrid] + 
	  BaryonField[iBy  ][igrid]*BaryonField[iBy  ][igrid] + 
	  BaryonField[iBz  ][igrid]*BaryonField[iBz  ][igrid]  ;
	curlB2 = curlBx*curlBx + curlBy*curlBy + curlBz*curlBz;

	l_res = sqrt(absB2)/max(sqrt(curlB2),tiny_number);


	int oldflag = FlaggingField[igrid];
	if (CellWidth[0][i]*RefineByResistiveLengthSafetyFactor > l_res) {	  
	  FlaggingField[igrid]++;
	}

      }
    }
  }

  /* Count number of flagged Cells. */

  int NumberOfFlaggedCells = 0;
  for (int i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }

  return NumberOfFlaggedCells;
 
}
