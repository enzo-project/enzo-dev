/***********************************************************************
/
/  GRID CLASS (SET SHEARING BOX BOUNDARY CONDITION)
/
/  written by: Fen Zhao
/  date:       2008
/  modified1: Peng Wang
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::SetShearingBoxExternalBoundaries()
{

  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
   
  /* Only in the X direction, set the shearing boundaries
     assumes the cellwidths are uniform in each direction; which is an okay assumption */    

  FLOAT Lx = DomainRightEdge[0] - DomainLeftEdge[0],
    Ly = DomainRightEdge[1] - DomainLeftEdge[1];

  int newj1, newj2; // the index of the two sheared cells

  /* We only need to calculate the shift index, j_shift, and interpolation 
   * coefficent, a, once for every time.
   * Interpolation coefficient: f_shear = a * f(newj1) + (1 - a) * f(newj2) 
   */

  FLOAT dis = AngularVelocity*VelocityGradient*Lx*Time;
  FLOAT dis_shear = dis - int(dis/Ly)*Ly;
  int j_shift = int(dis_shear/CellWidth[1][GridStartIndex[1]]);

  FLOAT newy = CellLeftEdge[1][GridStartIndex[1]+j_shift] + CellWidth[1][GridStartIndex[1]+j_shift];    
  float a = (newy - CellLeftEdge[1][GridStartIndex[1]] - dis_shear) / CellWidth[1][GridStartIndex[1]+j_shift];

  int Yactivesize = GridDimension[1] - 2*DEFAULT_GHOST_ZONES;

  /* set Shearing box for right edge */
    
  for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {

    newj1 = j + j_shift;
    if (newj1 > GridEndIndex[1])
      newj1 -= Yactivesize;
    
    if (newj1 == GridEndIndex[1])
      newj2 = GridStartIndex[1];
    else 
      newj2 = newj1 + 1;
    
    for (int field = 0; field < NumberOfBaryonFields; field++) {
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int shift = 0; shift < DEFAULT_GHOST_ZONES; shift++) {
	  float val1 = BaryonField[field][GetIndex(GridStartIndex[0]+shift, newj1, k)];
	  float val2 = BaryonField[field][GetIndex(GridStartIndex[0]+shift, newj2, k)];
	  
	  BaryonField[field][GetIndex(GridEndIndex[0]+1+shift, j, k)] = val1 * a + (1.0 - a) * val2;	  
	}
      }
    }

  } // set left edge

  /* set Shearing box for left edge */
  
  for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {

    newj1 = j - j_shift;
    if (newj1 < GridStartIndex[1])
      newj1 += Yactivesize;
    
    if (newj1 == GridEndIndex[1])
      newj2 = GridStartIndex[1];
    else 
      newj2 = newj1 + 1;
    
    for (int field = 0; field < NumberOfBaryonFields; field++) {
      for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	for (int shift = 0; shift < DEFAULT_GHOST_ZONES; shift++) {
	  float val1 = BaryonField[field][GetIndex(GridEndIndex[0]-shift, newj1, k)];
	  float val2 = BaryonField[field][GetIndex(GridEndIndex[0]-shift, newj2, k)];
	  
	  BaryonField[field][GetIndex(GridStartIndex[0]-1-shift, j, k)] = val1 * a + (1.0 - a) * val2;	  
	}
      }
    }

  } // set right edge

  /* Add shearing velocity to vy and correct for total Energy */

  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int shift = 0; shift < DEFAULT_GHOST_ZONES; shift++){

	/* x left edge */

	int igrid = GetIndex(GridStartIndex[0]-1-shift, j, k);

	float rho = BaryonField[iden][igrid];
	float vx  = BaryonField[ivx][igrid];  
	float vy  = BaryonField[ivy][igrid]; 
	float vz  = BaryonField[ivz][igrid];
	float v2  = vx*vx + vy*vy + vz*vz;
	
	BaryonField[ietot][igrid] -= 0.5*v2;
	BaryonField[ivy][igrid] += AngularVelocity*VelocityGradient*Lx;	  
	v2 = vx*vx + pow(BaryonField[ivy][igrid],2) + vz*vz;
	BaryonField[ietot][igrid] += 0.5*v2;
	
	/* x right edge */
	
	igrid = GetIndex(GridEndIndex[0]+1+shift, j, k);
	
	rho = BaryonField[iden][igrid];
	vx  = BaryonField[ivx][igrid];  
	vy  = BaryonField[ivy][igrid]; 
	vz  = BaryonField[ivz][igrid];
	v2  = vx*vx + vy*vy + vz*vz;
	
	BaryonField[ietot][igrid] -= 0.5*v2;
	BaryonField[ivy][igrid] -= AngularVelocity*VelocityGradient*Lx;
	v2 = vx*vx + pow(BaryonField[ivy][igrid],2) + vz*vz;
	BaryonField[ietot][igrid] += 0.5*v2;
    
      }
    }
  }

  /*for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
    int igrid = GridStartIndex[0]-1 + j*GridDimension[0];
    printf("%g ", BaryonField[ivy][igrid]);
  }
  printf("\n");*/

  /*printf("vy:\n");
  for (int i = 0; i < GridDimension[0]; i++) {
    int igrid = i + GridStartIndex[1]*GridDimension[0];
    printf("%g ", BaryonField[ivy][igrid]);
  }
  printf("\n");*/

  return SUCCESS;
}

