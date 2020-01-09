/***********************************************************************
/  GRID CLASS (INITIALIZE THE GRID FOR THE SINE WAVE TEST)
/
/  written by: JC Passy
/  date:       June, 2013
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "phys_constants.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::TestGravitySineWaveInitializeGrid(float Amplitude,
                                            float Period,
                                            float Angle,
                                            int PoissonSineWaveSubgridsAreStatic,
                                            int TotalRefinement,
                                            int grid_num)
{
  /* declarations */

  int dim, i, j, k, size, field, vel;

  /* Create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1)
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  if (WritePotential)
    FieldType[NumberOfBaryonFields++] = GravPotential;
  
  /* Set the subgrid static flag. */

  SubgridsAreStatic = PoissonSineWaveSubgridsAreStatic;

  /* Return if this doesn't concern us. */
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  
  /* compute size of fields */
  size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
    
  /* allocate fields */
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];

  /* Set velocities to 0 BaryonField[1,2,3] */
  for (dim = 0; dim < GridRank; dim++)
    for (i = 0; i < size; i++)
      BaryonField[vel+dim][i] = 0.0;
      
  int index;
  float x,y,z;
  float xtilted,dx,xtmp,mean_rho;

  int m;

  for (k = 0; k < GridDimension[2]; k++) 
    for (j = 0; j < GridDimension[1]; j++) 
      for (i = 0; i < GridDimension[0]; i++) {
          
	index = i + GridDimension[0]*(j + GridDimension[1]*k);

	/* Compute position */
          
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	
	/* Density  BaryonField[0]*/
	
	xtilted = x*cos(Angle*pi/180.0) + y*sin(Angle*pi/180.0);
	dx = CellWidth[0][i];
	
	/* Compute mean rho (volume-averaged) only for Angle = 0.0 */
	mean_rho = 0.0;
	if (Angle > tiny_number)
	  mean_rho = Amplitude*(2.0 +sin(2*pi*xtilted/Period));
	else {
	  for (m = -30; m < 30; m++) {
	    xtmp = x + m/60*CellWidth[0][i];
	    mean_rho = mean_rho + Amplitude*(2.0 +sin(2*pi*xtmp/Period)) ;
	  }
	  mean_rho = mean_rho/60.0;
	}
	
	BaryonField[0][index] = mean_rho;
		
	/* Total SPECIFIC Energy BaryonField[1], don't forget the macroscopic kinetic energic */
	BaryonField[1][index] = 1.0;
      }
  
  /* Set internal energy if necessary. */
  /* Since Velocities are set to 0, this is same as total energy */
  if (DualEnergyFormalism)  
    for (i = 0; i < size; i++) 
      BaryonField[2][i] = BaryonField[1][i];

  return SUCCESS;
}
