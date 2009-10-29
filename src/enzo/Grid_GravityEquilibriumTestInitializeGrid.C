/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GRAVITY EQUILIBRIUM TEST)
/
/  written by: Greg Bryan
/  date:       June, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
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
 
/* function prototypes */
 
int grid::GravityEquilibriumTestInitializeGrid(
				   float ScaleHeight)
{
  /* declarations */
 
  int dim, field, i, index, j, k;
  float density, density_old, pressure;
 
  /* error check */
 
  /* create fields */
 
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  int vel = NumberOfBaryonFields - 1;
  if (GridRank > 1)
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* Determine the size of the fields. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Allocate space for the fields. */
 
  for (field = 0; field < NumberOfBaryonFields; field++)
    BaryonField[field] = new float[size];
 
  /* set density, total energy and velocity in problem dimension */
 
  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
 
    /* Set density for exponential atmosphere. */
 
    density = 1.0*PEXP(CellLeftEdge[0][i]/ScaleHeight);
    if (i == GridStartIndex[0])
      density_old = density;
 
    /* Set pressure to 1 at bottom of atmosphere, otherwise compute hydrostatic
       equilibrium. */
 
    if (i == GridStartIndex[0])
      pressure = 1.0;
    else
      pressure = pressure + 0.5*UniformGravityConstant*
	                    (CellWidth[0][i-1]*density_old +
			     CellWidth[0][i  ]*density    );
    density_old = density;
 
//    printf("%"ISYM" %"GSYM" %"GSYM"\n", i, density, pressure);
 
    /* Loop over this level, set density, energy. */
 
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	index = (k*GridDimension[1] + j)*GridDimension[0] + i;
	BaryonField[0][index] = density;
	BaryonField[1][index] = pressure/(density*(Gamma-1.0));
	if (DualEnergyFormalism)
	  BaryonField[2][index] = pressure/(density*(Gamma-1.0));
      }
 
  } // end loop over i slices
 
  /* Set velocities to zero. */
 
  for (field = vel; field < vel+GridRank; field++)
    for (i = 0; i < size; i++)
      BaryonField[field][i] = 0.0;
 
  return SUCCESS;
}
