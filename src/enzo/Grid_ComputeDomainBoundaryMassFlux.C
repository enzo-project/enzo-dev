/**********************************************************************
/
/  GRID CLASS (COMPUTE DOMAIN BOUNDARY MASS FLUX)
/
/  written by: Andrew Emerick
/  date:       January, 2017
/  modified1:
/
/  PURPOSE:
/    For density field and all species fields, compute
/    mass flux on domain boundary to count mass that
/    has left the domain. Stores this in global counter
/
/  RETURNS:
/     SUCCESS or FAIL
***********************************************************************/

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

int grid::ComputeDomainBoundaryMassFlux(float *allgrid_BoundaryMassFluxContainer)
{

  /* Return if this doesn't concern us */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  PrepareBoundaryMassFluxFieldNumbers(); //  this should be put somewhere else

  /* Is this grid on the edge of the domain */
  int GridOffsetLeft[MAX_DIMENSION], GridOffsetRight[MAX_DIMENSION];
  for (int dim = 0; dim < MAX_DIMENSION; dim++){
    if (dim < this->GridRank){
      GridOffsetLeft[dim] = nint((this->GridLeftEdge[dim] - DomainLeftEdge[dim])/
                                  this->CellWidth[dim][0]);
      GridOffsetRight[dim] = nint((this->GridRightEdge[dim] - DomainRightEdge[dim])/
                                  this->CellWidth[dim][0]);
    } else{
      GridOffsetLeft[dim] = 0; GridOffsetRight[dim] = 0;
    }
  }

  int NumberOfBoundaryMassFields   = 1; // density always

  if (MultiSpecies > 0)
    NumberOfBoundaryMassFields += 5; // HI, HII, HeI, HeII, HeIII

  if (MultiSpecies > 1)
    NumberOfBoundaryMassFields += 3; // H2, H2I, HM

  if (MultiSpecies > 2)
    NumberOfBoundaryMassFields += 3; // D, D+, HD

  if (MultiMetals > 0)
    NumberOfBoundaryMassFields += 1; // metallicity

  if (MultiMetals > 1 && STARMAKE_METHOD(INDIVIDUAL_STAR)){
    NumberOfBoundaryMassFields += (StellarYieldsNumberOfSpecies);
    if (MultiSpecies > 0) NumberOfBoundaryMassFields -= 2; // don't double count H and He from MultiSpecies
  }

  for (int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i ++){
    grid_BoundaryMassFluxContainer[i] = 0.0;
  }

  /* conversion factor to go from flux to actual mass */
  float conversion;

  /* now loop over each dimension to check if at boundary - sum if we are */
  for(int dim = 0; dim < MAX_DIMENSION; dim++){

    int size = 1;
    for (int j = 0; j < GridRank; j++)
      size *= BoundaryFluxes->LeftFluxEndGlobalIndex[dim][j] -
                BoundaryFluxes->LeftFluxStartGlobalIndex[dim][j] + 1;

    /* dx**2 * dt * dx**3  (convert from dt * flux (density / area)  to mass) */
    conversion = POW( this->CellWidth[dim][0] , 5);

    if (GridOffsetLeft[dim] != 0 && GridOffsetRight[dim] != 0) continue; // not at any domain boundary

    /* if here, we are at a domain boundary */
    for (int i = 0; i < NumberOfBoundaryMassFields; i++){
      float left_mass = 0.0, right_mass = 0.0;

      int field_num = BoundaryMassFluxFieldNumbers[i];
      for (int index = 0; index < size; index ++){
        if (GridOffsetLeft[dim] == 0)
          left_mass = (-1.0 * this->BoundaryFluxes->LeftFluxes[field_num][dim][index], 0.0);

        if (GridOffsetRight[dim] == 0)
          right_mass += max(this->BoundaryFluxes->RightFluxes[field_num][dim][index], 0.0);
      }

      // only keep outflow (left < 0 = outflow, multiply by -1)
      grid_BoundaryMassFluxContainer[i] += (left_mass + right_mass)*conversion ; //max( -1.0*left_mass, 0.0) + max(right_mass, 0.0);

      // add into global container on this processor
      allgrid_BoundaryMassFluxContainer[i] += grid_BoundaryMassFluxContainer[i];

      if (grid_BoundaryMassFluxContainer[i] > 0) printf("field_num = %"ISYM" mass = %"ESYM"\n", field_num, grid_BoundaryMassFluxContainer[i]);

    } // end loop over fields
  } // end loop over dim

  return SUCCESS;
}
