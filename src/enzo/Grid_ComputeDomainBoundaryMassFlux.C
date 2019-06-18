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
#include "phys_constants.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);


int grid::ComputeDomainBoundaryMassFlux(float *allgrid_BoundaryMassFluxContainer,
                                        TopGridData *MetaData)
{

  /* Return if this doesn't concern us */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  PrepareBoundaryMassFluxFieldNumbers(); //  this should be put somewhere else

  /* get units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass unit

  /* get colour fields */
  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum,
    MetalIaNum, MetalIINum;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
              MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");


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

  if (MetalNum != -1){
    NumberOfBoundaryMassFields++;   // metallicity
    if (MultiMetals || TestProblemData.MultiMetals) {
      NumberOfBoundaryMassFields += 2;  // ExtraType0 and ExtraType1
    }
  }

  if (MetalIaNum       != -1) NumberOfBoundaryMassFields++;
  if (MetalIINum       != -1) NumberOfBoundaryMassFields++;
  if (SNColourNum      != -1) NumberOfBoundaryMassFields++;
  if (MBHColourNum     != -1) NumberOfBoundaryMassFields++;
  if (Galaxy1ColourNum != -1) NumberOfBoundaryMassFields++;
  if (Galaxy2ColourNum != -1) NumberOfBoundaryMassFields++;

  for (int i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i ++){
    grid_BoundaryMassFluxContainer[i] = 0.0;
  }

  /* conversion factor to go from flux to actual mass */
  float conversion;

  /* now loop over each dimension to check if at boundary - sum if we are */
  for(int dim = 0; dim < MAX_DIMENSION; dim++){

    if (GridOffsetLeft[dim] != 0 && GridOffsetRight[dim] != 0) continue; // not at any domain boundary

    int size = 1;
    for (int j = 0; j < GridRank; j++)
      size *= BoundaryFluxes->LeftFluxEndGlobalIndex[dim][j] -
                BoundaryFluxes->LeftFluxStartGlobalIndex[dim][j] + 1;

    /* (convert from dt * flux (density / area)  to mass in solar masses) */
    conversion = POW( this->CellWidth[dim][0] , 3) * MassUnits / SolarMass;

    /* if here, we are at a domain boundary */
    if (MetaData->LeftFaceBoundaryCondition[dim] == outflow) {
      for (int i = 0; i < NumberOfBoundaryMassFields; i++){
        float left_mass = 0.0;

        int field_num = BoundaryMassFluxFieldNumbers[i];
        for (int index = 0; index < size; index ++){
          if (GridOffsetLeft[dim] == 0)
            left_mass += max(-1.0 * this->BoundaryFluxes->LeftFluxes[field_num][dim][index], 0.0);
        }

        grid_BoundaryMassFluxContainer[i] += (left_mass) * conversion;
        // add into global container on this processor
        allgrid_BoundaryMassFluxContainer[i] += grid_BoundaryMassFluxContainer[i];
      }

    } // if left

    if (MetaData->RightFaceBoundaryCondition[dim] == outflow) {
      for (int i = 0; i < NumberOfBoundaryMassFields; i++){
        float right_mass = 0.0;

        int field_num = BoundaryMassFluxFieldNumbers[i];
        for (int index = 0; index < size; index++){
          if (GridOffsetRight[dim] == 0)
            right_mass += max(this->BoundaryFluxes->RightFluxes[field_num][dim][index], 0.0);
        }

        grid_BoundaryMassFluxContainer[i] += (right_mass)*conversion;
        // add into global container on this processor
        allgrid_BoundaryMassFluxContainer[i] += grid_BoundaryMassFluxContainer[i];
      }
    } // if right

//      if (grid_BoundaryMassFluxContainer[i] > 0) printf("field_num = %"ISYM" mass = %"ESYM"\n", field_num, grid_BoundaryMassFluxContainer[i]);

  } // end loop over dim

  return SUCCESS;
}
