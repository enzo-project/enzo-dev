/***********************************************************************
/
/  GRID CLASS (CREATE / CLEAN UP TEMPERATURE FIELD FOR H2 SHIELDING)
/
/  written by: John Regan (using Ji-hoon Kim template)
/  date:       April 2014
/  modified1:
/
/  PURPOSE: create or clean up temperature field for the use 
/           in H2 Shielding in Grid_WalkPhotonPackage.C 
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

int grid::InitializeTemperatureFieldForH2Shield() 
{
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (RadiativeTransferH2ShieldType != 1)
    return SUCCESS;

  int TemperatureField, size = 1;
  float *temperature;

  for (int dim=0; dim<GridRank; dim++) size *= GridDimension[dim];
  temperature = new float[size];
  
  if (this->ComputeTemperatureField(temperature) == FAIL) {
    fprintf(stderr, "Error in grid->ComputeTemperatureField.\n");
    return FAIL;
  }

  TemperatureField = this->GetTemperatureFieldNumberForH2Shield();
  
  BaryonField[TemperatureField] = new float[size];
  for (int i = 0; i < size; i++)
    BaryonField[TemperatureField][i] = temperature[i];

  delete [] temperature;

  return SUCCESS;
}


int grid::FinalizeTemperatureFieldForH2Shield() 
{
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (RadiativeTransferH2ShieldType != 1)
    return SUCCESS;

  int TemperatureField;
  TemperatureField = this->GetTemperatureFieldNumberForH2Shield();

  if (BaryonField[TemperatureField] != NULL) {
    delete [] BaryonField[TemperatureField];
    BaryonField[TemperatureField] = NULL;
  }

  return SUCCESS;
}


int grid::GetTemperatureFieldNumberForH2Shield()
{
  return NumberOfBaryonFields+1;
}

