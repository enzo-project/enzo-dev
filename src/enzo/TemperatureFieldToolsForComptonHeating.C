/***********************************************************************
/
/  GRID CLASS (CREATE / CLEAN UP TEMPERATURE FIELD FOR COMPTON HEATING)
/
/  written by: Ji-hoon Kim
/  date:       December, 2009
/  modified1:
/
/  PURPOSE: create or clean up temperature field for the use 
/           in Compton heating in Grid_WalkPhotonPackage.C 
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

int grid::InitializeTemperatureFieldForComptonHeating() 
{
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (RadiationXRayComptonHeating == FALSE)
    return SUCCESS;

  int TemperatureField, size = 1;
  float *temperature;

  for (int dim=0; dim<GridRank; dim++) size *= GridDimension[dim];
  temperature = new float[size];
  
  if (this->ComputeTemperatureField(temperature) == FAIL) {
    fprintf(stderr, "Error in grid->ComputeTemperatureField.\n");
    return FAIL;
  }

  TemperatureField = this->GetTemperatureFieldNumberForComptonHeating();
  
  BaryonField[TemperatureField] = new float[size];
  for (int i = 0; i < size; i++)
    BaryonField[TemperatureField][i] = temperature[i];

  delete [] temperature;

  return SUCCESS;
}


int grid::FinalizeTemperatureFieldForComptonHeating() 
{
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (RadiationXRayComptonHeating == FALSE)
    return SUCCESS;

  int TemperatureField;
  TemperatureField = this->GetTemperatureFieldNumberForComptonHeating();

  if (BaryonField[TemperatureField] != NULL) {
    delete [] BaryonField[TemperatureField];
    BaryonField[TemperatureField] = NULL;
  }

  return SUCCESS;
}


int grid::GetTemperatureFieldNumberForComptonHeating()
{
  return NumberOfBaryonFields+1;
}

