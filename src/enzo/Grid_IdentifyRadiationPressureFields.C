/***********************************************************************
/
/  GRID CLASS (IDENTIFY THE RADIATION PRESSURE FIELDS)
/
/  written by: John H. Wise (adopted from Grid_IdentifySpeciesFields)
/  date:       December, 2005
/  modified1:
/
/  PURPOSE: Return indeces of the relevant radiation pressure fields
/
/  NOTE: 
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

/* function prototypes */

int FindField(int f, int farray[], int n);

int grid::IdentifyRadiationPressureFields(int &RPresNum1, int &RPresNum2,
					  int &RPresNum3)
{
  
  RPresNum1 = RPresNum2 = RPresNum3 = 0;

  if (RadiativeTransfer && RadiationPressure) {

    RPresNum1 = FindField(RadPressure0, FieldType, NumberOfBaryonFields); 
    RPresNum2 = FindField(RadPressure1, FieldType, NumberOfBaryonFields); 
    RPresNum3 = FindField(RadPressure2, FieldType, NumberOfBaryonFields); 

    if (RPresNum1 < 0 || RPresNum2 < 0 || RPresNum3 < 0) {
      ENZO_VFAIL("Could not identify a RadiationPressureField! RadPressure0 = %"ISYM", RadPressure1 = %"ISYM", RadPressure2 = %"ISYM"\n",RPresNum1, RPresNum2, RPresNum3)

    }
  }

  return SUCCESS;
}
