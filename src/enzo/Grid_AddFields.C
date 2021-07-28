/***********************************************************************
/
/  GRID CLASS (ADD RADIATION FIELDS WHEN RESTARTING FROM NON-RT DATA)
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "StarParticleData.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

int grid::AddFields(int TypesToAdd[], int NumberOfFields)
{

  int dim, i, j, n, size = 1;
  float value;

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (i = 0; i < NumberOfFields; i++) {

    n = NumberOfBaryonFields + i;
    FieldType[n] = TypesToAdd[i];

    if (ProcessorNumber == MyProcessorNumber) {
      if (BaryonField[n] != NULL) {
	ENZO_VFAIL("BaryonField[%"ISYM"] already assigned?\n", n)

      }
      BaryonField[n] = new float[size];

      // added conditional for using a metallicity floor with rad-hydro
      if (TypesToAdd[i] == Metallicity ){
        value = MechStarsMetallicityFloor*0.02; // parameter is in zsun
          for (j = 0; j < size; j++)  
            BaryonField[n][j] = value * BaryonField[0][j];        
      }
      else{
        value = (TypesToAdd[i] == SNColour ||
	       TypesToAdd[i] == MetalSNIaDensity || TypesToAdd[i]==MetalSNIIDensity) ? tiny_number : 0.0;
      
      for (j = 0; j < size; j++)
	      BaryonField[n][j] = value;
          }
      } // ENDIF this processor

  } // ENDFOR i

  NumberOfBaryonFields += NumberOfFields;

  return SUCCESS;

}
