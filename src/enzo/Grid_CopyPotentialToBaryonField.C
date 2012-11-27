/***********************************************************************
/
/  GRID CLASS (COPY GRAVITATIONAL POTENTIAL FIELD TO BARYON FIELD)
/
/  written by: Alexei Kritsuk
/  date:       Aug 2001
/  modified1:  Robert Harkness
/  date:       June 2004
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
 
// Copy the potential field to baryon field
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int FindField(int field, int farray[], int numfields);
 
int grid::CopyPotentialToBaryonField()
{

  if (CopyGravPotential == FALSE)
    return SUCCESS;

  // Return if this doesn't concern us
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  // Find the field
 
  int field = FindField(GravPotential, FieldType, NumberOfBaryonFields);
 
  // Check to make sure BaryonField "GravPotential" exists
 
  if (BaryonField[field] == NULL) {
    ENZO_FAIL("GravPotential field missing.\n");
  }
 
  // Check to make sure PotentialField exists
 
  if (PotentialField == NULL) {
    ENZO_FAIL("PotentialField missing.\n");
  }
 
  /* Well, it appears that currently GravitatingMassField is larger
     than active BaryonField by 2*max(BufferSize, NumberOfGhostZones) = 12
     zones. BufferSize = GRAVITY_BUFFER_SIZE*RefinementFactor = 6.
     In other words, GravitatingMassField has 6 ghost zones, compared to 3
     for a BaryonField. That is why we use shift 6.
     See Grid_InitializeGravitatingMassField.C for details.
  */
 
  int BaryonFieldBufferSize = NumberOfGhostZones;
  int GravityBufferSize = GRAVITY_BUFFER_SIZE;
  int DimTemp, BufferSize;
  int dim;
  int Off[3];
 
  for (dim = 0; dim < GridRank; dim++)
  {
     DimTemp = GridEndIndex[dim] - GridStartIndex[dim] + 1;
     BufferSize = (GravitatingMassFieldDimension[dim] - DimTemp)/2;
     Off[dim] = (GravitatingMassFieldDimension[dim] - GridDimension[dim])/2;
//     fprintf(stderr, "CPOT (%"ISYM") %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n", dim, GridDimension[dim], GridStartIndex[dim], GridEndIndex[dim], GravitatingMassFieldDimension[dim], BufferSize, Off[dim]);
  }
 
  int i, j, k;
  int index;
  int jj = 0;
 
  float maxPot=-1e30, minPot=1e30;
  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++)
    {
      index = (((k+Off[2])*GravitatingMassFieldDimension[1]) + 
	       (j+Off[1]))*GravitatingMassFieldDimension[0] + Off[0];
 
      for (i = 0; i < GridDimension[0]; i++, index++)
      {
	//	BaryonField[field][jj++] = GravitatingMassField[index] + 1; // use this for debugging 
	BaryonField[field][jj++] = PotentialField[index];
	// debuggin:
	maxPot = max(maxPot,PotentialField[index]);
	minPot = min(minPot,PotentialField[index]);
      }
 
    }
  }
 
//  fprintf(stderr, "STUFF field %"ISYM"  elements %"ISYM"  %16.8e  %16.8e\n", field, crap, big, low);
//  if (debug1) printf("Potential minimum: %g \t maximum: %g\n", minPot, maxPot);


  return SUCCESS;
}
