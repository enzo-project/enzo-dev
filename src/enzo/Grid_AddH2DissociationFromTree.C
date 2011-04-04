/***********************************************************************
/
/  ADD H2 DISSOCIATION EMISSION FROM STAR PARTICLES FROM A TREE
/
/  written by: John Wise
/  date:       March, 2011
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
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"
#include "Star.h"

int FindSuperSourceByPosition(FLOAT *pos, SuperSourceEntry **result,
			      int DEBUG);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::AddH2DissociationFromTree(void)
{

  int i, j, k, index, dim, ci;
  FLOAT pos[MAX_DIMENSION];
  FLOAT radius2;
  FLOAT innerFront, outerFront, innerFront2, outerFront2;
  double Luminosity[MAX_ENERGY_BINS];
  float energies[MAX_ENERGY_BINS], kdiss_r2;

  const double pc = 3.086e18, clight = 3e10;
  double H2Luminosity, H2ISigma = 3.71e-18;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Exit if there are no sources */

  if (SourceClusteringTree == NULL)
    return SUCCESS;

  this->DebugCheck((char*) "Grid_AddH2Dissociation");

  /* Get photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, kphHeIINum, 
				  kdissH2INum);

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, 
    TimeUnits, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Absorb the unit conversions into the cross-section
  H2ISigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);

  // Dilution factor (prevent breaking in the rate solver near the star)
  float dilutionRadius = 10.0 * pc / (double) LengthUnits;
  float dilRadius2 = dilutionRadius * dilutionRadius;

  // Convert from #/s to RT units
  double LConv = (double) TimeUnits / pow(LengthUnits,3);
  double LConv_inv = 1.0 / LConv;

  /* Find sources in the tree that contribute to the cells */

  SuperSourceEntry *Leaf, *LastLeaf;
  FLOAT dx;
  float factor = LConv_inv * H2ISigma / (4.0 * M_PI);
  int DEBUG;
  int dbg[3] = {34,36,44};

  if (ProblemType == 50) {

    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      pos[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	  pos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];

	  /* Find the super source closest to this position.  Then go
	     back up the tree, looking for the siblings of the
	     ancestors.  These leafs will be the sources that
	     contribute to the position. */


	  if (k==dbg[2] && j==dbg[1] && i==dbg[0])
	    DEBUG = TRUE;
	  else
	    DEBUG = FALSE;

	  FindSuperSourceByPosition(pos, &Leaf, DEBUG);

	  // Include deepest leaf first
	  radius2 = 0.0;
	  for (dim = 0; dim < MAX_DIMENSION; dim++) {
	    dx = Leaf->Position[dim] - pos[dim];
	    radius2 += dx*dx;
	  }
	  H2Luminosity = Leaf->LWLuminosity / radius2;

	  if (k==dbg[2] && j==dbg[1] && i==dbg[0])
	    printf("Leaf %d: L = %g, radius = %g\n", 
		   Leaf->LeafID, Leaf->LWLuminosity, sqrt(radius2));


	  // Now go up the tree, finding siblings of ancestors
	  LastLeaf = Leaf;
	  Leaf = Leaf->ParentSource;
	  while (Leaf != NULL) {
	    ci = -1;
	    if (Leaf->ChildSource[0] == LastLeaf &&
		Leaf->ChildSource[1] != NULL)
	      ci = 1;
	    else if (Leaf->ChildSource[1] == LastLeaf &&
		     Leaf->ChildSource[0] != NULL)
	      ci = 0;

	    if (ci >= 0) {
	      radius2 = 0.0;
	      for (dim = 0; dim < MAX_DIMENSION; dim++) {
		dx = Leaf->ChildSource[ci]->Position[dim] - pos[dim];
		radius2 += dx*dx;
	      }
	      if (k==dbg[2] && j==dbg[1] && i==dbg[0])
		printf("Leaf %d: ci=%d (leaf %d), L = %g %g, radius = %g\n", 
		       Leaf->LeafID, ci, Leaf->ChildSource[ci]->LeafID,
		       Leaf->ChildSource[ci]->LWLuminosity, H2Luminosity,
		       sqrt(radius2));
	      H2Luminosity += Leaf->ChildSource[ci]->LWLuminosity / radius2;
	    } // ENDIF (c >= 0)

	    LastLeaf = Leaf;
	    Leaf = Leaf->ParentSource;

	  } // ENDWHILE Leaf

	  BaryonField[kdissH2INum][index] += H2Luminosity * factor;

	} // ENDFOR i
      } // ENDFOR j
    } // ENDFOR k

      

  } // ENDIF ProblemType == 50

  else {

    ENZO_FAIL("Not implemented for star particles yet.")

  } // ENDELSE ProblemType == 50

  return SUCCESS;

}
